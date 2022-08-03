%Model code for 'Expectations, Infections, and Economic Activity' by 
%M. Eichenbaum, M. Godinho de Matos, F. Lima, S. Rebelo and M. Trabandt

%April 2022. mathias.trabandt@gmail.com

%Baseline estimated model

close all;  clc; %clear all;

rng('default');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set switches for desired calculations
mpar.do_estimation=0;

mpar.start_mode_finder=0;
mpar.load_last_estimation_step=0;
mpar.est_maxiter=1000;
mpar.do_checkplots=0;
mpar.grdpts=5;%checkplot must be odd number of gridpoints
mpar.use_parallel=1;

mpar.do_mcmc=0;
mpar.do_load_mcmc_results=0;

mpar.nbsim=1; %number of grid points (assets) above and below median assets grid point
%for epi simulation (subset of mpar.nb);
%centered around median assets below
%using a subset of grid points around median assets for the simulation speeds up computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mpar.infect_surprise=1; %if =1, surprise second wave
mpar.infect_surprise_week=17; %week when first wave ends; next wave starts in mpar.infect_surprise_week+1 (unexpectedly)


mpar.timestart=5;  %start of simulation (week since March 1)

fminbnd_options=optimset('Display','off','MaxFunEvals',100,'TolX',1e-12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Infections, containment and cfr's for simulation%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get weekly new deaths and calculate implied infections
get_weeklynewdeaths;
Ivec_RAW=newdeaths_weekly.NewDeaths(3:end)/0.011/10000000; %weekly infections as share of initial population
Ivec=[Ivec_RAW',0];       %add zero as last observation to solve the model with perfect foresight and given the i=0 value functions as a start for backward induction
Ivec=Ivec(mpar.timestart:end); %start simulation at week mpar.timestart

%get containment measure
get_contain;
muvec=Ivec*0;
muvec(1:end)=contain_weekly.contain(mpar.timestart:end-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define numerical parameters
mpar.nb   = 150;        %number of grid points, assets
mpar.minb = 50;         %lowest point on the asset grid; chosen such that log grid for b is close to 75k
mpar.maxb = 150000;     %highest point on the asset grid

mpar.crit_run1=1e-7;    %precision of value function iteration
mpar.crit = 1e-3;       %precision up to which to solve the value function with broyden
mpar.broydencritx=1e-6; %min dx in broyden
mpar.broydeniter=200;   %max iter broyden

mpar.do_sol_sim=1;      %solve and simulate model or load from disk otherwise
                        %be careful: solving the model can take substantial
                        %time (value function iteration)

mpar.loadVini=1;        %use old solution as initial guess when solving the model

%Parameters for SIR model calibration
mpar.opts_fsolve=optimoptions('fsolve','Display','off','TolFun',1e-7,'MaxFunctionEvaluations',100000,'MaxIterations',4000); %options for fsolve
mpar.HH=200;            %number of periods in SIR model simulation

mpar.scale1=1000;       %scale pi1 for numerical solver

%get initial infections in the data
par.i_ini=Ivec(mpar.timestart+1);        %Initial seed of total infections for basic SIR model; Start in week 5 of march 2020 in line with treatment of infections above
if par.i_ini==0,disp('Initial infections data are zero. Setting it to 0.1%'); par.i_ini=0.001;end

%% Parameters

%% Define economic parameters
par.betta=0.97^(1/52);  %discount factor
par.deltao=1/(13*52);   %prob of dying naturally old (average life expectancy old...)
par.deltay=1/(51*52);   %prob of dying naturally young
par.nu = 1/(28*52);     %prob. of becoming old/aging

par.days=14;
par.pidy=7*0.001/par.days;        %Weekly probability of dying, young
par.piry=7*1/par.days-par.pidy;   %Weekly probability of recovering, young
par.pido=7*0.035/par.days;        %Weekly probability of dying, old
par.piro=7*1/par.days-par.pido;   %Weekly probability of recovering, old

par.pi1=0.588/1000;  %Transmission-function: consumption-based infections 
par.pi2=1.0418;      %Transmission-function: general infections
par.Rnottarget=2.5;  %target value for R0. pi2 adjusted in simulate_model to set R0

%par.Rnottarget=2.5;
%par.pi1=0.588/1000
%par.pi2=par.Rnottarget*(par.piry*par.sy_ss+par.piro*(1-par.sy_ss)+par.pidy*par.sy_ss+par.pido*(1-par.sy_ss))-(par.pi1*par.cry_ss*par.sy_ss+par.pi1*par.cro_ss*(1-par.sy_ss));
%par.pi1_shr_target=(par.pi1*(par.cry_ss*par.sy_ss+par.cro_ss*(1-par.sy_ss)))/(par.pi1*(par.cry_ss*par.sy_ss+par.cro_ss*(1-par.sy_ss))+par.pi2) %take weighted avg of old and young consumption


%Estimated parameters for epi and containment simulation
par.muvec_scale= 0.128350297916290;          %scaling factor containment measure
par.gyoung=0.074207800064035;                %const. gain param
par.gold=0.131968838244525;                  %const. gain param
par.pidy_ini=0.066599825791575*0.9;          %initial value const. gain learning, young
par.pido_ini=0.371397964352873*0.9;          %initial value const. gain learning, old
       

%further model parameters
par.w=19000/52; %net of taxes income per year
par.sy_ss=0.7;  %share of young population

par.alf=2;      %EZ risk aversion (Albuquerque, Eichenbaum, Luo, Rebelo, benchmark model)
par.rho=1/1.5;  %Inverse IES (Albuquerque, Eichenbaum, Luo, Rebelo, benchmark model)

par.kap=(1-par.rho)/(1-par.alf); %aux. parameter
par.gam=1/(1-par.rho); %aux. parameter

par.r=1.01^(1/52)-1; %real interest rate

par.theta=4;    %slope in bequest utility (in terms of assets)
par.xi=120;     %intercept in bequest utility;
par.z=2;        %constant in utility function (outside EZ aggregator)

%targets:
%vol: 890.000
%ratio of cons young and old 1.18 at median assets
%weigthed avg savings rate is 6.7% at median assets

par.median_assets=75000; %median assets (used to find appropriate idx in grid.b below)

%% Grid
grid.b = exp(linspace(log(mpar.minb),log(mpar.maxb),mpar.nb));

[~,par.grid_b_median_assets_idx]=min(abs(grid.b-par.median_assets)); %get index for median assets
mpar.grid_sim_idx=[par.grid_b_median_assets_idx-mpar.nbsim:1:par.grid_b_median_assets_idx+mpar.nbsim]; %indexes for simulation subset of grid points

%% Define utility function
util  = @(c) (1-par.betta)*(c)^(1-par.rho);

%% Calculate Income Levels
[meshes.b]= ndgrid(grid.b);
Y=par.w+(1+par.r)*meshes.b;

%% Initialize array to store value functions by types
U=[];

%%Definitions of flags that will be used in code later on
%flag=1 - ro
%flag=2 - ry
%flag=3 - io
%flag=4 - iy
%flag=5 - so
%flag=6 - sy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SOLVE VALUE FUNCTION FOR I=0 for all t     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('Solving Value Functions');
disp('-----------------------');
if mpar.do_sol_sim==1
    tstart=tic;
    V_ini=2000+zeros(mpar.nb,6);
    for flag=1:1:6
        if mpar.loadVini==1
            load V_ini;
            V=V_ini(:,flag);
        else
            V=V_ini(:,flag);
        end
        I=0;mu=0;solflag=1;
        dist=9999;count=1;tic;
        dist_V = @(V)(V-VFI_update_spline(V,Y,util,par,mpar,grid,U,flag,solflag,I,mu,[],[],[],[],[],[],fminbnd_options));
        while dist(count)>mpar.crit_run1
            count = count+1;
            DV    = dist_V(V);
            dist(count) = max(abs(DV(:))); % Calculate distance between old guess and update
            V     = V-DV;
            %dist(count)
        end
        Vsol=V;
        %[Vsol,~,count,dist] = broyden(dist_V,V,mpar.crit,mpar.broydencritx,mpar.broydeniter);
        [~,bprimesol]=VFI_update_spline(Vsol,Y,util,par,mpar,grid,U,flag,solflag,I,mu,[],[],[],[],[],[],fminbnd_options);
        csol=par.w+(1+par.r)*grid.b(:)-bprimesol;
        disp(['Flag ',num2str(flag),': Number of evaluations: ',num2str(count-1),'. Time: ',num2str(toc,2), ' sec.']);
        
        if     flag==1, U.ro=Vsol; bprime.ro=bprimesol;c.ro=csol;distF.ro=dist;
        elseif flag==2, U.ry=Vsol; bprime.ry=bprimesol;c.ry=csol;distF.ry=dist;
        elseif flag==3, U.io=Vsol; bprime.io=bprimesol;c.io=csol;distF.io=dist;
        elseif flag==4, U.iy=Vsol; bprime.iy=bprimesol;c.iy=csol;distF.iy=dist;
        elseif flag==5, U.so=Vsol; bprime.so=bprimesol;c.so=csol;distF.so=dist;
        elseif flag==6, U.sy=Vsol; bprime.sy=bprimesol;c.sy=csol;distF.sy=dist;
        end
    end
    disp(' ');
    disp(['Total time in minutes: ',num2str((toc(tstart))/60,3)]);
    
    V_ini=[U.ro U.ry U.io U.iy U.so U.sy];
    save V_ini V_ini;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Command Window Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp pi1 and pi2
    disp([' pi1 = ',num2str(par.pi1),', pi2 = ',num2str(par.pi2)]);
    par.cro_ss=c.ro(par.grid_b_median_assets_idx); par.cry_ss=c.ry(par.grid_b_median_assets_idx); %steady state consumption old (zero assets)
    
    %Calc. value of life and print out results to report in paper
    disp(' ');
    muc_ro=(U.ro(1)-par.z)^par.rho*(1-par.betta)*c.ro(1)^(-par.rho);
    VoL_old_poorest_millions=U.ro(1)/muc_ro/1000000;
    muc_ro=(U.ro(end)-par.z)^par.rho*(1-par.betta)*c.ro(end)^(-par.rho);
    VoL_old_richest_millions=U.ro(end)/muc_ro/1000000;
    
    muc_ry=(U.ry(1)-par.z)^par.rho*(1-par.betta)*c.ry(1)^(-par.rho);
    VoL_young_poorest_millions=U.ry(1)/muc_ry/1000000;
    muc_ry=(U.ry(end)-par.z)^par.rho*(1-par.betta)*c.ry(end)^(-par.rho);
    VoL_young_richest_millions=U.ry(end)/muc_ry/1000000;
    
    muc_ro=(U.ro(par.grid_b_median_assets_idx)-par.z)^par.rho*(1-par.betta)*c.ro(par.grid_b_median_assets_idx)^(-par.rho);
    VoL_old_median_assets_millions=U.ro(par.grid_b_median_assets_idx)/muc_ro/1000000;
    muc_ry=(U.ry(par.grid_b_median_assets_idx)-par.z)^par.rho*(1-par.betta)*c.ry(par.grid_b_median_assets_idx)^(-par.rho);
    VoL_young_median_assets_millions=U.ry(par.grid_b_median_assets_idx)/muc_ry/1000000;
    
    cro_median_assets=c.ro(par.grid_b_median_assets_idx);
    cry_median_assets=c.ry(par.grid_b_median_assets_idx);
    
    cro_min=min(c.ro);cry_min=min(c.ry);
    cry_cro_ratio=cry_median_assets/cro_median_assets;
    savings_rate_ro=(par.w+par.r*grid.b(par.grid_b_median_assets_idx)-c.ro(par.grid_b_median_assets_idx))/(par.w+par.r*grid.b(par.grid_b_median_assets_idx));
    savings_rate_ry=(par.w+par.r*grid.b(par.grid_b_median_assets_idx)-c.ry(par.grid_b_median_assets_idx))/(par.w+par.r*grid.b(par.grid_b_median_assets_idx));
    avg_savings_rate=par.sy_ss*savings_rate_ry+(1-par.sy_ss)*savings_rate_ro;
    
    disp(['Median assets on grid: ' num2str(grid.b(par.grid_b_median_assets_idx))]);
    disp(['VoL old, med. assets, mill = ',num2str(VoL_old_median_assets_millions),', VoL young, med. assets, mill =  ',num2str(VoL_young_median_assets_millions)]);
    disp(['cro, med. assets = ',num2str(cro_median_assets),', cry, med. assets =  ',num2str(cry_median_assets)]);
    disp(['savings_rate_ro, med. assets = ',num2str(savings_rate_ro),', savings_rate_ry, med. assets =  ',num2str(savings_rate_ry)]);
    disp(['avg_savings_rate, med. assets = ',num2str(avg_savings_rate),', cry_cro_ratio, med. assets =  ',num2str(cry_cro_ratio)]);
    disp(['min cro = ',num2str(cro_min),', min cry =  ',num2str(cry_min)]);
    disp(' ');
    
    par.pi1_shr_ini=(par.pi1*(par.cry_ss*par.sy_ss+par.cro_ss*(1-par.sy_ss)))/(par.pi1*(par.cry_ss*par.sy_ss+par.cro_ss*(1-par.sy_ss))+par.pi2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%SIMULATION OF EPIDEMIC                     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %above we have solved for the value functions for I=0, i.e. no infections.
    %Now we simulate a perfect foresight simulation starting with the terminal
    %period and recursing backbard in time. In each period, agents use the next
    %period value function form the previous period. The code also allows
    %to set a switch such that the second and third waves are unexpected
    %(i.e. no perfect foresight)
    
    %Calculate Sy, So etc which is needed later to calulate per capita cons
    %initial conditions
    Iy=Ivec'*par.sy_ss;Io=Ivec'*(1-par.sy_ss);
    
    Dy(1)=0;Do(1)=0;Ry(1)=0;Ro(1)=0;
    
    %iterate on SIR model equations
    for j=1:1:numel(Ivec)-1
        par.pidy=7*0.001/par.days;        %Weekly probability of dying, young
        par.piry=7*1/par.days-par.pidy;   %Weekly probability of recovering, young
        par.pido=7*0.035/par.days;        %Weekly probability of dying, old
        par.piro=7*1/par.days-par.pido;   %Weekly probability of recovering, old
        
        Ry(j+1,1)=Ry(j)*(1-par.deltay-par.nu)+par.piry*Iy(j);
        Ro(j+1,1)=Ro(j)*(1-par.deltao)+par.piro*Io(j)+par.nu*Ry(j);
        
        Dy(j+1,1)=Dy(j)+par.pidy*Iy(j);
        Do(j+1,1)=Do(j)+par.pido*Io(j);
    end
    Sy=par.sy_ss-Ry-Dy-Iy;
    So=(1-par.sy_ss)-Ro-Do-Io;
    
    
    %initialize value functions etc over time.
    Umat.ro=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));Umat.ry=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));
    Umat.io=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));Umat.iy=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));
    Umat.so=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));Umat.sy=zeros(numel(U.ro(mpar.grid_sim_idx)),numel(Ivec));
    
    Umat.ro(:,end)=U.ro(mpar.grid_sim_idx);Umat.ry(:,end)=U.ry(mpar.grid_sim_idx);
    Umat.io(:,end)=U.io(mpar.grid_sim_idx);Umat.iy(:,end)=U.iy(mpar.grid_sim_idx);
    Umat.so(:,end)=U.so(mpar.grid_sim_idx);Umat.sy(:,end)=U.sy(mpar.grid_sim_idx);
    
    bprimemat.ro=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));bprimemat.ry=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));
    bprimemat.io=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));bprimemat.iy=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));
    bprimemat.so=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));bprimemat.sy=zeros(numel(bprime.ro(mpar.grid_sim_idx)),numel(Ivec));
    
    bprimemat.ro(:,end)=bprime.ro(mpar.grid_sim_idx);bprimemat.ry(:,end)=bprime.ry(mpar.grid_sim_idx);
    bprimemat.io(:,end)=bprime.io(mpar.grid_sim_idx);bprimemat.iy(:,end)=bprime.iy(mpar.grid_sim_idx);
    bprimemat.so(:,end)=bprime.so(mpar.grid_sim_idx);bprimemat.sy(:,end)=bprime.sy(mpar.grid_sim_idx);
    
    solflag=0; %simulate rather than solve the value function for i=0 for all t.
    
    %data consumption young and old (see xls file) (monthly; starts in
    %march 20; goes thru april 21)
    dat_young=100*[-0.127	-0.322	-0.211	-0.012	0.046	0.046	0.024	-0.084	-0.059	-0.127	-0.296	-0.212	-0.083	-0.072];
    dat_old=100*[-0.136	-0.419	-0.254	-0.074	-0.027	-0.040	-0.031	-0.119	-0.109	-0.184	-0.363	-0.254	-0.099	-0.085];
    
    %standard deviations consumption young and old (see xls file)
    dat_stderr_young=[0.004577447	0.005042758	0.004922426	0.004931662	0.005567094	0.005077663	0.005096216	0.005207857	0.005260694	0.005419005	0.00568363	0.005701155	0.005936621	0.005972153];
    dat_stderr_old  =[0.004376884	0.005031653	0.004743569	0.004813666	0.005235872	0.004921257	0.004925781	0.005106109	0.005132598	0.005329912	0.005757338	0.005721398	0.005909867	0.005981088];
    par.dat_stderr_young=dat_stderr_young;
    par.dat_stderr_old=dat_stderr_old;
    
    %precision matrix
    VVvec=[dat_stderr_young(1:end)';dat_stderr_old(1:end)'].^2;
    VVmat=diag(VVvec);
    invVVmat=inv(VVmat);
    logdetVVmat=log(det(VVmat));
    par.invVVmat=invVVmat;
    par.logdetVVmat=logdetVVmat;
    
    %model implied data
    dt_cons = datetime('01/Mar/2020') : datetime('15/May/2021');
    
    %estimated parameter names
    par.guess_name=[{'mu_scale   '};{'pid0_young'};{'pid0_old  '};{'g_young    '};{'g_old      '}]; %];
    
    %initial guess     
     guess=[ 0.1
             0.1
             0.1
             0.1
             0.1
             ];
        
        
%     guess=[    0.128350297916290
%                0.066599825791575*0.9
%                0.371397964352873*0.9
%                0.074207800064035
%                0.131968838244525
%             ];
        
       
    %these are lower and upper bounds for fmincon mode-finder; 
    %all priors are uniform -- these are the LB and UB for the priors too
    par.LB=[0;0;0;0;0]; 
    par.UB=[1;7/par.days;7/par.days;1;1]; 
    
    
       
    if mpar.do_estimation==1
        tic
        disp(' ');
        disp('Starting Estimation');
        disp('-------------------');
        
        if mpar.use_parallel==0
            %opts_fmincon=optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-7,'StepTolerance',1e-7,'TolFun',1e-7,'MaxFunctionEvaluations',2,'UseParallel',false,'FiniteDifferenceStepSize',1e-7,'Algorithm','sqp'); %options for fmincon w/o parallel comp.
            opts_fmincon=optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',mpar.est_maxiter,'UseParallel',false,'algorithm','sqp','ConstraintTolerance',1e-12,'StepTolerance',1e-12,'TolFun',1e-12,'FiniteDifferenceStepSize',1e-4); %options for fmincon w/o parallel comp.
        elseif mpar.use_parallel==1
            %opts_fmincon=optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-12,'StepTolerance',1e-12,'TolFun',1e-12,'MaxFunctionEvaluations',3000,'UseParallel',true,'FiniteDifferenceStepSize',1e-8,'Algorithm','sqp');
            opts_fmincon=optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',mpar.est_maxiter,'UseParallel',true,'algorithm','sqp','ConstraintTolerance',1e-12,'StepTolerance',1e-12,'TolFun',1e-12,'FiniteDifferenceStepSize',1e-4); %options for fmincon w/o parallel comp.
        end
        
        if mpar.load_last_estimation_step==1
            disp('Loading last mode finder step');
            load last_estimation_step
            guess=sol;
        end
        
        if mpar.start_mode_finder==1
               
            [sol,post_value,exitflag,OUTPUT,LAMBDA,GRAD,HESSIAN_FMINCON] = fmincon(@simulate_model,guess,[],[],[],[],par.LB,par.UB,[],opts_fmincon,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
            
            disp('Computing Hessian: this may take a while...');
            disp('-------------------------------------------');
            hessian_mat=reshape(hessian('simulate_model',sol,eps,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options),numel(sol),numel(sol));
            %hessian_mat=HESSIAN_FMINCON;
            
            invhessian=inv(hessian_mat);  %inverse hessian
            
            save last_estimation_step sol hessian_mat invhessian;%save solution for possible use in subsequent maximization
        else
            load last_estimation_step;
        end
        
        
        %get parameter at mode
        par.muvec_scale= sol(1); %scaling factor containment measure
        par.pidy_ini=sol(2);     %initial value const. gain learning, young
        par.pido_ini=sol(3);     %initial value const. gain learning, old
        par.gyoung=sol(4);       %const. gain param, young
        par.gold=sol(5);         %const. gain param, old       
      
        
        disp(' ');
        disp('Simulating with estimated parameters');
        disp('------------------------------------');
        [neglogpost,cmat,comat,cymat,cons_monthly,cfryoung_constgain,cfrold_constgain]=simulate_model(sol,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
        
        %parameter posterior mode and posterior value
        para_mode=sol;
        logpost_mode=-neglogpost;
        
        %calculate laplace
        Laplace = 1/2*numel(sol)*log(2*pi) + (logpost_mode)- 1/2*log(det(invhessian))
        
        disp(' ');
        disp('Estimated Parameter(s)');
        disp('----------------------');
        for iii=1:1:numel(sol)
            disp([par.guess_name(iii),sol(iii)])
        end
        
        %standard deviations and parameter correlations
        [para_standard_deviations,para_correlation_matrix]=Cov2Corr(invhessian);
        para_standard_deviations=para_standard_deviations'
        
        % Posterior correlations
        ExtremeCorrBound=0.15;
        if  ~isnan(ExtremeCorrBound)
            %[Stds,CorrMat]=Cov2Corr(InvHessian); % Computing correlations and standard deviations from InvHessian.
            tril_para_correlation_matrix=tril(para_correlation_matrix,-1);
            [RowIndex,ColIndex]=find(abs(tril_para_correlation_matrix)>ExtremeCorrBound);
            ExtremeCorrParams=cell(length(RowIndex),3);
            for i=1:length(RowIndex)
                ExtremeCorrParams{i,1}=char(par.guess_name(RowIndex(i)));
                ExtremeCorrParams{i,2}=char(par.guess_name(ColIndex(i)));
                ExtremeCorrParams{i,3}=tril_para_correlation_matrix(RowIndex(i),ColIndex(i));
            end
        end
        disp(' ');
        disp(['Correlations of Parameters (at Posterior Mode) > abs(',num2str(ExtremeCorrBound),')']);
        disp(ExtremeCorrParams)
        
        %prior posterior plots
        figure('name','Prior-Posterior (Laplace)');
        ia=2;ib=3;
        trunc = 1e-12; steps = 2^15-1;       
                   
        for zz=1:1:numel(sol)            
                a = sol(zz);b = para_standard_deviations(zz);infbound = norminv(trunc,a,b);
                supbound = norminv(1-trunc,a,b); stepsize = (supbound-infbound)/steps;
                gg = infbound:stepsize:supbound;
                post_dens = normpdf(gg,a,b);
            
            if zz==1 %muscale -- uniform prior
                gg_prior=par.LB(zz):stepsize:par.UB(zz);
                prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
                
            elseif zz==2 %pidy_ini young -- uniform prior
                gg_prior=par.LB(zz):stepsize:par.UB(zz);
                prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
                                       
            elseif zz==3 %pido_xini old -- uniform prior
                gg_prior=par.LB(zz):stepsize:par.UB(zz);
                prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
                
            elseif zz==4 %gyoung  -- uniform prior
                gg_prior=par.LB(zz):stepsize:par.UB(zz);
                prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
                
            elseif zz==5 %gold  -- uniform prior
                gg_prior=par.LB(zz):stepsize:par.UB(zz);
                prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
            end
            subplot(ia,ib,zz)
            plot(gg_prior,prior_dens,'k--','linewidth',1.5); hold on
            plot(gg,post_dens,'r-','linewidth',2); hold on; vline(sol(zz));
            axis tight
            title(char(par.guess_name(zz)),'Interpreter','none');
        end
        suptitle('Priors vs. Posteriors (Laplace Approx.)');
        legend1=legend('Prior','Posterior (Laplace)');
        set(legend1,...
            'Position',[0.418594533212181 0.906613971512973 0.174041297935103 0.0217391304347826],...
            'Orientation','horizontal');
        legend boxoff;orient landscape;
        print -dpdf -fillpage prior_posterior_laplace
        
        %check plots
        if mpar.do_checkplots
            figure('name','Check Plots');
            disp('Computing checkplots. This may take a while...');
            ia=2;ib=3;
            for zz=1:1:numel(sol)
                check_grid=linspace(sol(zz)-2*para_standard_deviations(zz),sol(zz)+2*para_standard_deviations(zz),mpar.grdpts);
                for ww=1:1:numel(check_grid)
                    check_sol=sol;check_sol(zz)=check_grid(ww);
                    [neglogpost_check]=simulate_model(check_sol,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
                    logpost_check(ww)=-neglogpost_check;
                end
                subplot(ia,ib,zz)
                plot(check_grid,logpost_check,'k-','linewidth',2); hold on;vline(sol(zz));
                title(char(par.guess_name(zz)),'Interpreter','none');
            end
            suptitle('Log-Posterior Slices (Check plot)');
            orient landscape; print -dpdf -fillpage checkplots
        end
        
        disp(['Time in hrs: ',num2str(toc/3600)]);
        
        
         if mpar.do_mcmc==1
             disp('Starting MCMC. This may take a while....');
             mcmc;
         elseif mpar.do_load_mcmc_results==1
             disp('Loading MCMC results.');
             load mcmc_results;
         else
            disp('Taking posterior mode based on optimizer...');
            post_mode_mcmc=sol;
         end
        
        par.muvec_scale= post_mode_mcmc(1);             %scaling factor containment measure
        par.pidy_ini=post_mode_mcmc(2);                 %initial value const. gain learning, young
        par.pido_ini=post_mode_mcmc(3);                 %initial value const. gain learning, old
        par.gyoung=post_mode_mcmc(4);                   %const. gain param, young
        par.gold=post_mode_mcmc(5);                     %const. gain param, old
        
        disp('Simulating');
        [neglogpost,cmat,comat,cymat,cons_monthly,cfryoung_constgain,cfrold_constgain]=simulate_model(post_mode_mcmc,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);     
        
    else
        disp(' ');
        disp('Simulating with given parameters');
        disp('--------------------------------');
        tic
        %load all_results;
        %put final estimation results here (in terms of parameters)
        %par.guess_name=[{'mu_scale   '};{'pid_young'};{'pid_old  '};{'g_young    '};{'g_old      '}]; %];
        para=[ 
               0.128619166504464
               0.058564912204341
               0.311710237019398
               0.072570861873682
               0.132476252860019
             ];
        %profile on;
        tic;
        [neglogpost,cmat,comat,cymat,cons_monthly,cfryoung_constgain,cfrold_constgain]=simulate_model(para,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
        toc;
        %profile viewer;
    end
    
    %put data in time table format
    tty_data = timetable(cons_monthly.Cy.Time(1:end-1),dat_young','VariableNames',{'Consy_data'});
    tto_data = timetable(cons_monthly.Cy.Time(1:end-1),dat_old','VariableNames',{'Conso_data'});
    
    save all_results;
else
    load all_results;
end

%plot all results
%plot_results;
 






