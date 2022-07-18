clear all; clc; close all; tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters, calibration targets and other settings%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.96^(1/52);  %Weekly household discount factor
pid=7*0.005/18;     %Weekly probability of dying
pir=7*1/18-pid;     %Weekly probability of recovering
phii=0.8;           %Productivity of infected people

deltav=0/52;        %Weekly probability of discovering a vaccine
deltac=0/52;        %Weekly probability of discovering a treatment
kappa=0;            %Slope of pid-function in endog. pid scenario 
                    %(medical preparedness scenario)
tau = 2;                    
          
load('country_data')
load('US_params')

cbarUS = cbar;
worked = NaN(size(hours));

WAHvec = NaN(size(hours));
cbarvec = NaN(size(hours));
thetavec = NaN(size(hours));
Avec = NaN(size(hours));
bedsvec = NaN(size(hours));

US = find(strcmp(iso,'USA'));

pid = pid/etavec(US);

for i = 1:length(hours)
delete(strcat(char(iso(i)),'.mat'))
beds = BED(i)/BED(US)*0.00042;
n_target = hours(i)/52;
inc_target = GDPpc(i)/52;
cbar = PLI(i)*cbarUS;
wah = WAH(i);
eta = etavec(i);


%Calibation targets for shares of pis-terms in T-function in SIR model
pis3_shr_target=2/3;                   %share of T_0 jump due general infections
pis1_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to consumption-based infections
pis2_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to work-based infections
RplusD_target=0.60;                    %total share of people infected and then either recovered or dead after epidemic

pop_ini=1;          %Initial population
i_ini=0.001;        %Initial infected

HH=260;             %Number of periods to solve and simulate the model
GG=104;             %Number of periods to repay loan

%containment policy
muc=zeros(HH,1);    %exogenous path for muc over time.
Loan = zeros(HH,1); %if you want e.g. a containment policy of
                    
                    
use_parallel=1;     %when optimal policy is computed, use_parallel=1 uses 
                    %parallel computing to maximize PV utility using fmincon.
                    
err = 0;                

clear getU
%nonlinear solver and minimizer settings
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve
opts_fsolve_fmincon=optimoptions('fsolve','Display','off','TolFun',1e-9); %options for fsolve used opt. policy calcs. (fmincon)
if use_parallel==0
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',2000,'FiniteDifferenceStepSize',1e-2); %options for fmincon w/o parallel comp.
    options_fminsearch = optimset('Display','iter');
elseif use_parallel==1
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-6,'MaxFunctionEvaluations',25000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);
    opts_fminunc=optimoptions('fminunc','Display','iter','TolFun',1e-6,'MaxFunctionEvaluations',15000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steady State Calculations%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A=inc_target/n_target;  %calculate parameter A such that income is equal to
                        %desired target, c=inc_target=A*n_target
                        


theta = 4/((2*n_target-cbar/A)^2-(cbar/A)^2);%calculate disutility of labor parameter theta
                        %so that pre-infection steady state labor is equal to
                        %desired target 

%steady states
nrss=(cbar/A+sqrt((cbar/A)^2+4/theta))/2;           %labor recovered (same as post-infection steady state)
crss=A*nrss;                    %consumption recovered
urss=log(crss-cbar)-theta/2*nrss^2;  %utility recovered
Urss=1/(1-betta)*urss;          %PV utility recovered
UrssConsUnits=Urss*crss;        %PV utility in cons. units (Urss*Marg.Util.Cons); value of life
niss=(cbar/(A*phii)+sqrt((cbar/(A*phii))^2+4/theta))/2;%(1/theta)^(1/2);           %labor infected
ciss=phii*A*niss;               %consumption infected
uiss=log(ciss-cbar)-theta/2*niss^2;  %utility infected
Uiss=1/(1-(1-deltac)*betta*(1-pir-pid*eta))*(uiss...
+(1-deltac)*betta*pir*Urss+deltac*betta*Urss);  %PV utility infected

if ciss<cbar
    err = 1;
end


if err == 0



%initial guess of vectors of ns, ni and nr to solve nonlinear
%equilibrium equations
n_vec_guess=nrss*ones(3*HH,1); %guess of vectors for ns,ni,nr



%Given either optimal path for muc or exogenous path for muc,
%solve nonlinear equilibrium model equations (i.e. adjust guesses ns,nr,ni)
[n_vec,fval,exitflag]=fsolve(@get_err,n_vec_guess,opts_fsolve,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid*eta,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa,cbar,Loan,wah,beds,tau,eta);
if exitflag~=1
    error('Fsolve could not solve the model');
end

if ~isnan(n_vec)
    
%get allocations given either exogenous or optimal path for muc at ns,ni,nr
%solution
[err,I,S,R,D,T,Pop,cs,ns,Us,RnotSIRmacro,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U,pid_endo] = get_err(n_vec,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid*eta,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa,cbar,Loan,wah,beds,tau,eta);
disp(['Max. abs. error in equilib. equations:',num2str(max(abs(err)))]);
disp(' ');


%plotting
ia=2;ib=2;fsize=12;
horz=HH;
time=0:1:horz-1;
%output some data used in paper and robustness table
aggCons_trough_percent=min((100*(aggC-crss)/crss));
aggCons_avg_first_year_percent=mean((100*(aggC(1:52)-crss)/crss));
terminal_one_minus_susceptibles_percent=100*(1-S(end));
peak_infection_percent=max(100*I);
terminal_death_share_percent=100*D(end);
terminal_number_deaths_US_millions=terminal_death_share_percent/100*330;
PVLoan = 0;
for j = 1:GG
    if Loan(j)>0
        PVLoan = PVLoan+1/rint^(j-1)*Loan(j);
    end
end
save(char(iso(i)))
worked(i) = 1;
WAHvec(i) = wah;
cbarvec(i) = cbar;
thetavec(i) = theta;
Avec(i) = A;
bedsvec(i) = beds;
end
end
end
worked(US) = 1;
save worked worked

save params WAHvec cbarvec thetavec Avec bedsvec

toc;
