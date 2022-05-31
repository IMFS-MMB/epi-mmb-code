clear 

setpath_dynare

global M_ oo_ changes_names changes_vec

modnam = 'two_sector_vcu';
modnamstar= 'two_sector_vcu_binding';

modnam_one_sector = 'one_sector_vcu';
modnamstar_one_sector= 'one_sector_vcu_binding';

capital_irreversibility_switch = 1;

constraint = 'in<-in_ss'; 
constraint_relax = 'lambdai<-0.0000001';

maxiter = 20;

n_periods_average = 6;

load_from_disk = 1;

if ~load_from_disk
    
    trans_prob0_vec = [0.011:0.00125:0.05725];
    for ibase = 1:length(trans_prob0_vec)

        eval(['dynare ',modnam,' noclearall'])

        paramfile_two_sector_vcu;
        
        
        
        % length of simulation
    months = 24;
    months_total = months;
    
    nperiods = months*30;             % duration until arrival of vaccine (measured in days)
    nperiods_total = months_total*30; % total length of simulation
    nperiods_plot = 360;
    
    % obtain contact matrices run US contacts
    UScontacts
    Jgroups = 4;
    share_working_while_infected = 0.4;
    
    if Jgroups == 1
        contacts_ = contacts1_;
        Ninit = Ninit_agg_1;
    elseif Jgroups == 3
        contacts_ = contacts3_;
        Ninit = Ninit_agg_3;
    elseif Jgroups == 4
        contacts_ = contacts4_;
        Ninit = Ninit_agg_4;
    elseif Jgroups == 5
        contacts_ = contacts5_;
        Ninit = Ninit_agg_5;
    end
    
    % epidemiological parameters
    gammap = 0.2*ones(Jgroups,1);    % transition rate from I to R
    thetap = 0.1*ones(Jgroups,1);    % transition rate from R to C and D
    deltap = 0.01*ones(Jgroups,1);   % death rate
    omegap = 0.0*deltap;             % slope of death rate
    trans_prob0 = trans_prob0_vec(ibase);             % initial transmission rate
    %trans_prob0 = 0.066;            % initial transmission rate
    irate = 0.00001;
    
    % contact restrictions (no restrictions if set to 1)
    res_.school = ones(Jgroups,nperiods_total);
    res_.work   = ones(Jgroups,nperiods_total);
    res_.other  = ones(Jgroups,nperiods_total);
    res_.home   = ones(Jgroups,nperiods_total);
    
    % transmission probability
    trans_prob  = trans_prob0*ones(1,nperiods_total);
    
    % basic contact rate
    betap0 = trans_prob0*(contacts_.home+contacts_.school+contacts_.work+contacts_.other);
    
    %% compute path of contact rate
    for t = 1:1:nperiods
        betamat(:,:,t) =  trans_prob(t)...
            *(contacts_.home.*(res_.home(:,t)*res_.home(:,t)')...
            + contacts_.school.*(res_.school(:,t)*res_.school(:,t)')...
            + contacts_.work.*(res_.work(:,t)*res_.work(:,t)')...
            + contacts_.other.*(res_.other(:,t)*res_.other(:,t)'));
        
        % for R0 computation check: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(09)62126-7/fulltext
        % should be the spectral radius of betap/gammap or max(abs(eig(betap(:,:,1)/gammap(1))))
        R0vec(t) = max(abs(eig(betamat(:,:,t)/gammap(1))));
    end
    
    betap_base_vec(ibase) = R0vec(1); 
    
    %% solve epidemiological model
    [Smat,Imat,Rmat,Dmat,Cmat,Nmat]...
        = SIRmodel_Jgroups_func(betamat,nperiods,nperiods_total,irate,Ninit,gammap,thetap,deltap,omegap);
    
    
    
    
    
    N1 = Ninit(2);
    N2 = Ninit(3);
    
    % take off infected, dead, and people under lockdown from the labor force
    l1_path_days = - (Ninit(2)-Nmat(2,:))...
        - Imat(2,:)*(1-share_working_while_infected); %...
    %  - L_1_vec.*(N_1_vec- I_1_vec*(1-share_working_while_infected));
    l2_path_days = - (Ninit(3)-Nmat(3,:))...
        - Imat(3,:)*(1-share_working_while_infected); %...
    %  - L_2_vec.*(N_2_vec- I_2_vec*(1-share_working_while_infected));
    pop_path_days = sum(-(Ninit*ones(1,nperiods_total)-Nmat));
    
    % temporal aggregation of SIR model from days to months
    l1_path = mean(reshape(l1_path_days,30,months_total))';
    l2_path = mean(reshape(l2_path_days,30,months_total))';
    pop_path = mean(reshape(pop_path_days,30,months_total))';
    
    
    l1_path_lag = [0; l1_path(1:end-1)];
    l2_path_lag = [0; l2_path(1:end-1)];
    pop_path_lag = [0; pop_path(1:end-1)];
    
    % define matrix needed to calculate the sizes of the shocks
    amat = diag(rhop_l*ones(months,1));
    
    % translate the path for a into the path for the shocks m
    m_path = l1_path-amat*l1_path_lag;
    n_path = l2_path-amat*l2_path_lag;
    
    
    % define matrix needed to calculate the sizes of the shocks
    amat = diag(rhop_l*ones(months,1));
    
    % translate the path for a into the path for the shocks m
    m_path = l1_path-amat*l1_path_lag;
    n_path = l2_path-amat*l2_path_lag;
    pop_shock_path = pop_path - amat*pop_path_lag;
    
    % put the shocks on the vector of predetermined conditions
    m_pos_vec = zeros(months,1);
    n_pos_vec = zeros(months,1);
    pop_shock_pos_vec = zeros(months,1);
    
    eval(['dynare ',modnam,' noclearall'])
    
    endo_names = cellstr(M_.endo_names);
    for m_indx = 1:months
        this_m = ['m',num2str(m_indx)];
        this_n = ['n',num2str(m_indx)];
        this_pop_shock = ['pop_shock',num2str(m_indx)];
        m_pos_vec(m_indx)=find(strcmp(this_m, endo_names));
        n_pos_vec(m_indx)=find(strcmp(this_n, endo_names));
        pop_shock_pos_vec(m_indx) = find(strcmp(this_pop_shock, endo_names));
    end
    
    init = zeros(M_.endo_nbr,1);
    init(m_pos_vec)=m_path;
    init(n_pos_vec)=n_path;
    init(pop_shock_pos_vec) = pop_shock_path;

%         % call SIR model
% 
%         months = 6;
%         months_lockdown = 0;
% 
%         effectiveness_share = 0.8;  % irrelevant without a lockdown
%         share1_lockdown = 0;
%         share2_lockdown = 0;
%         share3_lockdown = 0;   
% 
%         share_working_while_infected = 0.4;
% 
%         N1 = l1_ss;
%         N2 = l2_ss;
% 
%         % run as daily model
%         betap_base = betap_base_vec(ibase);
% 
%         gammap1 = 1/20;    % recovery rate
%         gammap2 = gammap1;
%         gammap3 = gammap1;
% 
%         omegap1 = 0;       % death rate -- abstract from the effects of deaths for now
%         omegap2 = omegap1;
%         omegap3 = omegap1;
% 
% 
%         [S1,R1,I1,S2,R2,I2,S3,I3,R3]=SIRmodel_3groups_func(months,...
%                 months_lockdown,effectiveness_share,share1_lockdown,share2_lockdown,share3_lockdown,...
%                 N1,N2,betap_base,gammap1,gammap2,gammap3,omegap1,omegap2,omegap3);
% 
%         % temporal aggregation of SIR model from days to month
%         I1_month = mean(reshape(I1',30,months));
%         I2_month = mean(reshape(I2',30,months));
%         I3_month = mean(reshape(I3',30,months));
% 
%         
%         % path for the exogenous variable
%         paramfile_two_sector_vcu
% 
%         % take off the infected folks from the labor force
%         l1_path =  - I1_month'*(1-share_working_while_infected);
%         l2_path =  - I2_month'*(1-share_working_while_infected);
% 
%         l1_path_lag = [0; l1_path(1:end-1)];
%         l2_path_lag = [0; l2_path(1:end-1)];
% 
%         % define matrix needed to calculate the sizes of the shocks
%         amat = diag(rhop_l*ones(months,1));
% 
%         % translate the path for a into the path for the shocks m
%         m_path = l1_path-amat*l1_path_lag;
%         n_path = l2_path-amat*l2_path_lag;
% 
% 
%         % put the shocks on the vector of predetermined conditions
%         m_pos_vec = zeros(months,1);
%         n_pos_vec = zeros(months,1);
% 
%         endo_names = cellstr(M_.endo_names);
%         for m_indx = 1:months
%             this_m = ['m',num2str(m_indx)];
%             this_n = ['n',num2str(m_indx)];
%             m_pos_vec(m_indx)=find(strcmp(this_m, endo_names));
%             n_pos_vec(m_indx)=find(strcmp(this_n, endo_names));
%         end
%         init = [];
%         init = zeros(M_.endo_nbr,1);
%         init(m_pos_vec)=m_path;
%         init(n_pos_vec)=n_path;

        % now get the response to the initial conditions
        shockssequence = 0;  
        irfshock = '';
        nperiods = months;

        changes_names = char('rhop');
        changes_vec =       [-3/2];

        if capital_irreversibility_switch

        [zdatalinear, zdatapiecewise, zdatass, oobase_ Mbase_] = ...
        solve_one_constraint(modnam,modnamstar,...
        constraint, constraint_relax,...
        shockssequence,irfshock,nperiods,maxiter,init)

        zdatalinear = zdatapiecewise;

        else

        [zdatalinear, zdatass, oobase_, Mbase_ ] = ...
            solve_no_constraint(modnam,...
            shockssequence,irfshock,nperiods,init);
        end

        %% now do it again but for a one sector model
        N1 = 0;
        N2 = l1_ss+l2_ss;

        

        l2_path =  l1_path+l2_path;
        l2_path_lag = [0; l2_path(1:end-1)];


        n_path = l2_path-amat*l2_path_lag;

        eval(['dynare ',modnam_one_sector,' noclearall'])

        endo_names = cellstr(M_.endo_names);
        for m_indx = 1:months

            this_n = ['n',num2str(m_indx)];
            n_pos_vec(m_indx)=find(strcmp(this_n, endo_names));
            pop_shock_pos_vec(m_indx) = find(strcmp(this_pop_shock, endo_names));
        end

        init = zeros(M_.endo_nbr,1);
        init(n_pos_vec)=n_path;
        init(pop_shock_pos_vec) = pop_shock_path;

        if capital_irreversibility_switch

        [zdatalinear2, zdatapiecewise2, zdatass2, oobase2_ Mbase2_] = ...
        solve_one_constraint(modnam_one_sector,modnamstar_one_sector,...
        constraint, constraint_relax,...
        shockssequence,irfshock,nperiods,maxiter,init)

        zdatalinear2 = zdatapiecewise2;

        else

        [zdatalinear2, zdatass2, oobase2_, Mbase2_ ] = ...
            solve_no_constraint(modnam_one_sector,...
            shockssequence,irfshock,nperiods,init);
        end

        prefix1 = 's1_';
        prefix2 = 's2_';
        for this_var = 1:Mbase_.endo_nbr
            eval([prefix1,deblank(Mbase_.endo_names(this_var,:)),'_irf=zdatalinear(:,this_var);'])
            eval([prefix1,deblank(Mbase_.endo_names(this_var,:)),'_ss=zdatass(this_var);'])   
        end 

        for this_var = 1:Mbase2_.endo_nbr
            eval([prefix2,deblank(Mbase2_.endo_names(this_var,:)),'_irf=zdatalinear2(:,this_var);'])
            eval([prefix2,deblank(Mbase2_.endo_names(this_var,:)),'_ss=zdatass2(this_var);'])   
        end
        s1_cumy(ibase,1) = mean(s1_y_irf(1:n_periods_average))/s1_y_ss;
        s2_cumy(ibase,1) = mean(s2_y_irf(1:n_periods_average))/s2_y_ss;   
        %s1_maxI(ibase,1) = max(I1_month+I2_month+I3_month);
    end   
    

    
    

%save sensitivity_results_six_months s1_cumy s2_cumy betap_base_vec

save sensitivity_results_six_months_chip_8_10 betap_base_vec s1_cumy s2_cumy

else

%%

figure; subplot(2,1,1)
load sensitivity_results_six_months

plot(betap_base_vec,-100*s1_cumy,'xb')
hold on;plot(betap_base_vec,-100*s2_cumy,'or')
xlabel('R0')
ylabel('Average output loss over 6 months')
title('Initial R0 and Cumulative Output Loss')

legend('Two-Sector Model','One-Sector Model')

subplot(2,1,2)

load sensitivity_results_six_months_chip_8_10

plot(betap_base_vec,-100*s1_cumy,'xb')
hold on;plot(betap_base_vec,-100*s2_cumy,'or')
xlabel('R0')
ylabel('Average output loss over 6 months')
title('Initial R0 and Cumulative Output Loss: A Lower Minimum Scale For Sector 1')

printpref
print -depsc2 figure_beta_sensitivity.eps

end
