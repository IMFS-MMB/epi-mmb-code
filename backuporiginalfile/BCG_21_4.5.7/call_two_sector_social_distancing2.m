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


slide_switch = 0;


eval(['dynare ',modnam,' noclearall'])

for case_switch = 1
    
    paramfile_two_sector_vcu;
      
    months = 36;
    months_lockdown = 3;
        
    nperiods_plot = 360;
    
    rhop0 = 0.1/30;
    
    %trans_prob0 = 0.011                % R0 = 0.9;
    trans_prob0 = 0.0255;             % initial transmission rate --- R0 = 2
    
    %trans_prob0  = 0.04;             % 3.1
    %trans_prob0 = 0.05725;           % initial transmission rate --- R0 = 4.5

       
    
    share_working_while_lockdown_1 = 0.15;
    share_working_while_lockdown_2 = 0.4;
    
    N1 = l1_ss;
    N2 = l2_ss;
    
    
    if case_switch ==1
        % first case --- let folks who can work from home
        effectiveness_share = 1;  % irrelevant without a lockdown
        share1_lockdown = 0.15;
        share2_lockdown = 0.40;
        share3_lockdown = share1_lockdown*(l1_ss/(l1_ss+l2_ss))+share2_lockdown*(l2_ss/(l1_ss+l2_ss));

        months_lockdown = 12;
        T_vac=months*30+1;
        school_open = 0;
        close_other = 0;
              

        
    elseif case_switch == 4
        % Short, too aggressive policy
        months = 24;
        effectiveness_share = 1;
        share1_lockdown = 0.20;
        share2_lockdown = 0.9;
        months_lockdown = 3;
        share3_lockdown = 1;
        T_vac=months*30+1;
        school_open = 0;
        close_other = 0;
        
    elseif case_switch == 5  % waiting for a vaccine, baseline
        % nonwork lockdown
        months = 24;
        effectiveness_share = 1;
        share1_lockdown = 0.15;
        share2_lockdown = 0.40;
        months_lockdown = 12;
        share3_lockdown = .4;
        T_vac=9*30+1;
        school_open = 0; % close schools
        close_other = 1; % more aggressive but not crazy lockdown elsewhere
        
        
    elseif case_switch == 6  % waiting for a vaccine, reopen schools
        % nonwork lockdown
        months = 24;
        effectiveness_share = 1;
        share1_lockdown = 0.35;
        share2_lockdown = 0.90;
        months_lockdown = 12;
        share3_lockdown = share1_lockdown*(l1_ss/(l1_ss+l2_ss))+share2_lockdown*(l2_ss/(l1_ss+l2_ss))+.1;
        T_vac=9*30+1;
        school_open = 1; % close schools
        close_other = 1; % more aggressive but not crazy lockdown elsewhere
        
    elseif case_switch == 7 % waiting for a vaccine
        % LOWER effectiveness of the lockdown
        months = 24;
        effectiveness_share = 0.8;
        share1_lockdown = 0.19;
        share2_lockdown = 0.50;
        months_lockdown = 12;
        share3_lockdown = share1_lockdown*(l1_ss/(l1_ss+l2_ss))+share2_lockdown*(l2_ss/(l1_ss+l2_ss))+.1;
        T_vac=9*30+1;
        school_open = 0; % close schools
        close_other = 1; % more aggressive but not crazy lockdown elsewhere
        
        
        
    elseif case_switch == 8
        % waiting for a vaccine
        % Initial R0=4
        trans_prob0 = 0.0255*1.7;             % initial transmission rate --- R0 = 2
    
        effectiveness_share = 1;
        share1_lockdown = 0.4;
        share2_lockdown = 0.9;
        months_lockdown = 15;
        share3_lockdown = share1_lockdown*(l1_ss/(l1_ss+l2_ss))+share2_lockdown*(l2_ss/(l1_ss+l2_ss))+.1;
        T_vac=9*30+1;
        school_open = 0; % close schools
        close_other = 1; % more aggressive but not crazy lockdown elsewhere
        
 
        
    end
    
    

    nperiods = months*30;        
    
    
    
    % obtain contact matrices run US contacts
    UScontacts
    Jgroups = 4;
    share_working_while_infected = 0.5;
    
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
    %omegap = 0.0*deltap;             % slope of death rate
    
    irate = 0.00001;
    
    
   
    effectiveness_share_home = effectiveness_share;
    effectiveness_share_work = effectiveness_share;
    effectiveness_share_school = effectiveness_share;
    effectiveness_share_other = effectiveness_share;

    % contact restrictions (no restrictions if set to 0) -- first
    % initialize to no restrictions then override the period in which
    % restrictoins are imposed.
    res_.school = zeros(Jgroups,nperiods);
    res_.work   = zeros(Jgroups,nperiods);
    res_.other  = zeros(Jgroups,nperiods);
    res_.home   = zeros(Jgroups,nperiods);
    
   
    res_.work(2,1:months_lockdown*30) =  share1_lockdown;
    res_.work(3,1:months_lockdown*30) =  share2_lockdown;
    
    % school_open
    if school_open
        res_.school(1,1:months_lockdown*30) = 0;
    else
        res_.school(1,1:months_lockdown*30) = share3_lockdown;
    end
    
    res_.other(1,1:months_lockdown*30) = share3_lockdown;
    
    if close_other
        res_.other(2,1:months_lockdown*30) = share3_lockdown;
        res_.other(3,1:months_lockdown*30) = share3_lockdown;
    else
        res_.other(2,1:months_lockdown*30) = share1_lockdown;
        res_.other(3,1:months_lockdown*30) = share2_lockdown;
    end
    res_.other(4,1:months_lockdown*30) = share3_lockdown;
    
    
    
    
    
    % transmission probability
    trans_prob  = trans_prob0*ones(1,nperiods);
    
    % basic contact rate
    betap0 = trans_prob0*(contacts_.home+contacts_.school+contacts_.work+contacts_.other);
    
   
    %% compute path of contact rate
    for t = 1:1:nperiods
        betamat(:,:,t) =  trans_prob(t)...
            *(contacts_.home.*(      (1-effectiveness_share_home*res_.home(:,t))*...
                                     (1-effectiveness_share_home*res_.home(:,t)') )...
            + contacts_.school.*(    (1-effectiveness_share_school*res_.school(:,t))*...
                                     (1-effectiveness_share_school*res_.school(:,t)') )...
            + contacts_.work.*(      (1-effectiveness_share_work*res_.work(:,t))*...
                                     (1-effectiveness_share_work*res_.work(:,t)') )...
            + contacts_.other.*(     (1-effectiveness_share_other*res_.other(:,t))*...
                                     (1-effectiveness_share_other*res_.other(:,t)') ) ...
              );
                                
        
        % for R0 computation check: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(09)62126-7/fulltext
        % should be the spectral radius of betap/gammap or max(abs(eig(betap(:,:,1)/gammap(1))))
        R0vec(t) = max(abs(eig(betamat(:,:,t)/gammap(1))));
    end
    
    %% solve epidemiological model
%     [Smat,Imat,Rmat,Dmat,Cmat,Nmat]...
%         = SIRmodel_Jgroups_func(betamat,nperiods,nperiods,irate,Ninit,gammap,thetap,deltap,omegap);
    

    
    
    rhop = zeros(Jgroups,nperiods);
    rhop(:,T_vac:nperiods) = rhop0*ones(Jgroups,nperiods-(T_vac-1));

    [Smat,Vmat,Imat,Rmat,Dmat,Cmat,Nmat]...
         = SIRmodel_Jgroups_vaccine_func(betamat,rhop,nperiods,irate, Ninit, gammap, thetap, deltap);
    
    
    figure
    subplot(2,1,1)
     plot(Rmat')
     subplot(2,1,2)
     
    plot(Dmat')
    
    
   
    
    % temporal aggregation of SIR model from days to months
    
    
    R1_lockdown_month = mean(reshape(Rmat(2,:)',30,months));
    R2_lockdown_month = mean(reshape(Rmat(3,:)',30,months));
    
    N1 = Ninit(2);
    N2 = Ninit(3);
    
    Dead1_lockdown_month = mean(reshape(Dmat(2,:)',30,months));
    Dead2_lockdown_month = mean(reshape(Dmat(3,:)',30,months));
    
    N1_month = mean(reshape(Nmat(2,:)',30,months));
    N2_month = mean(reshape(Nmat(3,:)',30,months));
    
    pop_path_days = -sum(Dmat);
    pop_path= mean(reshape(pop_path_days,30,months))';
  
    % path for the exogenous variable
    paramfile_two_sector_vcu
    
    % take off the infected folks from the labor force
    
    path_lockdown = [ones(months_lockdown,1); zeros(months-months_lockdown,1)];
    
    l1_path_lockdown = path_lockdown.*(-max(share1_lockdown-share_working_while_lockdown_1,0).*N1_month'...
        -min(share1_lockdown,share_working_while_lockdown_1)*(1-share_working_while_infected)*R1_lockdown_month'...
        -(1-share1_lockdown)*(1-share_working_while_infected)*R1_lockdown_month')...  % this term is 0 if the lockdown is on otherwise, it is the only term that survives and all the other terms are 0
        -(1-path_lockdown)*(1-share_working_while_infected).*R1_lockdown_month'...
        -Dead1_lockdown_month';
    
    l2_path_lockdown = path_lockdown.*(-max(share2_lockdown-share_working_while_lockdown_2,0).*N2_month'...
        -min(share2_lockdown,share_working_while_lockdown_2)*(1-share_working_while_infected)*R2_lockdown_month'...
        -(1-share2_lockdown)*(1-share_working_while_infected).*R2_lockdown_month')...
        -(1-path_lockdown)*(1-share_working_while_infected).*R2_lockdown_month'...  % this term is 0 if the lockdown is on otherwise, it is the only term that survives and all the other terms are 0
        -Dead2_lockdown_month';
    
    % should this be all the dead or only the dead among the working group?
    %pop_path = -Dead1_month'-Dead2_month';
    
    %l1_path_martin =  - I1_lockdown_month'*(1-share_working_while_infected)-(share1_lockdown-share_working_while_lockdown_1)*(N1-(1-share_working_while_infected)*I1_lockdown_month').*path_lockdown;
    %l2_path_martin =  - I2_lockdown_month'*(1-share_working_while_infected)-(share2_lockdown-share_working_while_lockdown_2)*(N2-(1-share_working_while_infected)*I2_lockdown_month').*path_lockdown;
    
    
    l1_path_lockdown_lag = [0; l1_path_lockdown(1:end-1)];
    l2_path_lockdown_lag = [0; l2_path_lockdown(1:end-1)];
    pop_path_lag = [0; pop_path(1:end-1)];
    
    % define matrix needed to calculate the sizes of the shocks
    amat = diag(rhop_l*ones(months,1));
    
    % translate the path for a into the path for the shocks m
    m_path = l1_path_lockdown-amat*l1_path_lockdown_lag;
    n_path = l2_path_lockdown-amat*l2_path_lockdown_lag;
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
    
    % now get the response to the initial conditions
    shockssequence = 0;
    irfshock = '';
    
    
    
    % [QUESTION FOR MARTIN what is this about?]
    %changes_names = char('rhop');
    %changes_vec =       [-3/2];
    
    if capital_irreversibility_switch
        
        [zdatalinear, zdatapiecewise, zdatass, oobase_ Mbase_] = ...
            solve_one_constraint(modnam,modnamstar,...
            constraint, constraint_relax,...
            shockssequence,irfshock,months,maxiter,init)
        
        zdatalinear = zdatapiecewise;
        
    else
        
        [zdatalinear, zdatass, oobase_, Mbase_ ] = ...
            solve_no_constraint(modnam,...
            shockssequence,irfshock,months,init);
    end
    
    
    %% now repeat without lockdown
    
    months_lockdown = 0;
    
    gammap = 0.2*ones(Jgroups,1);    % transition rate from I to R
    thetap = 0.1*ones(Jgroups,1);    % transition rate from R to C and D
    deltap = 0.01*ones(Jgroups,1);   % death rate
    omegap = 0.0*deltap;             % slope of death rate
    
    
    % contact restrictions (no restrictions if set to 0) -- first
    % initialize to no restrictions then override the period in which
    % restrictoins are imposed.
    res_.school = zeros(Jgroups,nperiods);
    res_.work   = zeros(Jgroups,nperiods);
    res_.other  = zeros(Jgroups,nperiods);
    res_.home   = zeros(Jgroups,nperiods);

    
    % transmission probability
    trans_prob  = trans_prob0*ones(1,nperiods);

    % basic contact rate
    betap0 = trans_prob0*(contacts_.home+contacts_.school+contacts_.work+contacts_.other);
    
   
    %% compute path of contact rate
    for t = 1:1:nperiods
        betamat(:,:,t) =  trans_prob(t)...
            *(contacts_.home.*(      (1-effectiveness_share_home*res_.home(:,t))*...
                                     (1-effectiveness_share_home*res_.home(:,t)') )...
            + contacts_.school.*(    (1-effectiveness_share_school*res_.school(:,t))*...
                                     (1-effectiveness_share_school*res_.school(:,t)') )...
            + contacts_.work.*(      (1-effectiveness_share_work*res_.work(:,t))*...
                                     (1-effectiveness_share_work*res_.work(:,t)') )...
            + contacts_.other.*(     (1-effectiveness_share_other*res_.other(:,t))*...
                                     (1-effectiveness_share_other*res_.other(:,t)') ) ...
              );
                                
        
        % for R0 computation check: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(09)62126-7/fulltext
        % should be the spectral radius of betap/gammap or max(abs(eig(betap(:,:,1)/gammap(1))))
        R0vec(t) = max(abs(eig(betamat(:,:,t)/gammap(1))));
    end
    
    %% solve epidemiological model
   
    rhop = zeros(Jgroups,nperiods);
    %rhop(:,T_vac:nperiods) = rhop0*ones(Jgroups,nperiods-(T_vac-1));

    [Smat,Vmat,Imat,Rmat,Dmat,Cmat,Nmat]...
         = SIRmodel_Jgroups_vaccine_func(betamat,rhop,nperiods,irate, Ninit, gammap, thetap, deltap);
    
    
    figure
    subplot(2,1,1)
    plot(Rmat')
     
    subplot(2,1,2) 
    plot(Dmat')
    
    
    % temporal aggregation of SIR model from days to months
    
    
    R1_month = mean(reshape(Rmat(2,:)',30,months));
    R2_month = mean(reshape(Rmat(3,:)',30,months));
    
    N1 = Ninit(2);
    N2 = Ninit(3);
    
    Dead1_month = mean(reshape(Dmat(2,:)',30,months));
    Dead2_month = mean(reshape(Dmat(3,:)',30,months));
    
    N1_month = mean(reshape(Nmat(2,:)',30,months));
    N2_month = mean(reshape(Nmat(3,:)',30,months));
    
    pop_path_days = -sum(Dmat);
    pop_path= mean(reshape(pop_path_days,30,months))';
  
    % path for the exogenous variable
    paramfile_two_sector_vcu
    
    % take off the infected folks from the labor force
    
    
    path_lockdown = [ones(months_lockdown,1); zeros(months-months_lockdown,1)];
    
    l1_path = path_lockdown.*(-max(share1_lockdown-share_working_while_lockdown_1,0).*N1_month'...
        -min(share1_lockdown,share_working_while_lockdown_1)*(1-share_working_while_infected)*R1_month'...
        -(1-share1_lockdown)*(1-share_working_while_infected)*R1_month')...  % this term is 0 if the lockdown is on otherwise, it is the only term that survives and all the other terms are 0
        -(1-path_lockdown)*(1-share_working_while_infected).*R1_month'...
        -Dead1_month';
    
    l2_path = path_lockdown.*(-max(share2_lockdown-share_working_while_lockdown_2,0).*N2_month'...
        -min(share2_lockdown,share_working_while_lockdown_2)*(1-share_working_while_infected)*R2_month'...
        -(1-share2_lockdown)*(1-share_working_while_infected).*R2_month')...
        -(1-path_lockdown)*(1-share_working_while_infected).*R2_month'...  % this term is 0 if the lockdown is on otherwise, it is the only term that survives and all the other terms are 0
        -Dead2_month';
    
    % should this be all the dead or only the dead among the working group?
    %pop_path = -Dead1_month'-Dead2_month';
    
    %l1_path_martin =  - I1_lockdown_month'*(1-share_working_while_infected)-(share1_lockdown-share_working_while_lockdown_1)*(N1-(1-share_working_while_infected)*I1_lockdown_month').*path_lockdown;
    %l2_path_martin =  - I2_lockdown_month'*(1-share_working_while_infected)-(share2_lockdown-share_working_while_lockdown_2)*(N2-(1-share_working_while_infected)*I2_lockdown_month').*path_lockdown;
    
    
    l1_path_lag = [0; l1_path(1:end-1)];
    l2_path_lag = [0; l2_path(1:end-1)];
    pop_path_lag = [0; pop_path(1:end-1)];
    
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
    
    % now get the response to the initial conditions
    shockssequence = 0;
    irfshock = '';
    nperiods = months;
    
    
    
    
    if capital_irreversibility_switch
        
        [zdatalinear2, zdatapiecewise2, zdatass2, oobase2_ Mbase2_] = ...
            solve_one_constraint(modnam,modnamstar,...
            constraint, constraint_relax,...
            shockssequence,irfshock,nperiods,maxiter,init)
        
        zdatalinear2 = zdatapiecewise2;
        
    else
        
        [zdatalinear2, zdatass2, oobase2_, Mbase2_ ] = ...
            solve_no_constraint(modnam,...
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
    
    
    
    

        
        %% Aggregate figure
        titlelist = char('Group-1 Recovering','Group-2 Recovering',...
            'Group-1 Deceased, cumulative','Group-2 Deceased, cumulative',...
            'Output, per capita','Capacity Utilization',...
            'Consumption, per capita','Investment, per capita');
        
        
        line1 = 100*[ R1_lockdown_month'/s1_l1_ss, R2_lockdown_month'/s2_l2_ss, ...
            Dead1_lockdown_month'/s1_l1_ss, Dead2_lockdown_month'/s2_l2_ss, ...
            s1_ypc_irf/s1_ypc_ss, s1_u_irf/s1_u_ss,...
            s1_cpc_irf/s1_cpc_ss, s1_inpc_irf/s1_inpc_ss];
        
        line2 = 100*[ R1_month'/s2_l1_ss, R2_month'/s2_l2_ss,...
            Dead1_month'/s1_l1_ss, Dead2_month'/s2_l2_ss, ...
            s2_ypc_irf/s2_y_ss, s2_u_irf/s2_u_ss,...
            s2_cpc_irf/s2_c_ss, s2_inpc_irf/s2_in_ss];
        
        percent = '% dev. from steady state';
        ppt = 'Perc. Point';
        level = 'Level';
        percent_group = '% of group size';
        
        ylabels = char(percent_group,percent_group,...
            percent_group, percent_group,...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Do Nothing','Social Distancing');
        makechart(titlelist,legendlist,figlabel,ylabels,line2,line1);
        subplot(4,2,7)
        xlabel('Months')
        subplot(4,2,8)
        xlabel('Months')
       
        
        
        makechart(titlelist,legendlist,figlabel,ylabels,...
            month2quarter([zeros(2,size(line2,2));line2]),...
            month2quarter([zeros(2,size(line1,2));line1]));
        subplot(4,2,7)
        xlabel('Quarters')
        subplot(4,2,8)
        xlabel('Quarters')
        

        
        
        printpref
        
        eval(['print -depsc2 figure1_social_distancing_vs_do_nothing_case',num2str(case_switch),'.eps'])
        
        
        %% slides
        
        %% Aggregate figure
        titlelist = char('Group-1 Recovering','Group-2 Recovering',...
            'Group-1 Deceased, cumulative','Group-2 Deceased, cumulative');
        
        
        line1 = 100*[ R1_lockdown_month'/s1_l1_ss, R2_lockdown_month'/s2_l2_ss, ...
            Dead1_lockdown_month'/s1_l1_ss, Dead2_lockdown_month'/s2_l2_ss];
        
        line2 = 100*[ R1_month'/s2_l1_ss, R2_month'/s2_l2_ss,...
            Dead1_month'/s1_l1_ss, Dead2_month'/s2_l2_ss];
        
        percent = '% dev. from steady state';
        ppt = 'Perc. Point';
        level = 'Level';
        percent_group = '% of group size';
        
        ylabels = char(percent_group,percent_group,...
            percent_group, percent_group);
        figlabel = '';
        legendlist = char('Do Nothing','Social Distancing');
        makechart(titlelist,legendlist,figlabel,ylabels,...
            month2quarter([zeros(2,size(line2,2));line2]),...
            month2quarter([zeros(2,size(line1,2));line1]));
        subplot(2,2,3)
        xlabel('Quarters')
        subplot(2,2,4)
        xlabel('Quarters')
       
        
        
        titlelist = char('Output, per capita','Capacity Utilization',...
            'Consumption, per capita','Investment, per capita');
        
        
        line1 = 100*[ s1_ypc_irf/s1_ypc_ss, s1_u_irf/s1_u_ss,...
            s1_cpc_irf/s1_cpc_ss, s1_inpc_irf/s1_inpc_ss];
        
        line2 = 100*[ s2_ypc_irf/s2_y_ss, s2_u_irf/s2_u_ss,...
            s2_cpc_irf/s2_c_ss, s2_inpc_irf/s2_in_ss];
        
        percent = '% dev. from steady state';
        ppt = 'Perc. Point';
        level = 'Level';
        percent_group = '% of group size';
        
        ylabels = char(percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Do Nothing','Social Distancing');
        makechart(titlelist,legendlist,figlabel,ylabels,...
            month2quarter([zeros(2,size(line2,2));line2]),...
            month2quarter([zeros(2,size(line1,2));line1]));
        subplot(2,2,3)
        xlabel('Quarters')
        subplot(2,2,4)
        xlabel('Quarters')
       
        %% sectoral detail
        
        titlelist = char('Group-1 Recovering','Group-2 Recovering',...
            'Labor, Sector 1', 'Labor, Sector 2',...
            'Value Added, Sector 1', 'Value Added, Sector 2');
        
        
        line1 = 100*[ R1_lockdown_month'/s1_l1_ss, R2_lockdown_month'/s1_l2_ss, ...
            s1_l1_irf/s1_l1_ss, s1_l2_irf/s1_l2_ss,...
            s1_v1_irf/s1_v1_ss, s1_v2_irf/s1_v2_ss...
            ];
        
        line2 = 100*[ R1_month'/s2_l1_ss, R2_month'/s2_l2_ss, ...
            s2_l1_irf/s2_l1_ss, s2_l2_irf/s2_l2_ss,...
            s2_v1_irf/s2_v1_ss, s2_v2_irf/s2_v2_ss...
            ];
        
        
        ylabels = char(percent_group,percent_group,...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Do Nothing','Social Distancing');
        makechart(titlelist,legendlist,figlabel,ylabels,line2,line1);
        subplot(3,2,5)
        xlabel('Months')
        subplot(3,2,6)
        xlabel('Months')
        
          makechart(titlelist,legendlist,figlabel,ylabels,...
            month2quarter([zeros(2,size(line2,2));line2]),...
            month2quarter([zeros(2,size(line1,2));line1]));
        subplot(3,2,5)
        xlabel('Quarters')
        subplot(3,2,6)
        xlabel('Quarters')
        
        
        printpref
        eval(['print -depsc2 figure2_social_distancing_vs_do_nothing_case',num2str(case_switch),'.eps'])
        

end



