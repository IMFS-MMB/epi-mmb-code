%%% Epi-mmb note: We chang in=1e-12 in one_sector_vcu_binding.mod and two_sector_vcu_binding.mod 
%%% so that the result can be replicated with Dynare 4.6 and 5.1 . 

%setpath_dynare
addpath('toolkit_files')  
global M_ oo_ changes_names changes_vec

modnam = 'two_sector_vcu';
modnamstar= 'two_sector_vcu_binding';

modnam_one_sector = 'one_sector_vcu';
modnamstar_one_sector= 'one_sector_vcu_binding';

capital_irreversibility_switch = 1;

%case_switch = 1;  % 1, baseline, 2 cheapest, 3 shift lockdown to outside the labor force, 4 too strict
constraint = 'in<-in_ss';
constraint_relax = 'lambdai<-0.0000001';

maxiter = 20;


slide_switch = 0;



jv_table = readtable('JFV_Jones_DetailedResults_transformed.csv');

us_jv_table = jv_table(strmatch('United States',jv_table{:,'region'},'exact'),:);

mydates_R = us_jv_table{:,'mydate'};
R = us_jv_table{:,'R'};

%{
figure
plot(datetime(datestr(mydates_R)),movmean(R,7),'k-','LineWidth',2)
%}
pos_apr_1 = find(datetime(datestr(mydates_R))==datetime(datestr('1-Apr-2020')));


observed_R0 = movmean(R,7);

observed_R0 = observed_R0(pos_apr_1:end);


%eval(['dynare ',modnam,' noclearall'])

for case_switch = 1:1
    paramfile_two_sector_vcu;
      
    months = 24;
    months_lockdown = 3;
    nperiods = months*30;            
    nperiods_plot = 360;
    
    rhop0 = 0.2/30;
    
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
        T_vac=nperiods+1;
        school_open = 0;
        close_other = 0;
        
    
    end
    
    

       
    
    
    
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
    
    %irate = 1-0.999116189042276;
    helper=load('inf_ini.mat');
    irate=helper.helper;
    
    
   
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
    
    share1_add_lockdown = 0.1;
    share2_add_lockdown = 0.26;
    share_add_lockdown = share1_add_lockdown*(l1_ss/(l1_ss+l2_ss))+share2_add_lockdown*(l2_ss/(l1_ss+l2_ss));
  
    res_.work(2,1:15) = linspace(0,1,15)*(share1_add_lockdown+share1_lockdown);
    res_.work(3,1:15) = linspace(0,1,15)*(share2_add_lockdown+share2_lockdown);
    
    res_.work(2,16:75) = share1_lockdown+share1_add_lockdown - linspace(0,1,60)*share1_add_lockdown;
    res_.work(3,16:75) = share2_lockdown+share2_add_lockdown - linspace(0,1,60)*share2_add_lockdown;
    
    % start the sim in March --- Keep schools open in March
    res_.school(1,1:months_lockdown*30) = 0;
    % close them through the end of July
    res_.school(1,15:5*30) = 1;
    
    % reopen in August at 50%
    res_.school(1,5*30+1:end) = 0.5;
    
    res_.other(1,1:months_lockdown*30) = share3_lockdown;
    res_.other(2,1:months_lockdown*30) = share1_lockdown;
    res_.other(3,1:months_lockdown*30) = share2_lockdown;
    res_.other(4,1:months_lockdown*30) = share3_lockdown;

    

    
    % transmission probability
    trans_prob = trans_prob0*ones(1,nperiods);
%    trans_prob(20:120) = 0.7*trans_prob0;
%    trans_prob(120:300) = 0.7*trans_prob0;
%    trans_prob(301:330) = 0.8*trans_prob0;
%    trans_prob(331:360) = trans_prob0;
    
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
        R0vecOrig(t) = max(abs(eig(betamat(:,:,t)/gammap(1))));
        betamatOrig = betamat;
        
        R0vec(t) = R0vecOrig(t);
        if (t>15 & t<=length(observed_R0)+15)
           scale_factor = R0vec(t)/observed_R0(t-15);
           betamat(:,:,t) =  betamat(:,:,t)/scale_factor;
           R0vec(t) = max(abs(eig(betamat(:,:,t)/gammap(1)))); 
        elseif t>length(observed_R0)+15
            
           scale_factor = R0vec(t)/observed_R0(end);
           betamat(:,:,t) =  betamat(:,:,t)/scale_factor;
           R0vec(t) = max(abs(eig(betamat(:,:,t)/gammap(1)))); 
        end
    end
    
    %% solve epidemiological model
%     [Smat,Imat,Rmat,Dmat,Cmat,Nmat]...
%         = SIRmodel_Jgroups_func(betamat,nperiods,nperiods,irate,Ninit,gammap,thetap,deltap,omegap);
    %%
    %{
    figure
    
    if slide_switch
       subplot(2,2,1)
    else
       subplot(4,1,1)
    end
    pos_mar_16 = find(datetime(datestr(mydates_R))==datetime(datestr('16-Mar-2020')));
    
    nobs = length(mydates_R(pos_mar_16:end));
    plot(datetime(datestr(mydates_R(pos_mar_16:end))),R0vecOrig(1:nobs),'r--','linewidth',1.5)
    hold on
    
    plot(datetime(datestr(mydates_R(pos_mar_16:end))),movmean(R(pos_mar_16:end),7),'k-','linewidth',1.5)
    
    legend('Reproduction, model implied','Reproduction rate, observed')
   
    %}
    
    % input for the vaccine (assume it is irrelevant)
    rhop = zeros(Jgroups,nperiods);
    rhop(:,T_vac:nperiods) = rhop0*ones(Jgroups,nperiods-(T_vac-1));

    [Smat,Vmat,Imat,Rmat,Dmat,Cmat,Nmat]...
         = SIRmodel_Jgroups_vaccine_func(betamat,rhop,nperiods,irate, Ninit, gammap, thetap, deltap);
    
    %[Smat,Vmat,Imat,Rmat,Dmat,Cmat,Nmat]...
    %     = SIRmodel_Jgroups_vaccine_func(betamatOrig,rhop,nperiods,irate, Ninit, gammap, thetap, deltap);
    
     
    % make plots  
    %{
    
    months_dates = datetime(datestr(datenum(2020,3:27,1)));
    
    if slide_switch
       subplot(2,2,2)
    else
       subplot(4,1,2)
    end
    plot(months_dates,day2month([ones(15,1); sum(Smat)'])*100,'k','linewidth',1.5); hold on
   
    plot(months_dates,day2month([zeros(15,1); sum(Cmat)'])*100,'r--','linewidth',1.5); 
    
    legend('Susceptible','Cured')
    ylabel('Percent of population')
    xtickformat('MMM')
    xlim([months_dates(1) months_dates(8)])
    
    if slide_switch
       subplot(2,2,3)
    else
       subplot(4,1,3)
    end
    plot(months_dates,day2month([zeros(15,1); sum(Imat)']*100000),'k','linewidth',1.5); hold on
    
    plot(months_dates,day2month([zeros(15,1); sum(Rmat)']*100000),'r--','linewidth',1.5)
    legend('Infected','Recovering')
    ylabel('Average occurrences per 100,000')
    xtickformat('MMM')
    xlim([months_dates(1) months_dates(8)])
    
    if slide_switch
       subplot(2,2,4)
    else
       subplot(4,1,4)
    end
    plot(months_dates,day2month([zeros(15,1); sum(Dmat)']*330000000),'k','linewidth',1.5)
    legend('Deceased')
%  
    ylabel('Total')
    xtickformat('MMM')
    xlim([months_dates(1) months_dates(8)])
        %}
    %% temporal aggregation of SIR model from days to months
    
    
    R1_lockdown_month = mean(reshape(Rmat(2,:)',30,months));
    R2_lockdown_month = mean(reshape(Rmat(3,:)',30,months));
    
    

    pop_path_days = [zeros(1,15),-sum(Dmat(:,1:end-15))];
    pop_path= mean(reshape(pop_path_days,30,months))';
  
    % path for the exogenous variable
    paramfile_two_sector_vcu
    
    % take off the infected folks from the labor force
    
    


    l1_path_lockdown_day = -max(res_.work(2,:)'-share_working_while_lockdown_1,0).*(Nmat(2,:)'-Dmat(2,:)')...
        -min(res_.work(2,:)',share_working_while_lockdown_1)*(1-share_working_while_infected).*Rmat(2,:)'...
        -(1-res_.work(2,:)')*(1-share_working_while_infected).*Rmat(2,:)'...  
        -Dmat(2,:)';

    l1_path_lockdown = day2month([zeros(15,1);l1_path_lockdown_day(1:end-15)]);

    l2_path_lockdown_day = -max(res_.work(3,:)'-share_working_while_lockdown_2,0).*(Nmat(3,:)'-Dmat(3,:)')...
        -min(res_.work(3,:)',share_working_while_lockdown_2)*(1-share_working_while_infected).*Rmat(3,:)'...
        -(1-res_.work(3,:)')*(1-share_working_while_infected).*Rmat(3,:)'...  
        -Dmat(3,:)';

    l2_path_lockdown = day2month([zeros(15,1);l2_path_lockdown_day(1:end-15)]);

%{
    figure
    
    subplot(2,1,1)
    plot(res_.work(2,:))
    hold on
    plot(res_.work(3,:))
    
    subplot(2,1,2)
    
    plot(months_dates(1:24),-(l1_path_lockdown+l2_path_lockdown)/(l1_ss+l2_ss)*100); hold on
    
    ur_table = readtable('UNRATE.xls'); 

    plot(ur_table.observation_date,ur_table.UNRATEDIFF)
    
    legend('Shock Path','Unemployment Rate (deviation from February)')
 %}   
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
    
    
    
    if capital_irreversibility_switch
        
        [zdatalinear, zdatapiecewise, zdatass, oobase_ Mbase_] = ...
            solve_one_constraint(modnam,modnamstar,...
            constraint, constraint_relax,...
            shockssequence,irfshock,months,maxiter,init);
        
        zdatalinear = zdatapiecewise;
        
    else
        
        [zdatalinear, zdatass, oobase_, Mbase_ ] = ...
            solve_no_constraint(modnam,...
            shockssequence,irfshock,months,init);
    end
    
    
    N1 = 0;
    N2 = l1_ss+l2_ss;
    
    
    
    l2_path = l1_path_lockdown+l2_path_lockdown;
    
    l2_path_lag = [0; l2_path(1:end-1)];
    
    
    n_path = l2_path-amat*l2_path_lag;
    
    
    eval(['dynare ',modnam_one_sector,' noclearall'])
    
    endo_names = cellstr(M_.endo_names);
    for m_indx = 1:months
        
        this_n = ['n',num2str(m_indx)];
        n_pos_vec(m_indx)=find(strcmp(this_n, endo_names));
        this_pop_shock = ['pop_shock',num2str(m_indx)];
        pop_shock_pos_vec(m_indx) = find(strcmp(this_pop_shock, endo_names));
    end
    
    init = zeros(M_.endo_nbr,1);
    init(n_pos_vec)=n_path;
    init(pop_shock_pos_vec) = pop_shock_path;
    
    
    
    if capital_irreversibility_switch
        
        [zdatalinear2, zdatapiecewise2, zdatass2, oobase2_ Mbase2_] = ...
            solve_one_constraint(modnam_one_sector,modnamstar_one_sector,...
            constraint, constraint_relax,...
            shockssequence,irfshock,months,maxiter,init);
        
        zdatalinear2 = zdatapiecewise2;
        
    else
        
        [zdatalinear2, zdatass2, oobase2_, Mbase2_ ] = ...
            solve_no_constraint(modnam_one_sector,...
            shockssequence,irfshock,months,init);
    end
    
    
    
    
    
    prefix1 = 's1_';
    prefix2 = 's2_';
    for this_var = 1:Mbase_.endo_nbr
        eval([prefix1,deblank(Mbase_.endo_names{this_var,:}),'_irf=zdatalinear(:,this_var);'])
        eval([prefix1,deblank(Mbase_.endo_names{this_var,:}),'_ss=zdatass(this_var);'])
    end
    
    for this_var = 1:Mbase2_.endo_nbr
        eval([prefix2,deblank(Mbase2_.endo_names{this_var,:}),'_irf=zdatalinear2(:,this_var);'])
        eval([prefix2,deblank(Mbase2_.endo_names{this_var,:}),'_ss=zdatass2(this_var);'])
    end
    
    
    
    %{
    if slide_switch
        
        
        %% Aggregate figure
        titlelist = char(...
            'Output, Total','Capacity Utilization',...
            'Consumption','Investment');
        
        
        line1 = 100*[ ...
            s1_y_irf/s1_y_ss, s1_u_irf/s1_u_ss,...
            s1_c_irf/s1_c_ss, s1_in_irf/s1_in_ss];
        
        line2 = 100*[...
            s2_y_irf/s2_y_ss, s2_u_irf/s2_u_ss,...
            s2_c_irf/s2_c_ss, s2_in_irf/s2_in_ss];
        
        percent = '% dev. from steady state';
        ppt = 'Perc. Point';
        level = 'Level';
        percent_group = '% of group size';
        
        ylabels = char(...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Two-Sector Model','One-Sector Model');
        makechart(titlelist,legendlist,figlabel,ylabels,month2quarter([zeros(2,size(line1,2));line1]),month2quarter([zeros(2,size(line2,2));line2]));
        subplot(2,2,3)
        xlabel('Quarters')
        subplot(2,2,4)
        xlabel('Quarters')
        
        
        printpref_slide
        
        
        %% sectoral detail
        
        titlelist = char(...
            'Labor, Sector 1', 'Labor, Sector 2',...
            'Value Added, Sector 1', 'Value Added, Sector 2');
        
        
        line1 = 100*[ ...
            s1_l1_irf/s1_l1_ss, s1_l2_irf/s1_l2_ss,...
            s1_v1_irf/s1_v1_ss, s1_v2_irf/s1_v2_ss...
            ];
        
  
        
        
        ylabels = char(percent_group,percent_group,...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = '';
        makechart(titlelist,legendlist,figlabel,ylabels,month2quarter([zeros(2,size(line1,2));line1]));
        subplot(2,2,3)
        xlabel('Quarters')
        subplot(2,2,4)
        xlabel('Quarters')
        
        printpref_slide
        
        
    else % figures for paper
        
        %% Aggregate figure
        titlelist = char('Recovering (per 100,000)','Deaths (total)',...
            'Output, Total','Capacity Utilization',...
            'Consumption','Investment');
        
        
        dead_cum = day2month(sum(Dmat)');
       
        line1 = 100*[day2month(sum(Rmat)'*1000), dead_cum*3300000, ...
            s1_y_irf/s1_y_ss, s1_u_irf/s1_u_ss,...
            s1_c_irf/s1_c_ss, s1_in_irf/s1_in_ss];
        
        line2 = 100*[day2month(sum(Rmat)'*1000), dead_cum*3300000,...
            s2_y_irf/s2_y_ss, s2_u_irf/s2_u_ss,...
            s2_c_irf/s2_c_ss, s2_in_irf/s2_in_ss];
        
        percent = '% dev. from steady state';
        ppt = 'Perc. Point';
        level = 'Level';
        percent_group = '% of group size';
        
        ylabels = char('','',...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Two-Sector Model','One-Sector Model');
        makechart(titlelist,legendlist,figlabel,ylabels,line1,line2);
        subplot(3,2,5)
        xlabel('Months')
        subplot(3,2,6)
        xlabel('Months')
        
        
        % create inset figure with quarterly path for GDP
        subplot(3,2,3) 
        
        
        h = axes('Parent', gcf, 'Position', [0.27   0.475    0.1708    0.1144]);
    
        %p = get(gca, 'Position');
        %h = axes('Parent', gcf, 'Position', [p(1)+.07 p(2)+.0 p(3)-.1 p(4)-.1]);
        
        s1_y_irf_q = month2quarter([zeros(2,1);s1_y_irf/s1_y_ss]*100);
        
        s2_y_irf_q = month2quarter([zeros(2,1);s2_y_irf/s2_y_ss]*100);
        
        
        plot(h, [1:length(s1_y_irf_q)]', s1_y_irf_q, 'b-','linewidth',2)
        
        xlim([1 8])
        
        hold on 
        plot(h, [1:length(s2_y_irf_q)]', s2_y_irf_q,'r--','linewidth',2);
        
        
        box on
        
        xlabel('Quarters')
        
        
        
        
        printpref
        
        %eval(['print -depsc2 figure1_sectoral_differences',num2str(case_switch),'.eps'])
        %%
        
        makechart(titlelist,legendlist,figlabel,ylabels,month2quarter([zeros(2,size(line1,2));line1]),month2quarter([zeros(2,size(line2,2));line2]));
        subplot(3,2,5)
        xlabel('Quarters')
        subplot(3,2,6)
        xlabel('Quarters')
        
        
        %% sectoral detail
        
        titlelist = char('Group-1 Recovering','Group-2 Recovering',...
            'Labor, Sector 1', 'Labor, Sector 2',...
            'Value Added, Sector 1', 'Value Added, Sector 2');
        
        
        
        line1 = 100*[ R1_lockdown_month'/s1_l1_ss, R2_lockdown_month'/s1_l2_ss, ...
            s1_l1_irf/s1_l1_ss, s1_l2_irf/s1_l2_ss,...
            s1_v1_irf/s1_v1_ss, s1_v2_irf/s1_v2_ss...
            ];
        
        
        
        
        ylabels = char(percent_group,percent_group,...
            percent,percent,...
            percent,percent);
        figlabel = '';
        legendlist = char('Two-Sector Model');
        makechart(titlelist,legendlist,figlabel,ylabels,line1);
        subplot(3,2,5)
        xlabel('Months')
        subplot(3,2,6)
        xlabel('Months')
        
        
        
        
        printpref
        %eval(['print -depsc2 figure2_social_distancing_vs_do_nothing_case',num2str(case_switch),'.eps'])
        
        makechart(titlelist,legendlist,figlabel,ylabels,month2quarter([zeros(2,size(line1,2));line1]));
        subplot(3,2,5)
        xlabel('Quarters')
        subplot(3,2,6)
        xlabel('Quarters')
        
        
    end
    %}
end




    
