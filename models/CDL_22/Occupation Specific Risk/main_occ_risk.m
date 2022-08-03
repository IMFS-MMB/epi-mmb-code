% File to create Figure 10 in paper where each occupation has their own
% transmission probability to infect in market work

clc;
clear all;

fig_dir = fileparts('C:\Users\17786\OneDrive\Desktop\Data\Figures\');
% Specify directory where figure will get saved
main_dir = fileparts('C:\Users\17786\OneDrive\Desktop\Data\');
main_results = fullfile(main_dir,'results.mat');
%% Inputting data on weights and elasticities
% Output formulation in a flexible way but for calibration we use weights
% and shares such that output is a linear sum across all occupations
% One can change the weights and shares to get a more flexible output
% formulation

%Sector specific weights and elasticities to get GDP
sec_param = csvread('gdp_sectors_parameters.csv',1,0); % gives both sector weights and elascticities
sec_weight = sec_param(:,1); % sector weights
sec_elas = sec_param(1,2); % sector elasticities
n_sec = size(sec_weight,1); % number of sectors

%Occupation specific weights and elasticities to get sectoral GDP
occ_elas = csvread('sectors_occ_elasticities.csv',1,1); % gives both sector weights and elascticities
%row are sectors
occ_weight = csvread('sectors_occ_weights.csv',1,1); % gives both sector weights and elascticities
%row are sectors and columns are occupations
n_occ = size(occ_weight,2);

periods = 250; % number of periods to simulate
tau_m = zeros(periods+1,n_occ); % tax on market wage to generate lockdowns
% Inputting wages and risk factors

[ID,NOC,NOC_short,employment,~,~,market_wage,dummy_home,rel_wage,home_wage,Riskiness,Relativehomewage] = readvars('wage_schedule.xlsx');
NOC = string(NOC);
NOC_short = string(NOC_short);

ini_pop = employment'./sum(employment'); % Initial population in each occupation

wage_sm = repmat(market_wage',periods+1,1); % wage of each occupation at market over time
wage_sh = repmat(home_wage',periods+1,1); % wage of each occupation at home over time
theta_ss =  length(find(dummy_home(:,1)==0))/size(dummy_home(:,1),1);

gamma_Ih = 0.8; %relative productivity of infected to susceptible for working at home
gamma_Im = 0.8; %relative productivity of infected to susceptible for working at market
[wage_r,pol_mkt_r,wage_im,wage_ih,wage_i,pol_mkt_i,theta_r,theta_i] = steady_state_cutoff(periods,ini_pop,wage_sm,wage_sh,tau_m,gamma_Im,gamma_Ih);

risk_prop = Riskiness'; 
pol_mkt_s = pol_mkt_r;
wage_s = wage_r;

% Government revenue at steady state
G_ss = sum(tau_m.*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r==2) ,2);
C_ss = sum((1-tau_m).*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r==2) + wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r==1) ,2)+  G_ss ;
%consumption at steady state

%output at steady state
y_sec = NaN*ones(periods+1,n_sec);
y_occ = NaN*ones(periods+1,n_sec,n_occ); % output at time t, in sector j in occupation i
for t = 1:periods+1
    for j = 1:n_sec
    y_occ(t,j,:) = wage_r(t,:).*ini_pop;
    y_sec(t,j) = sum(occ_weight(j,:).*(reshape(y_occ(t,j,:),1,[]).^(occ_elas(j))) ).^(1/occ_elas(j));
    end
end
Y_ss = sum( repmat(sec_weight',periods+1,1).*(y_sec.^sec_elas),2).^(1./sec_elas);

% remember that Y_ss is not equal to C_ss so need to solve for that
% constant such that they equalize

opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve

betta=0.96^(1/52);  %Weekly household discount factor
pid = (7*0.005)/18; % Weekly probability of dying
pir=7/18-pid;     %Weekly probability of recovering
epsilon = 0.001;  %initial infected

%Calibation targets for shares of pis-terms in T-function in SIR model
pis3_shr_target=2/3;                   %share of T_0 jump due general infections
pis1_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to consumption-based infections
pis2_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to work-based infections
RplusD_target=0.60;                    %total share of people infected and then either recovered or dead after epidemic

go_calibrate_pis;

% storing infection parameters
pi_par(1) = pir;
pi_par(2) = pi;
pi_par(3) = varepsilon;
pi_par(4) = epsilon;
pi_par(5) = pic;
pi_par(6) = pid;

pol_mkt_s_exo = pol_mkt_r;
wage_s_exo = wage_r;
[~,~,~,I_exo,~,~,S_exo,~,~,R_exo,D_exo,~,~,~,~,~,~,~,C_exo,G_exo,Y_exo,theta_total_exo] ...
    = sir_dynamics(pol_mkt_s_exo,wage_s_exo,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,...
    pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given NO endogeneous responses
C_dev_exo = (log(C_exo) - log(C_ss))*100;

pol_mkt_s_init = pol_mkt_r;
wage_s_init = wage_r;

[pol_mkt_s_opt,~,RnotSIR_endog,I_endog,~,~,S_endog,~,~,R_endog,D_endog,~,~,~,~,Cs_endog,Ci_endog,Cr_endog,C_endog,G_endog,Y_endog,theta_total_endog,error_endog] ...
    = policy_susceptible(pol_mkt_s_init,wage_s_init,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog = (log(C_endog)- log(C_ss))*100; %deviation of total consumption

%% Changing tau's to create various lockdown scenarios (Figure 5,6)
tau_m_val = [0.05; 0.12; 0.20];
dur_tax = [36;12;52];
tau_size = size(tau_m_val,1);
tau_m_pol = zeros(periods+1,n_occ,tau_size);

%pre-allocate endogenous and exogenous output
C_ss_tau = NaN*ones(periods+1,tau_size);
G_ss_tau = NaN*ones(periods+1,tau_size);

theta_total_exo_tau = NaN*ones(periods+1,tau_size);
pol_mkt_s_tau = NaN*ones(periods+1,n_occ,tau_size);
pol_mkt_r_tau = NaN*ones(periods+1,n_occ,tau_size);
pol_mkt_i_tau = NaN*ones(periods+1,n_occ,tau_size);

I_exo_tau = NaN*ones(periods+1,tau_size);
S_exo_tau = NaN*ones(periods+1,tau_size);
R_exo_tau = NaN*ones(periods+1,tau_size);
D_exo_tau = NaN*ones(periods+1,tau_size);
C_exo_tau = NaN*ones(periods+1,tau_size);
G_exo_tau = NaN*ones(periods+1,tau_size);
Y_exo_tau = NaN*ones(periods+1,tau_size);
C_dev_exo_tau = NaN*ones(periods+1,tau_size);
C_pv_exo_tau = NaN*ones(1,tau_size);

pol_mkt_s_opt_tau = NaN*ones(periods+1,n_occ,tau_size);
theta_opt_tau = NaN*ones(periods+1,tau_size);
theta_total_endog_tau = NaN*ones(periods+1,tau_size);
I_endog_tau = NaN*ones(periods+1,tau_size);
S_endog_tau = NaN*ones(periods+1,tau_size);
R_endog_tau = NaN*ones(periods+1,tau_size);
D_endog_tau = NaN*ones(periods+1,tau_size);
C_endog_tau = NaN*ones(periods+1,tau_size);
G_endog_tau = NaN*ones(periods+1,tau_size);
Y_endog_tau = NaN*ones(periods+1,tau_size);
C_dev_endog_tau = NaN*ones(periods+1,tau_size);
C_pv_endog_tau = NaN*ones(1,tau_size);

for ii = 1:tau_size
    tau_m_pol(1:dur_tax(ii),find(pol_mkt_r(1,:)==2),ii) = tau_m_val(ii);
    
        [wage_r,pol_mkt_r_tau(:,:,ii),wage_im,wage_ih,wage_i,pol_mkt_i_tau(:,:,ii),theta_r,theta_i] ...
        = steady_state_cutoff(periods,ini_pop,wage_sm,wage_sh,tau_m_pol(:,:,ii),gamma_Im,gamma_Ih);
    
    pol_mkt_s_tau(:,:,ii) = pol_mkt_r_tau(:,:,ii);
    wage_s_tau = wage_r;
    
    G_ss_tau(:,ii) = sum(tau_m_pol(:,:,ii).*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_tau(:,:,ii)==2) ,2);
    C_ss_tau(:,ii) = sum((1-tau_m_pol(:,:,ii)).*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_tau(:,:,ii)==2) ...
        + wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_tau(:,:,ii)==1) ,2)+  G_ss ;
        
    [pol_mkt_s_opt_tau(:,:,ii),~,~,I_endog_tau(:,ii),~,~,S_endog_tau(:,ii),~,~,R_endog_tau(:,ii),D_endog_tau(:,ii),~,~,~,~,~,~,~,...
      C_endog_tau(:,ii),G_endog_tau(:,ii),Y_endog_tau(:,ii),theta_total_endog_tau(:,ii),error_tau(:,:,ii)] = policy_susceptible(pol_mkt_s_tau(:,:,ii)...
        ,wage_s_tau,sec_weight,sec_elas,occ_weight,occ_elas,tau_m_pol(:,:,ii),wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,...
        pol_mkt_r_tau(:,:,ii),pol_mkt_i_tau(:,:,ii),risk_prop,ini_pop,betta,pi_par,periods);

    % steady-state deviations given optimal responses
    C_dev_endog_tau(:,ii) = (log(C_endog_tau(:,ii))- log(C_ss))*100; %deviation of total consumption
    C_pv_endog_tau(:,ii) = round(sum(C_endog_tau(1:105,ii).*(betta.^(linspace(0,104,104+1)'))));
    
end
C_diff_tau = round(sum( (C_endog_tau(1:105,:) - repmat(C_endog(1:105),1,tau_size))*1000,1));

%% Creating table of results for which occupations choose to leave market
exit_baseline = pol_mkt_s_exo - pol_mkt_s_opt ;
baseline = sum(exit_baseline,1)';
exit_tau_1 = pol_mkt_s_exo - pol_mkt_s_opt_tau(:,:,1);
tax5 = sum(exit_tau_1,1)';
exit_tau_2 = pol_mkt_s_exo - pol_mkt_s_opt_tau(:,:,2);
tax12 = sum(exit_tau_2,1)';
T = table(ID,NOC,NOC_short,dummy_home,Relativehomewage,baseline,tax5,tax12);
tablefile = fullfile(main_dir,'endogenous_occupations_occ_risk.xlsx');
writetable(T,tablefile,'Sheet',1,'Range','A1:H41', 'WriteVariableNames', 1);

%% Figure 10 - Comparing baseline model with endogenous and exogenous model with occupation specific risk
constant_risk = load(main_results); % loading the baseline model results
constant_risk_endog = constant_risk.endog;

ia=2;ib=3;fsize=12;
horz=105;
time=0:1:horz-1;
ticks = 0:35:horz;
newcolors = [0.5 0.5 0.5; 0.25 0.25 0.25; 0 0 0];
figfile = fullfile(fig_dir,'Figure_10.pdf');

figure(10)
subplot(ia,ib,1)
plot(time,constant_risk_endog.I_endog(1:horz)*100,'-',time,I_exo(1:horz)*100,'--',time,I_endog(1:horz)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,constant_risk_endog.S_endog(1:horz)*100,'-',time,S_exo(1:horz)*100,'--',time,S_endog(1:horz)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,constant_risk_endog.R_endog(1:horz)*100,'-',time,R_exo(1:horz)*100,'--',time,R_endog(1:horz)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,constant_risk_endog.D_endog(1:horz)*100,'-',time,D_exo(1:horz)*100,'--',time,D_endog(1:horz)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,constant_risk_endog.theta_total_endog(1:horz)*100,'-',time,theta_total_exo(1:horz)*100,'--',time,theta_total_endog(1:horz)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,constant_risk_endog.C_dev_endog(1:horz),'-',time,C_dev_exo(1:horz),'--',time,C_dev_endog(1:horz),':','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Endogenous - Constant','Exogenous,Occupation-specific', 'Endogenous,Occupation-specific','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(10, 'Position', get(0, 'Screensize'));
orient landscape
print('-f10','-dpdf','-fillpage',figfile);

close all