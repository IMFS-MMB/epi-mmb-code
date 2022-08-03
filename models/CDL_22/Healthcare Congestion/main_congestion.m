% File to create Figure 16 in paper where infections create congestion in
% healthcare

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

risk_prop = -1.*(pol_mkt_r(1,:)==1); % Baseline model has constant risk of transmitting infection 
% through market activity
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
kappa = 1.585; % congestion parameter

go_calibrate_pis;

% storing infection parameters
pi_par(1) = pir;
pi_par(2) = pi;
pi_par(3) = varepsilon;
pi_par(4) = epsilon;
pi_par(5) = pic;
pi_par(6) = pid;
pi_par(7) = kappa;

pol_mkt_s_init = pol_mkt_r;
wage_s_init = wage_r;

[pol_mkt_s_opt,~,RnotSIR_endog,I_endog,~,~,S_endog,~,~,R_endog,D_endog,~,~,~,~,Cs_endog,Ci_endog,Cr_endog,C_endog,G_endog,Y_endog,theta_total_endog,error_endog] ...
    = policy_susceptible(pol_mkt_s_init,wage_s_init,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog = (log(C_endog)- log(C_ss))*100; %deviation of total consumption

%% Figure 16 - Comparing baseline model and model with congestion
baseline = load(main_results); % loading the baseline model results
baseline_endog = baseline.endog;

ia=2;ib=3;fsize=12;
horz=105;
time=0:1:horz-1;
ticks = 0:35:horz;
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_16.pdf');

figure(16);
subplot(ia,ib,1)
plot(time,baseline_endog.I_endog(1:horz)*100,'-',time,I_endog(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,baseline_endog.S_endog(1:horz)*100,'-',time,S_endog(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,baseline_endog.R_endog(1:horz)*100,'-',time,R_endog(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,baseline_endog.D_endog(1:horz)*100,'-',time,D_endog(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,baseline_endog.theta_total_endog(1:horz)*100,'-',time,theta_total_endog(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,baseline_endog.C_dev_endog(1:horz),'-',time,C_dev_endog(1:horz),'--','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline','Healthcare congestion','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(16, 'Position', get(0, 'Screensize'));
orient landscape
print('-f16','-dpdf','-fillpage',figfile);

close all