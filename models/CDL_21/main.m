% Main file to run baseline model and other robustness and lockdown
% exercises
% Generates Figures 4, 5, 6, 7, 13, 14, 15, 17, 18 Model
%clc;
%clear all;

%fig_dir = fileparts('C:\Users\clindema\Desktop\CDL_21\Codes\Figures');
% Specify directory where figures will get saved
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

[ID,NOC,NOC_short,employment,~,~,market_wage,dummy_home,rel_wage,home_wage,~,Relativehomewage] = readvars('wage_schedule.xlsx');
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

opts_fsolve=optimoptions('fsolve','Display','off','TolFun',1e-9); %options for fsolve

betta=0.96^(1/52);  %Weekly household discount factor
pid = (7*0.005)/18; % Weekly probability of dying
pir=7/18-pid;     %Weekly probability of recovering
%epsilon = 0.001;  %initial infected
helper=load('inf_ini.mat');
epsilon=helper.helper;

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

% deriving path of infections, consumption given NO endogeneous response
pol_mkt_s_exo = pol_mkt_r;
wage_s_exo = wage_r;
[~,~,~,I_exo,~,~,S_exo,~,~,R_exo,D_exo,~,~,~,~,~,~,~,C_exo,G_exo,Y_exo,theta_total_exo] ...
    = sir_dynamics(pol_mkt_s_exo,wage_s_exo,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,...
    pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given NO endogeneous responses
C_dev_exo = (log(C_exo) - log(C_ss))*100;
C_pv_exo = round(sum(C_exo(1:105,1).*(betta.^(linspace(0,104,104+1)'))));
%Storing in structure
exo.I_exo = I_exo;
exo.S_exo = S_exo;
exo.R_exo = R_exo;
exo.D_exo = D_exo;
exo.theta_total_exo = theta_total_exo;
exo.C_exo = C_exo;
exo.C_dev_exo = C_dev_exo;
exo.C_pv_exo = C_pv_exo;
save('results.mat','exo');

%% Baseline results

pol_mkt_s_init = pol_mkt_r;
wage_s_init = wage_r;

[pol_mkt_s_opt,~,RnotSIR_endog,I_endog,~,~,S_endog,~,~,R_endog,D_endog,~,~,~,~,Cs_endog,Ci_endog,Cr_endog,C_endog,G_endog,Y_endog,theta_total_endog,error_endog] ...
    = policy_susceptible(pol_mkt_s_init,wage_s_init,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog = (log(C_endog)- log(C_ss))*100; %deviation of total consumption
G_dev_endog = zeros(periods+1,1);
C_pv_endog = round(sum(C_endog(1:105,1).*(betta.^(linspace(0,104,104+1)'))));
C_diff = round(sum((C_endog(1:105)-C_exo(1:105))*1000));
%Storing in structure
endog.I_endog = I_endog;
endog.S_endog = S_endog;
endog.R_endog = R_endog;
endog.D_endog = D_endog;
endog.theta_total_endog = theta_total_endog;
endog.C_dev_endog = C_dev_endog;
endog.C_pv_endog = C_pv_endog;
endog.pol_mkt_s_opt = pol_mkt_s_opt;
save('results.mat','endog','-append');

%% Calculating consumption losses across income quantiles
wage = [wage_sh(1,1:4),wage_sm(1,5:end)];
% finding the quantile values and occupations in each quantile
quantile_vals = quantile(wage(1,:),[0.25,0.5,0.75]);
loc_25 = find(wage(1,:)<=quantile_vals(1));
loc_50 = find(wage(1,:)>quantile_vals(1) & wage(1,:)<=quantile_vals(2));
loc_75 = find(wage(1,:)>quantile_vals(2) & wage(1,:)<=quantile_vals(3));
loc_100 = find(wage(1,:)>quantile_vals(3));
% consumption at steady state for different quantiles
C_ss_25 = sum(wage(1,loc_25).*repmat(ini_pop(1,loc_25),periods+1,1),2);
C_ss_50 = sum(wage(1,loc_50).*repmat(ini_pop(1,loc_50),periods+1,1),2);
C_ss_75 = sum(wage(1,loc_75).*repmat(ini_pop(1,loc_75),periods+1,1),2);
C_ss_100 = sum(wage(1,loc_100).*repmat(ini_pop(1,loc_100),periods+1,1),2);
% consumption in baseline endogenous model for different quantiles
C_endog_25 = sum( Cs_endog(:,loc_25) + Ci_endog(:,loc_25) + Cr_endog(:,loc_25),2);
C_endog_50 = sum( Cs_endog(:,loc_50) + Ci_endog(:,loc_50) + Cr_endog(:,loc_50),2);
C_endog_75 = sum( Cs_endog(:,loc_75) + Ci_endog(:,loc_75) + Cr_endog(:,loc_75),2);
C_endog_100 = sum( Cs_endog(:,loc_100) + Ci_endog(:,loc_100) + Cr_endog(:,loc_100),2);
% percentage deviation from steady state
C_dev_endog_25 = (log(C_endog_25) - log(C_ss_25))*100;
C_dev_endog_50 = (log(C_endog_50) - log(C_ss_50))*100;
C_dev_endog_75 = (log(C_endog_75) - log(C_ss_75))*100;
C_dev_endog_100 = (log(C_endog_100) - log(C_ss_100))*100;

%% Role of Risk aversion (Figure 13)
pol_mkt_s_linear = pol_mkt_r;
wage_s_linear = wage_r;

[~,~,RnotSIR_endog,I_endog_linear,~,~,S_endog_linear,~,~,R_endog_linear,D_endog_linear,~,~,~,~,~,~,~,C_endog_linear,~,~,theta_total_endog_linear,~] ...
    = policy_susceptible_linear(pol_mkt_s_linear,wage_s_linear,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog_linear = (log(C_endog_linear)- log(C_ss))*100; %deviation of total consumption

%% Lower probabilitiy of infections (Figure 7)
% 10 percent reduction in all probabilities of infection
pi_par(2) = 0.9*pi;
pi_par(3) = 0.9*varepsilon;
pi_par(5) = 0.9*pic;

pol_mkt_s_low = pol_mkt_r;
wage_s_low = wage_r;

[~,~,~,I_endog_low,~,~,S_endog_low,~,~,R_endog_low,D_endog_low,~,~,~,~,~,~,~,C_endog_low,~,~,theta_total_endog_low,~] ...
    = policy_susceptible(pol_mkt_s_low,wage_s_low,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

C_dev_endog_low = (log(C_endog_low)- log(C_ss))*100; %deviation of total consumption

%% Saving Infection data
% This will be needed to generate Figure 8
Week = [1:1:251]';
T = table(Week,I_endog,S_endog,R_endog,D_endog,C_endog,theta_total_endog, ...
    I_endog_low,S_endog_low,R_endog_low,D_endog_low,C_endog_low,theta_total_endog_low);
writetable(T,'endogenous_infections.xlsx','Sheet',1,'Range','A1:M252', 'WriteVariableNames', 1);

%% Tax and Rebate targeted to some
tau = 0.05; %proportional tax rate 
pi_par(1) = pir;
pi_par(2) = pi;
pi_par(3) = varepsilon;
pi_par(4) = epsilon;
pi_par(5) = pic;
pi_par(6) = pid;

% Calculating occupations working at market who will be get government
% rebates
med_val = quantile(wage_sh(1,5:end),0.5);
below_med = find(wage_sh(1,5:end)<=med_val) + 4;
above_med = find(wage_sh(1,5:end)>med_val) + 4;

pol_mkt_s_rebate = pol_mkt_r;
wage_s_rebate = wage_r;

[~,~,RnotSIR_endog,I_endog_rebate,~,~,S_endog_rebate,~,~,R_endog_rebate,D_endog_rebate,~,~,~,~,~,~,~,C_endog_rebate,~,~,theta_total_endog_rebate,~] ...
    = policy_susceptible_rebate(pol_mkt_s_rebate,wage_s_rebate,sec_weight,sec_elas,occ_weight,occ_elas,tau,below_med,above_med,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog_rebate = (log(C_endog_rebate)- log(C_ss))*100; %deviation of total consumption

%% Vaccine and Treatment - Figure 15
delta_val = [1/26; 1/13]; % two possible values for probabilitiy that vaccine 
% is discovered
delta_size = size(delta_val,1);
% generating matrices where results will be stored
pol_mkt_s_opt_delta = NaN*ones(periods+1,n_occ,delta_size);
theta_total_endog_delta = NaN*ones(periods+1,delta_size);
I_endog_delta = NaN*ones(periods+1,delta_size);
S_endog_delta = NaN*ones(periods+1,delta_size);
R_endog_delta = NaN*ones(periods+1,delta_size);
D_endog_delta = NaN*ones(periods+1,delta_size);
C_endog_delta = NaN*ones(periods+1,delta_size);
C_dev_endog_delta = NaN*ones(periods+1,delta_size);

for ii = 1:2

delta_v = delta_val(ii); % Weekly probability of vaccines
delta_t = delta_val(ii); % Weekly probability of treatment

% storing infection parameters
pi_par(1) = pir;
pi_par(2) = pi;
pi_par(3) = varepsilon;
pi_par(4) = epsilon;
pi_par(5) = pic;
pi_par(6) = pid;
pi_par(7) = delta_v;
pi_par(8) = delta_t;

% Endogenous results with vaccine and treatment

pol_mkt_s_v = pol_mkt_r;
wage_s_v = wage_r;

[pol_mkt_s_opt_delta(:,:,ii),~,~,I_endog_delta(:,ii),~,~,S_endog_delta(:,ii),~,~,R_endog_delta(:,ii),D_endog_delta(:,ii),~,~,~,~,~,~,~,C_endog_delta(:,ii),~,~,theta_total_endog_delta(:,ii),~] ...
    = policy_susceptible_vaccine(pol_mkt_s_v,wage_s_v,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

% steady-state deviations given optimal responses
C_dev_endog_delta(:,ii) = (log(C_endog_delta(:,ii))- log(C_ss))*100; %deviation of total consumption

end

%% Changing gamma or productivity of infected (Figure 14)

pi_par(1) = pir;
pi_par(2) = pi;
pi_par(3) = varepsilon;
pi_par(4) = epsilon;
pi_par(5) = pic;
pi_par(6) = pid;
% Two possible values that productivity of infected will become
gamma_h_val = [0.6;0.9];
gamma_m_val = [0.6;0.9];
gamma_size = size(gamma_h_val,1);
%pre-allocate endogenous and exogenous output
C_ss_gamma = NaN*ones(periods+1,gamma_size);
G_ss_gamma = NaN*ones(periods+1,gamma_size);

pol_mkt_s_gamma = NaN*ones(periods+1,n_occ,gamma_size);
pol_mkt_r_gamma = NaN*ones(periods+1,n_occ,gamma_size);
pol_mkt_i_gamma = NaN*ones(periods+1,n_occ,gamma_size);

pol_mkt_s_opt_gamma = NaN*ones(periods+1,n_occ,gamma_size);
theta_opt_gamma = NaN*ones(periods+1,gamma_size);
theta_total_endog_gamma = NaN*ones(periods+1,gamma_size);
I_endog_gamma = NaN*ones(periods+1,gamma_size);
S_endog_gamma = NaN*ones(periods+1,gamma_size);
R_endog_gamma = NaN*ones(periods+1,gamma_size);
D_endog_gamma = NaN*ones(periods+1,gamma_size);
C_endog_gamma = NaN*ones(periods+1,gamma_size);
G_endog_gamma = NaN*ones(periods+1,gamma_size);
Y_endog_gamma = NaN*ones(periods+1,gamma_size);
C_dev_endog_gamma = NaN*ones(periods+1,gamma_size);
C_pv_endog_gamma = NaN*ones(1,gamma_size);

for ii = 1:gamma_size

    gamma_Ih = gamma_h_val(ii); %relative productivity of infected to susceptible for working at home
    gamma_Im = gamma_m_val(ii); %relative productivity of infected to susceptible for working at market

    [wage_r,pol_mkt_r_gamma(:,:,ii),wage_im,wage_ih,wage_i,pol_mkt_i_gamma(:,:,ii),theta_r,theta_i] ...
        = steady_state_cutoff(periods,ini_pop,wage_sm,wage_sh,tau_m,gamma_Im,gamma_Ih);
    
    pol_mkt_s_gamma(:,:,ii) = pol_mkt_r_gamma(:,:,ii);
    wage_s_gamma = wage_r;
    
    G_ss_gamma(:,ii) = sum(tau_m.*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_gamma(:,:,ii)==2) ,2);
    C_ss_gamma(:,ii) = sum((1-tau_m).*wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_gamma(:,:,ii)==2) ...
        + wage_r.*repmat(ini_pop,periods+1,1).*(pol_mkt_r_gamma(:,:,ii)==1) ,2)+  G_ss ;
            
    [pol_mkt_s_opt_gamma(:,:,ii),~,~,I_endog_gamma(:,ii),~,~,S_endog_gamma(:,ii),~,~,R_endog_gamma(:,ii),D_endog_gamma(:,ii),~,~,~,~,~,~,~,...
      C_endog_gamma(:,ii),G_endog_gamma(:,ii),Y_endog_gamma(:,ii),theta_total_endog_gamma(:,ii),error_gamma(:,:,ii)] = policy_susceptible(pol_mkt_s_gamma(:,:,ii)...
        ,wage_s_gamma,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,...
        pol_mkt_r_gamma(:,:,ii),pol_mkt_i_gamma(:,:,ii),risk_prop,ini_pop,betta,pi_par,periods);

    % steady-state deviations given optimal responses
    C_dev_endog_gamma(:,ii) = (log(C_endog_gamma(:,ii))- log(C_ss))*100; %deviation of total consumption
    C_pv_endog_gamma(:,ii) = round(sum(C_endog_gamma(:,ii).*(betta.^(linspace(0,periods,periods+1)'))));

end
C_diff_gamma = round(sum( (C_endog_gamma(1:105,:) - repmat(C_endog(1:105),1,2))*1000,1));

%% Changing tau's to create various lockdown scenarios (Figure 5,6)
gamma_Ih = 0.8; %relative productivity of infected to susceptible for working at home
gamma_Im = 0.8; %relative productivity of infected to susceptible for working at market

% three lockdown scenarios
tau_m_val = [0.05; 0.12; 0.20]; % wedge on market wage
dur_tax = [36;12;52]; % duration of lockdown
tau_size = size(tau_m_val,1);
tau_m_pol = zeros(periods+1,n_occ,tau_size);

%pre-allocate endogenous and exogenous output
C_ss_tau = NaN*ones(periods+1,tau_size);
G_ss_tau = NaN*ones(periods+1,tau_size);

pol_mkt_s_tau = NaN*ones(periods+1,n_occ,tau_size);
pol_mkt_r_tau = NaN*ones(periods+1,n_occ,tau_size);
pol_mkt_i_tau = NaN*ones(periods+1,n_occ,tau_size);

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
writetable(T,'endogenous_occupations.xlsx','Sheet',1,'Range','A1:H41', 'WriteVariableNames', 1);

Susceptibles=S_endog;
Infected=I_endog;
Recovered=R_endog;
Deaths=D_endog;
Labour=theta_total_endog;
Consumption=C_dev_endog;
Output=Y_endog;

Consumption_ss=0;%C_ss
Output_ss=Y_ss(1);
Labour_ss=theta_ss;
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');


%% Figures
%{
ia=2;ib=3;fsize=12;
horz=105;
time=0:1:horz-1;
ticks = 0:35:horz;

% Figure 4: baseline model exogenous versus endogenous response
str = {['Consumption gap is ' , num2str(C_diff)]};
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_4.pdf');

figure(4);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_exo(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_exo(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_exo(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_exo(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz)*100,'-',time,theta_total_exo(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_exo(1:horz),'--','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);
annotation('textbox', [0.77, 0.12, 0.1, 0.1], 'String',str,'EdgeColor','w');

leg1=legend('Endogeneous','Exogenous','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(4, 'Position', get(0, 'Screensize'));
orient landscape
print('-f4','-dpdf','-fillpage',figfile);
%% Lockdown Figures

newcolors = [0.5 0.5 0.5; 0.25 0.25 0.25; 0 0 0];

% Figure 5 - Comparing no lockdown with mild but longer and strong but short
% lockdowns
str = {['Consumption gap (\tau^m = ' num2str(tau_m_val(1)) ') is ' , num2str(C_diff_tau(1,1)) ] ...
    ['Consumption gap (\tau^m = ' num2str(tau_m_val(2)) ') is ' , num2str(C_diff_tau(1,2)) ]};
figfile = fullfile(fig_dir,'Figure_5.pdf');

figure(5);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_tau(1:horz,1)*100,'--',time,I_endog_tau(1:horz,2)*100,'-.','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_tau(1:horz,1)*100,'--',time,S_endog_tau(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_tau(1:horz,1)*100,'--',time,R_endog_tau(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_tau(1:horz,1)*100,'--',time,D_endog_tau(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz,1)*100,'-',time,theta_total_endog_tau(1:horz,1)*100,'--',...
    time,theta_total_endog_tau(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_tau(1:horz,1),'--',time,C_dev_endog_tau(1:horz,2),':','LineWidth',2);
box off;
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Baseline Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
colororder(newcolors);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);
annotation('textbox', [0.77, 0.1, 0.1, 0.1], 'String',str,'EdgeColor','w');

leg1=legend('Baseline',['\tau^{m} = ' num2str(tau_m_val(1)) ', duration ' num2str(dur_tax(1)) ' weeks'], ...
    ['\tau^{m} = ' num2str(tau_m_val(2)) ', duration ' num2str(dur_tax(2)) ' weeks'],'FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.179142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(5, 'Position', get(0, 'Screensize'));
orient landscape
print('-f5','-dpdf','-fillpage',figfile);

% Figure 6 - Comparing no lockdown with mild 1 year lockdown
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_6.pdf');

figure(6);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_tau(1:horz,3)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_tau(1:horz,3)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_tau(1:horz,3)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_tau(1:horz,3)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz,1)*100,'-',time,theta_total_endog_tau(1:horz,3)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);


subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_tau(1:horz,3),'--','LineWidth',2);
box off;
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Baseline Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
colororder(newcolors);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline',['\tau^{m} = ' num2str(tau_m_val(ii)) ', duration ' num2str(dur_tax(ii)) ' weeks'],'FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.179142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(6, 'Position', get(0, 'Screensize'));
orient landscape
print('-f6','-dpdf','-fillpage',figfile);

% Figure 7: Comparing baseline model with 10% reduction in all
% probabilities of infections
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_7.pdf');

figure(7);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_low(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_low(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_low(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_low(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz)*100,'-',time,theta_total_endog_low(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_low(1:horz),'--','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline','10% reduction in probability of infections','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(7, 'Position', get(0, 'Screensize'));
orient landscape
print('-f7','-dpdf','-fillpage',figfile);

% Figure 13: Comparing baseline model with model with risk neutral
% preferences
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_13.pdf');

figure(13);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_linear(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_linear(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_linear(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_linear(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz)*100,'-',time,theta_total_endog_linear(1:horz)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_linear(1:horz),'--','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline (Log utility)','Linear Utility','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(13, 'Position', get(0, 'Screensize'));
orient landscape
print('-f13','-dpdf','-fillpage',figfile);

% Figure 14 - Comparing baseline model with lower and higher productivity
% of infected
str = {['Consumption gap (\gamma_m = ' num2str(gamma_m_val(1)) ') is ' , num2str(C_diff_gamma(1,1)) ] ...
    ['Consumption gap (\gamma_m = ' num2str(gamma_m_val(2)) ') is ' , num2str(C_diff_gamma(1,2)) ]};
newcolors = [0.5 0.5 0.5; 0.25 0.25 0.25; 0 0 0];
figfile = fullfile(fig_dir,'Figure_14.pdf');

figure(14);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_gamma(1:horz,1)*100,'--',time,I_endog_gamma(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_gamma(1:horz,1)*100,'--',time,S_endog_gamma(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_gamma(1:horz,1)*100,'--',time,R_endog_gamma(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_gamma(1:horz,1)*100,'--',time,D_endog_gamma(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz,1)*100,'-',time,theta_total_endog_gamma(1:horz,1)*100,'--',...
    time,theta_total_endog_gamma(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);


subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_gamma(1:horz,1),'--',time,C_dev_endog_gamma(1:horz,2),':','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Baseline Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);
annotation('textbox', [0.76, 0.11, 0.1, 0.1], 'String',str,'EdgeColor','w');

leg1=legend('Baseline',['\gamma_I^{m} = \gamma_I^{h} = ' num2str(gamma_m_val(1))], ...
    ['\gamma_I^{m} = \gamma_I^{h} = ' num2str(gamma_m_val(2))],'FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.1979142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(14, 'Position', get(0, 'Screensize'));
orient landscape
print('-f14','-dpdf','-fillpage',figfile);

% Figure 15 - Comparing baseline model with introduction of vaccines and
% treatment
newcolors = [0.5 0.5 0.5; 0.25 0.25 0.25; 0 0 0];
figfile = fullfile(fig_dir,'Figure_15.pdf');

figure(15);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_delta(1:horz,1)*100,'--',time,I_endog_delta(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_delta(1:horz,1)*100,'--',time,S_endog_delta(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_delta(1:horz,1)*100,'--',time,R_endog_delta(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_delta(1:horz,1)*100,'--',time,D_endog_delta(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz,1)*100,'-',time,theta_total_endog_delta(1:horz,1)*100,'--',...
    time,theta_total_endog_delta(1:horz,2)*100,':','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);


subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_delta(1:horz,1),'--',time,C_dev_endog_delta(1:horz,2),':','LineWidth',2);
box off;
colororder(newcolors);
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Baseline Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline','\delta_{v} = \delta_{c} = 1/26' ,'\delta_{v} = \delta_{c} = 1/13','FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(15, 'Position', get(0, 'Screensize'));
orient landscape
print('-f15','-dpdf','-fillpage',figfile);

% Figure 17 - Comparing baseline model with model where government taxes
% all individuals but rebates poorer individuals
newcolors = [0.5 0.5 0.5; 0 0 0];
figfile = fullfile(fig_dir,'Figure_17.pdf');

figure(17);
subplot(ia,ib,1)
plot(time,I_endog(1:horz)*100,'-',time,I_endog_rebate(1:horz,1)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Infected, I','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,2)
plot(time,S_endog(1:horz)*100,'-',time,S_endog_rebate(1:horz,1)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Susceptibles, S','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,3)
plot(time,R_endog(1:horz)*100,'-',time,R_endog_rebate(1:horz,1)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Recovered, R','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,4)
plot(time,D_endog(1:horz)*100,'-',time,D_endog_rebate(1:horz,1)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Deaths, D','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

subplot(ia,ib,5)
plot(time,theta_total_endog(1:horz,1)*100,'-',time,theta_total_endog_rebate(1:horz,1)*100,'--','LineWidth',2);
box off;
colororder(newcolors);
title('Employment Market','FontSize',fsize);
ylabel('% Share of Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);


subplot(ia,ib,6)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_rebate(1:horz,1),'--','LineWidth',2);
box off;
title('Consumption','FontSize',fsize);
ylabel('% Dev. from Baseline Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
colororder(newcolors);
set(gca,'FontSize',fsize);
xticks(ticks);
xlim([0 horz]);

leg1=legend('Baseline',[num2str(tau*100) '% Proportional tax, rebate targeted'],'FontSize',16);
set(leg1,...
    'Position',[-0.0704335567178815 0.94078940309976 1.19107142857143 0.074047619047619],...
    'Orientation','horizontal');
legend boxoff;
set(17, 'Position', get(0, 'Screensize'));
orient landscape
print('-f17','-dpdf','-fillpage',figfile);

% Figure 18 - Consumption quartiles in model
horz = 105;
time=0:1:horz-1;
ticks = 0:35:horz;
fsize = 16;
newcolors = [0.5 0.5 0.5; 0 0 0; 0.1 0.1 0.1; 0.35 0.35 0.35; 0.2 0.2 0.2];
figfile = fullfile(fig_dir,'Figure_18_b.pdf');

figure(18)
plot(time,C_dev_endog(1:horz),'-',time,C_dev_endog_25(1:horz),'--',time,C_dev_endog_50(1:horz),':',time,C_dev_endog_75(1:horz),'-*' ...
    ,time,C_dev_endog_100(1:horz),'-o','LineWidth',2);
box off;
ylabel('% Deviation from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
leg1 = legend('Overall','Below 25th','Between 25th and 50th','Between 50th and 75th','Above 75th','Location','bestoutside','Orientation','horizontal');
xticks(ticks);
xlim([0 horz]);
legend boxoff;
colororder(newcolors);
set(leg1,'FontSize',fsize);
set(18, 'Position', get(0, 'Screensize'));
orient landscape
print('-f18','-dpdf','-fillpage',figfile);

close all
%}