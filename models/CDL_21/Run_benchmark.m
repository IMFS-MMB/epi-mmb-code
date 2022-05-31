%% epi-mmb variables 
main; 
%{
Consumption= C_dev_endog; %% Dev. from Initial Steady State'
Labor= theta_total_endog; % Levels
Output= (log(Y_endog)- log(Y_ss))*100; %% Dev. from Initial Steady State'
Susceptible= S_endog*100; % perc Share of Population
Infected= I_endog*100; % perc Share of Population
Recovered= R_endog*100; % perc Share of Population
Deaths= D_endog*100; % perc Share of Population
save('simulated_results.mat','Consumption','Labor','Output','Susceptible','Infected','Recovered','Deaths')
%}