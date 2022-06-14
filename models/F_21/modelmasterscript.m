%% modelmasterscript
%%  running the baseline simulations 
%cd main
addpath(horzcat(pwd,'\main'))
calibrate_model
multipliers_normal_times
%crisis_experiments

%% Save the variables

Consumption_quarterly= (pchi)*exp(irf_store.nopol(strmatch('Cb', M_.endo_names, 'exact'),1:end))+(1-pchi)*exp(irf_store.nopol(strmatch('Cs', M_.endo_names, 'exact'),1:end));
Consumption= repelem((pchi)*exp(irf_store.nopol(strmatch('Cb', M_.endo_names, 'exact'),1:end))+(1-pchi)*exp(irf_store.nopol(strmatch('Cs', M_.endo_names, 'exact'),1:end)),12); %convert to weekly from quarterly
Labour_quarterly = exp(N_a) + exp(N_n); 
Labour = repelem(exp(N_a) + exp(N_n),12); 
Output_quarterly = exp(irf_store.nopol(strmatch('GDP', M_.endo_names, 'exact'),1:end)); 
Output = repelem(exp(irf_store.nopol(strmatch('GDP', M_.endo_names, 'exact'),1:end)),12); 

Consumption_ss = pchi*detss.Cb+(1-pchi)*detss.Cs; 
Labour_ss = detss.N_a + detss.N_n; 
Output_ss= detss.GDP; 

%% as % deviation from Steady State %%
% Consumption_dev_ss=100*(Consumption/Consumption_ss -1); 
% Labour_dev_ss=100*(Labour/Labour_ss -1); 
% Output_dev_ss=100*(Output/Output_ss -1); 

save('simulated_results.mat','Consumption','Labour','Output');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss');

