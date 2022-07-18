%close all
%clearvars
%clc


saveflag    = 1;    % 1 to save results for each run

% Need to feed "par" containing parameter values and "MUs" containing
% containment policy and wage subsidies.

load HKK_parameters
helper=load('inf_ini.mat');
par_ERT0.I0=helper.helper;
par_ERT0.S0=1-par_ERT0.I0;

par_SIR0.I0=par_ERT0.I0;
par_SIR0.S0=par_ERT0.S0;

par_HKK0.I0=par_ERT0.I0;
par_HKK0.S0=par_ERT0.S0;

PSIs = [0.5:0.5:4];
% PSIs = [0.4:0.3:1.6];

npsi = length(PSIs);

for i = 1:npsi+1
    spec = i;
    
    if i <= npsi
        par = par_ERT0;
        
        par.psi = PSIs(i);
        par.chi = 1;
        par.chi_bar = par.chi;
        par.kappa = 0;
        par.chi_AB = 1;
    elseif  i == npsi+1
        % SIR version
        par = par_SIR0;
        
        par.chi = 1;
        par.chi_bar = par.chi;
        par.kappa = 0;
        par.chi_AB = 1;
    end
    % Solve for the imperfect information model
    imperfect_Solve_Compact
end

% save and present only case when psi = 0.5 (hh underestimate the severity
% of the pandemic)
Consumption = (C_v1 - par_v1.c_ss) ./ par_v1.c_ss;
Labour = (H_v1 - par_v1.n_ss) ./ par_v1.n_ss;
Output = (C_v1 - par_v1.c_ss) ./ par_v1.c_ss;
Susceptibles = S_v1 ;
Infected = I_v1 ;
Recovered = R_v1 ;
Deaths = D_v1 ;

Consumption_ss = 0;
Labour_ss = 0;
Output_ss = 0;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');

