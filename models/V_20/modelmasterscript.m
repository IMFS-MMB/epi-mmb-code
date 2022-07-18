% competitive equilibrium solution
%close all
%clearvars
%clc

Figure5_6_CompetitiveEq;
save allresults;

Consumption = cs_demand;    % aggregate social goods only
Labour = L1;            
Output = B .* (frac_alive*g');
Susceptibles = S_phi*g';
Infected = I_overall;
Recovered = (1-S_phi-I_phi)*lambda_R/(lambda_R+lambda_D)*g';
Deaths =  1-frac_alive*g';
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = NaN(length(Consumption),1);

Consumption_ss = cs_ss*g';
Labour_ss = 1;
Output_ss = B_ss;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
Interest_ss = 0;
Inflation_ss =0;
Investment_ss=0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');



