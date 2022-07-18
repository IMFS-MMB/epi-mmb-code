dynare VDS_21;
load simulated_results_base;
name = 'Consumption';
[z]=dyn2vec(M_, oo_, options_,name);
Consumption = z;
Consumption_ss = z(1);
name = 'Labour';
[z]=dyn2vec(M_, oo_, options_,name);
Labour = z;
Labour_ss = z(1);
name = 'Output';
[z]=dyn2vec(M_, oo_, options_,name);
Output = z;
Output_ss = z(1); 
name = 'Interest';
[z]=dyn2vec(M_, oo_, options_,name);
Interest = z;
Interest_ss = z(1); 
name = 'Inflation';
[z]=dyn2vec(M_, oo_, options_,name);
Inflation = z;
Inflation_ss = z(1); 
name = 'Investment';
[z]=dyn2vec(M_, oo_, options_,name);
Investment = z;
Investment_ss = z(1); 


load SIR;
Susceptibles = Y_SIR_eu(:,1);
Infected = Y_SIR_eu(:,2); 
Recovered =  Y_SIR_eu(:,3);
Deaths = NaN(length(Susceptibles),1);
Deaths_ss = 0;
Susceptibles_ss = 0;
Infected_ss = 0;
Recovered_ss = 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');

