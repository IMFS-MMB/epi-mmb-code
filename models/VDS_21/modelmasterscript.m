dynare VDS_21;
load simulated_results_base;
name = 'Consumption';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Consumption = z;
Consumption_ss = zss;
name = 'Labour';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Labour = z;
Labour_ss = zss;
name = 'Output';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Output = z;
Output_ss = zss; 
name = 'Interest';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Interest = z;
Interest_ss = zss; 
name = 'Inflation';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Inflation = z;
Inflation_ss = zss; 
name = 'In';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Investment = z;
Investment_ss = zss; 


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

