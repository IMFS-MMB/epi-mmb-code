main

Susceptibles=S_endog;
Infected=I_endog;
Recovered=R_endog;
Deaths=D_endog;
Labour=theta_total_endog;
Consumption=C_dev_endog./100;
Output=Y_endog;
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = NaN(length(Consumption),1);

Consumption_ss=0;%C_ss
Output_ss=Y_ss(1);
Labour_ss=theta_ss;
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;
Interest_ss = 0;
Inflation_ss = 0;
Investment_ss = 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');
