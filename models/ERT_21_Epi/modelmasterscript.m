sir_macro;

Consumption = aggC;
Labour = aggH;
Output = A.*aggH;
Susceptibles = S(1:HH); 
Infected = I(1:HH);
Recovered = R(1:HH);
Deaths = D(1:HH);
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = NaN(length(Consumption),1);

Consumption_ss = crss;
Labour_ss = nrss;
Output_ss = crss;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
Interest_ss = 0;
Inflation_ss = 0;
Investment_ss = 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');
