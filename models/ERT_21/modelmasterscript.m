sir_macro;

Consumption = aggC;
Labour = aggH;
Output = A.*aggH;
Susceptibles = S(1:HH); 
Infected = I(1:HH);
Recovered = R(1:HH);
Deaths = D(1:HH);

Consumption_ss = crss;
Labour_ss = nrss;
Output_ss = crss;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');
