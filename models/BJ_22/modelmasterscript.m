main;

Consumption = NaN(156,1);
Labour = C{1}.n+C{1}.h;
Output = C{1}.production;
Susceptibles = C{1}.susceptible; 
Infected = C{1}.exposed + C{1}.incubation + C{1}.s + C{1}.a;
Recovered = C{1}.r_a + C{1}.r_s;
Deaths = C{1}.d;
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = NaN(length(Consumption),1);

Consumption_ss = 0;
Labour_ss = C{1}.n(1)+C{1}.h(1);
Output_ss = C{1}.production(1);
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
Interest_ss = 0;
Inflation_ss = 0;
Investment_ss = 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');
