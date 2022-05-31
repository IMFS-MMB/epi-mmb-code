%clear
SIR_macro_decentralized_US
%sir_macro_decentralized_non_US
Consumption= ((aggC-crss)/crss);
Labour= ((aggH-nrss)/nrss) ;
Output= ((aggC-crss)/crss);
Susceptibles= S(1:size(Consumption));
Infected= I(1:size(Consumption));
Recovered= R(1:size(Consumption));
Deaths = D(1:size(Consumption));

Consumption_ss = 0;
Labour_ss = 0;
Output_ss = 0;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');
