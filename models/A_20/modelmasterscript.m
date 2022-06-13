%clear
i= 3; % Choose transmission of the disease, 7 options, R0 = [3.0,2.8,2.5,2.2,2.0,1.8,1.6];
basecasesimulations

Susceptibles= day2week(y(:,1));
Infected= day2week(y(:,3));
%Exposed= day2week(y(:,2));
Recovered= day2week(y(:,4));
Deaths = NaN(size(Susceptibles));
Consumption = NaN(size(Susceptibles));
Labour = NaN(size(Susceptibles));
Output = NaN(size(Susceptibles));


Consumption_ss = 0;
Labour_ss = 0;
Output_ss = 0;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');%,'Exposed'
