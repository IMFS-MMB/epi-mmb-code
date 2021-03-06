%close all
%clearvars
%clc


saveflag    = 1;    % 1 to save results for each run
helper=load('inf_ini.mat');
I0in=helper.helper;
main_uk

% save and present only case when psi = 0.5 (hh underestimate the severity
% of the pandemic)
%Consumption = (C_v1 - par_v1.c_ss) ./ par_v1.c_ss;
%Labour = (H_v1 - par_v1.n_ss) ./ par_v1.n_ss;
Output = day2week(uk_gdp);
Consumption = Output;
Labour = nan(size(Output));
Susceptibles = nan(size(Output));
Infected = day2week(uk_sir) ;
Recovered = nan(size(Output));
Deaths = day2week(uk_death) ;
Inflation= NaN(size(Susceptibles));
Interest= NaN(size(Susceptibles));
Investment = NaN(size(Susceptibles));



Consumption_ss = 0;%Consumption(1);
Labour_ss = Labour(1);
Output_ss = 0;%Output(1);
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;

Inflation_ss= 0;
Interest_ss= 0;
Investment_ss = 0;
save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');%,'Exposed'