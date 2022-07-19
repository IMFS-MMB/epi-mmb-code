%close all
%clear
%clc

% master;
% cd(DIR_SIMULATIONS);
T = readtable('sim_time_series.xlsx');

% Optimal allcation of economy with lockdown "commitment" 
opts = detectImportOptions('sim_time_series.xlsx');
opts.SelectedVariableNames = [4 6 7 9 10 19 22 25 28];
[lockC Y0 YC X0 XC SC IC RC DC] = readvars('sim_time_series.xlsx',opts);

n = size(YC);
Consumption_ss = Y0(1);     % "GDP with no pandemic"
Labour_ss = 1;              % effective labour with no pandemic
Output_ss = Y0(1);          
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
Interest_ss = 0;
Inflation_ss =0;
Investment_ss = X0(1);

Consumption = YC;        
Output = YC;            
Susceptibles = SC;
Infected = IC;
Recovered = RC;
Deaths = DC;
Labour = (ones(n(1),1) - lockC) .* (SC+ 0.5 * IC + RC) ; % definition of effective labour
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = XC;

% Exact replications of figures presented in the paper
%Consumption = YC/Y0(1);        % relative to GDP without pandemic 
%Output = YC/Y0(1);             % relative to GDP without pandemic 
%Susceptibles = SC;
%Infected = IC;
%Recovered = RC;
%Deaths = DC;
%Lockdown = lockC;              % "lockdown" share

%Consumption_ss = 0;     
%Labour_ss = 0;              
%Output_ss = 0;          
%Susceptibles_ss= 0;
%Infected_ss= 0;
%Recovered_ss= 0;
%Deaths_ss= 0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');

