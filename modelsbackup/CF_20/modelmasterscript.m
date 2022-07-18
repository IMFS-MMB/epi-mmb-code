%executed = system('R CMD BATCH C:\Users\jaspe\Documents\GitHub\exe-epi-mmb\models\CF_20\CF_20\Codes\CruciniOFlahertyCovidFiles/SIRMultipleLocations.R')
% clear all;clc;

helper=load('inf_ini.mat');
I0in=helper.helper;

try
   delete  initval.txt
end
save('initval.txt', 'I0in', '-ASCII','-append');

!Rscript SIRMultipleLocations.r
%T = readtable('dtresults.csv','NumHeaderLines',1);

dtresults = readtable('dtresults.csv');
dtparam = readtable('dtparams.csv');

Susceptibles = specificpositionsextract(dtresults.S_Agg);
Infected = specificpositionsextract(dtresults.I_Agg);
Deaths = specificpositionsextract(dtresults.D_Agg);
Recovered = specificpositionsextract(dtresults.R_Agg);

%Consumption = 100*dtresults.C_Agg/dtparam.values(dtparam.name == "c_rss") - 100;
Consumption = specificpositionsextract(dtresults.C_Agg);
%Labour = 100*dtresults.hours_Agg/dtparam.values(dtparam.name == "hours_rss") - 100;
Labour = specificpositionsextract(dtresults.hours_Agg);
Output = 10*Labour;

Consumption_ss = dtparam.values(dtparam.name == "c_rss");
Labour_ss = dtparam.values(dtparam.name == "hours_rss");
Output_ss = 10*dtparam.values(dtparam.name == "hours_rss");

Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');%,'Exposed'
