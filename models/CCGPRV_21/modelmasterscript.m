%%% This replicates the forced savings shock, for other shocks check det_simul_script_irfs2.dyn
% dynare 4.6.4 is needed 

dynare gemc console nointeractive
%clear
currentFolder = pwd; 
global M_ options_ oo_
%cd gemc\Output\

load gemc_simul_DUMMY_TUC.mat
name='LYOBS_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Output_q=(z(2:end)-zss); 
Output=repelem(Output_q,12);
Output_ss=0; 

name='LC_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Consumption_q=(z(2:end)-zss); 
Consumption=repelem(Consumption_q,12);
Consumption_ss=0; 

name='LN_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Labour_q=(z(2:end)-zss); 
Labour=repelem(Labour_q,12);
Labour_ss=0; 

name='PHICPI_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
%Inflation_q=(z(2:end)-zss); 
Inflation_q=(exp(z(2:end))-exp(zss)); 
Inflation=repelem(Inflation_q,12);
Inflation_ss=0; 

name='LI_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
Investment_q=(z(2:end)-zss); 
Investment=repelem(Investment_q,12);
Investment_ss=0;

name='INOM_EA';
[z,zss]=dyn2vec(M_, oo_, options_,name);
%Interest_q=4*(z(2:end)-zss); 
Interest_q=(exp(z(2:end))-exp(zss)); 
Interest=repelem(Interest_q,12);
Interest_ss=0;

% or 
% name='LYOBS_EA';
% [z,zss]=dyn2vec(M_, oo_, options_,name);
% Output_q=z(2:end); 
% Output=repelem(Output_q,12);
% Output_ss_q=zss; 
% Output_ss=repelem(Output_ss_q,12);
% 
% name='LC_EA';
% [z,zss]=dyn2vec(M_, oo_, options_,name);
% Consumption_q=z(2:end); 
% Consumption=repelem(Consumption_q,12);
% Consumption_ss_q=zss; 
% Consumption_ss=repelem(Consumption_ss_q,12);
% Consumptiontest=100*(Consumption_q(1:end)-Consumption_ss_q);
% 
% name='LN_EA';
% [z,zss]=dyn2vec(M_, oo_, options_,name);
% Labor_q=z(2:end); 
% Labor=repelem(Labor_q,12);
% Labor_ss_q=zss; 
% Labor_ss=repelem(Labor_ss_q,12);



Susceptibles = nan(size(Consumption));
Deaths = nan(size(Consumption));
Infected = nan(size(Consumption));
Recovered =  nan(size(Consumption));

Deaths_ss = nan;
Infected_ss= nan;
Recovered_ss= nan;
Susceptibles_ss= nan;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');
