%% epi-mmb variables 
% Figures are suppressed in the main file 
main;
%{
% Consumption= 
% Labor= 
% technically there are consumption and labor, but only for specific group
% and I have not yet come around to aggregate them 
Output = S_fig.gdp; 
%Susceptible
Infected = frac_age(1,1)*S_fig.M_i_all(1,:)+frac_age(1,2)*S_fig.M_i_all(2,:);
Recovered = frac_age(1,1)*S_fig.M_r(1,:)+frac_age(1,2)*S_fig.M_r(2,:);
Deaths = frac_age(1,1)*S_fig.M_d(1,:)+frac_age(1,2)*S_fig.M_d(2,:);
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');
%}