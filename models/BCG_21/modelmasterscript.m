%The model is solve using occbin where it  linearize everything, all irfs
%except SIRD varibles are expressed as percentage dev. from steady state
%(time 100 to get percent), SIRD in daily, Macro vars in monthly, quarterly IRFs reported the paper figure
%clear 
case_switch = 1;  % 1, baseline, 2 cheapest, 3 shift lockdown to outside the labor force, 4 too strict
call_two_sector_path % the results are provide for 1 sector case and 2 sector, the irfs are named differently, 1=2 sector, 2=1 sector


Consumption =  repelem(s1_c_irf./s1_c_ss,4); % convert to weekly from monthly IRFs
Labour = repelem((s1_l1_irf+s1_l2_irf)./(s1_l1_ss+s1_l2_ss),4); % convert to weekly from monthly IRFs
Output = repelem(s1_y_irf/s1_y_ss,4); % convert to weekly from monthly IRFs
Susceptibles =  day2week(sum(Smat)'); % convert to weekly from daily IRFs
Infected = day2week(sum(Imat)'); % convert to weekly from daily IRFs
Recovered =  day2week(sum(Cmat)'); % convert to weekly from daily IRFs
Deaths = day2week(sum(Dmat)'); % convert to weekly from daily IRFs
Susceptibles(97:end)=[];
Infected(97:end)=[];
Recovered(97:end)=[];
Deaths(97:end)=[];
Inflation= NaN(size(Susceptibles));
Interest= NaN(size(Susceptibles));
Investment = repelem(s1_in_irf/s1_in_ss,4);

Consumption_ss = 0;
Labour_ss = 0;
Output_ss = 0;
Susceptibles_ss= 0;
Infected_ss= 0;
Recovered_ss= 0;
Deaths_ss= 0;
Inflation_ss= 0;
Interest_ss= 0;
Investment_ss = 0;
save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');
%Paper original results are reported in quarter

%Consumption = month2quarter([zeros(2,1); s1_c_irf./s1_c_ss]); % convert to quartrer from monthly IRFs
%Labour = month2quarter([zeros(2,1);(s1_l1_irf+s1_l2_irf)./(s1_l1_ss+s1_l2_ss)]); % convert to quartrer from monthly IRFs
%Output = month2quarter([zeros(2,1);s1_y_irf/s1_y_ss]); % convert to quartrer from monthly IRFs
%Susceptibles = month2quarter([zeros(2,1); day2month(sum(Smat)')]); % convert to quartrer from monthly IRFs
%Infected = month2quarter([zeros(2,1); day2month(sum(Imat)')]); % convert to quartrer from monthly IRFs
%Recovered = month2quarter([zeros(2,1); day2month(sum(Rmat)')]); % convert to quartrer from monthly IRFs
%Deaths = month2quarter([zeros(2,1); day2month(sum(Dmat)')]); % convert to quartrer from monthly IRFs