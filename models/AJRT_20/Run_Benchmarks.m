Model_Params;

Model_Solution;

Consumption = aggC.*pop+aggCstar.*popstar;
Labour = aggL.*pop+aggLstar.*popstar;
Output = z.*aggL.*pop+zstar.*aggLstar.*popstar;
Susceptibles = (S+Sstar)./(pop0+popstar0); 
Infected = (I+Istar)./(pop0+popstar0); 
Recovered = (R+Rstar)./(pop0+popstar0); 
Deaths = pop0-Susceptibles-Infected-Recovered;
Interest=NaN(size(Susceptibles));
Inflation=NaN(size(Susceptibles));
Investment=NaN(size(Susceptibles));

Consumption_ss = (pop0+popstar0)*xrss;
Labour_ss = (pop0+popstar0)*lrss;
Output_ss = (pop0+popstar0)*z*lrss;
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;
Interest_ss=0;
Inflation_ss=0;
Investment_ss=0;

save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');


%save(['../Mat/', job_name, '_Benchmark.mat']);

%plot1;
%suptitle(['Benchmark ', Utext]);
%print(['../Figures/', job_name, '_Benchmark.pdf'],'-dpdf','-fillpage')
%close;
