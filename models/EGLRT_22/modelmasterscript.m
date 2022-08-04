young_old 

Consumption_month= par.sy_ss*cons_monthly.Cy.Consy+(1- par.sy_ss)*cons_monthly.Co.Conso; 
Consumption=repelem(Consumption_month,4)
Consumption(59:end)=[];
Consumption_ss=0; 

Infected= par.sy_ss*Iy+(1-par.sy_ss)*Io; 
Infected_ss=0; 

Susceptibles=par.sy_ss*Sy+(1-par.sy_ss)*So; 
Susceptibles_ss=0; 

Recovered=par.sy_ss*Ry+(1-par.sy_ss)*Ro;
Recovered_ss=0;

Deaths=par.sy_ss*Dy+(1-par.sy_ss)*Do; 
Deaths_ss=0; 

save('simulated_results.mat','Consumption','Deaths','Susceptibles','Infected','Recovered');
save('simulated_results_ss.mat','Consumption_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');