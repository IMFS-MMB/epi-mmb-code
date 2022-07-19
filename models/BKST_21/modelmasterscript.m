main

Susceptibles = frac_young.*(M_h(i_young,:)+M_fh(i_young,:))+frac_old.*(M_h(i_old,:)+M_fh(i_old,:)); 
Susceptibles=Susceptibles';
Infected = frac_young.*(M_i(i_young,:)+M_fi(i_young,:)+M_s(i_young,:))+frac_old.*(M_i(i_old,:)+M_fi(i_old,:)++M_s(i_old,:));
Infected=Infected';
%Infected = frac_young.*N_c_s(i_young,:)+frac_old.*N_c_s(i_old,:); %infected stock
Recovered = frac_young.*(M_r(i_young,:))+frac_old.*(M_r(i_old,:)); 
Recovered=Recovered';
Deaths = frac_young.*(M_d(i_young,:))+frac_old.*(M_d(i_old,:)); %note: Covid deaths only
Deaths=Deaths';
%Deaths = frac_young.*(M_d(i_young,:)+M_dn(i_young,:))+frac_old.*(M_d(i_old,:)+M_dn(i_old,:)); %note: Covid deaths and natural deaths

Consumption = frac_young.*(c_r(i_young,:).*M_r(i_young,:)+c_h(i_young,:).*M_h(i_young,:)+c_i(i_young,:).*(M_i(i_young,:)+M_s(i_young,:))+c_f(i_young).*(M_fh(i_young,:)+M_fi(i_young,:)))...
    +frac_old.*(c_r(i_old,:).*M_r(i_old,:)+c_h(i_old,:).*M_h(i_old,:)+c_i(i_old,:).*(M_i(i_old,:)+M_s(i_old,:))+c_f(i_old).*(M_fh(i_old,:)+M_fi(i_old,:)));
Consumption=Consumption';
Labour = frac_young.*((n_r(i_young,:)+v_r(i_young,:)).*M_r(i_young,:)+(n_h(i_young,:)+v_h(i_young,:)).*M_h(i_young,:)+(n_i(i_young,:)+v_i(i_young,:)).*(M_i(i_young,:)+M_s(i_young,:))+(n_f(i_young,:)+v_f(i_young,:)).*(M_fh(i_young,:)+M_fi(i_young,:)));
Labour=Labour';
Output = gdp';
Interest = NaN(length(Consumption),1);
Inflation = NaN(length(Consumption),1);
Investment = NaN(length(Consumption),1);


Consumption_ss = frac_young.*c_r(i_young,1)+frac_old.*c_r(i_old,1);
Labour_ss = frac_young.*(n_r(i_young,1)+v_r(i_young,1));
Output_ss = w*(n_r(i_young,1)+tau_r(i_young,1));
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;
Interest_ss = 0;
Inflation_ss =0;
Investment_ss=0;


save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered','Interest','Inflation','Investment');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss','Interest_ss','Inflation_ss','Investment_ss');

