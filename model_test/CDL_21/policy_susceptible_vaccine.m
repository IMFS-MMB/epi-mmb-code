% Function that does policy function iteration to get equilibrium for
% model with vaccine and treatment
function [pol_mkt_s_opt,wage_s_opt,RnotSIR,I,Ih,Im,S,Sh,Sm,R,D,V_r,V_i,V_sh,V_sm,Cs,Ci,Cr,C,G,Y,theta_total,error] ...
    = policy_susceptible_vaccine(pol_mkt_s_guess,wage_s_guess,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods)

stop=0;
iter = 0;% iterations
max_iter = 1000;%maximum iterations
tol = 0.00001;% tolerance level for error

n_occ = size(occ_weight,2); % number of occupations

while stop==0
    
    iter = iter + 1;
    
    [pol_mkt_s_opt,wage_s_opt,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = sir_dynamics_vaccine(pol_mkt_s_guess,wage_s_guess,...
        sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,...
        risk_prop,ini_pop,betta,pi_par,periods);
          
   error = abs(pol_mkt_s_opt - pol_mkt_s_guess);
   
   if ( max(max(error))<tol || iter>max_iter)    
       stop = 1;
   else
       pol_mkt_s_guess = pol_mkt_s_opt;
       wage_s_guess = wage_s_opt;
   end
   
   
end
% Generating final results when policy functions have converged
[pol_mkt_s_opt,wage_s_opt,RnotSIR,I,Ih,Im,S,Sh,Sm,R,D,V_r,V_i,V_sh,V_sm,Cs,Ci,Cr,C,G,Y,theta_total] ...
    = sir_dynamics_vaccine(pol_mkt_s_opt,wage_s_opt,sec_weight,sec_elas,occ_weight,occ_elas,tau_m,wage_sm,wage_sh,wage_r,...
    wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,risk_prop,ini_pop,betta,pi_par,periods);

end
 
