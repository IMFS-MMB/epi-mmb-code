% Function to get steady-state wages and proportion of occupations that
% work in market
function [wage_r,pol_mkt_r,wage_im,wage_ih,wage_i,pol_mkt_i,theta_r,theta_i] ...
    = steady_state_cutoff(periods,ini_pop,wage_sm,wage_sh,tau_m,gamma_Im,gamma_Ih)

n_occ = size(wage_sm,2);
wage_r = NaN*ones(periods+1,n_occ); % Optimal wage chosen by recovered
pol_mkt_r = NaN*ones(periods+1,n_occ); % Policy function of recovered of whether to choose to work at home or market
% If policy is 1 then work at home and 2 is work at market

wage_ih = gamma_Ih.*wage_sh; % wage of infected working at market
wage_im = gamma_Im.*wage_sm; % wage of infected working at home

wage_i = NaN*ones(periods+1,n_occ); % Optimal wage chosen by infected
pol_mkt_i = NaN*ones(periods+1,n_occ); % Policy function of infected of whether to choose to work at home or market
% If policy is 1 then work at home and 2 is work at market

for t = 1:1:periods+1
    for i = 1:n_occ
        [~,pol_mkt_r(t,i)] = max([wage_sh(t,i); (1-tau_m(t,i)).*wage_sm(t,i)]);
         wage_r(t,i) = wage_sh(t,i).*(pol_mkt_r(t,i)==1) + wage_sm(t,i).*(pol_mkt_r(t,i)==2); % Pre-tax wage for recovered
        [~,pol_mkt_i(t,i)] = max([wage_ih(t,i); (1-tau_m(t,i))*wage_im(t,i)]);
        wage_i(t,i) = wage_ih(t,i).*(pol_mkt_i(t,i)==1) + wage_im(t,i).*(pol_mkt_i(t,i)==2); % Pre-tax wage for infected
    end    
end

theta_r = 1 - sum(repmat(ini_pop,periods+1,1).*(pol_mkt_r==1),2); % finding the proportion of occupations who go to market for recovered

theta_i = 1 - sum(repmat(ini_pop,periods+1,1).*(pol_mkt_i==1),2); % finding the proportion of occupations who go to market for infected

end