/* 
New Keynesian model with SIR-type epidemic dynamics & stochastic vaccine arrival
The Limited Power of Monetary Policy in a Pandemic, Antoine Lepetit and Cristina Fuentes-Albero, May 2022
Version with Taylor Rule + effective lower bound 
*/

var t s i r d R0 tau w infl ffr ffr_shadow y lambda c c_s c_i c_r n_s n_i n_r V_s V_i V_r mu;
@#include "commonVar.mod"    // added by epi-mmb team

varexo mp_shock covid_19; 

parameters pi_s1 pi_s2 pi_s3 pi_r pi_d R0_ss share
           beta teta phip delta_infl delta_gap xi varphi sigma phi A infl_target flow_value_life delta_v subsidy
           ffr_init n_s_init n_i_init n_r_init w_init c_s_init c_i_init c_r_init V_s_init V_i_init V_r_init y_init lambda_init
           ffr_end n_s_end n_i_end n_r_end w_end c_s_end c_i_end c_r_end V_s_end V_i_end V_r_end y_end lambda_end
           s_end r_end d_end i_ini;


load inf_ini
i_ini=helper; 
//i_ini=0.005;


// SIR PARAMETERS // 

pi_d = 7*0.005/15;
pi_r = 7/15 - pi_d;
A = 1;
n_s_init = 1;
share = 1/4;
R0_ss = 1.45;

pi_s3 = R0_ss*(pi_r+pi_d)/(1 + 2*(share/(1-share))*(1+share/(1-share))/(1 - (share/(1-share))^2));

delta_v = 1/52;

// NK PARAMETERS //

beta = 0.9996;
teta = 6;
phip = 40926;
delta_infl = 1.5;
delta_gap = 0.5/52;
varphi = 0.5;
sigma = 2;
phi = 0.8;
infl_target = 1.02^(1/52);
subsidy = 1/(teta-1);

// INITIAL STEADY STATE //

ffr_init = infl_target/beta - 1;
w_init = 1;
xi = w_init/((n_s_init^(1/varphi+sigma))*(A^sigma));
c_s_init = A*n_s_init;
c_i_init = c_s_init;
c_r_init = c_s_init;
n_i_init = (phi*w_init/(xi*(c_i_init^sigma)))^varphi;
n_r_init = n_s_init;
y_init = A*n_s_init;
lambda_init = c_s_init^(-sigma);
flow_value_life = (1-beta)*10310 - (c_s_init^(1-sigma))/(1-sigma) + xi*(n_s_init^(1+1/varphi))/(1+1/varphi);
V_r_init = ((c_r_init^(1-sigma))/(1-sigma) - xi*(n_r_init^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_init*(w_init*n_s_init - c_s_init))/(1-beta);
V_i_init = ((c_i_init^(1-sigma))/(1-sigma) - xi*(n_i_init^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_init*(w_init*phi*n_i_init - c_i_init) + beta*pi_r*V_r_init)/(1-beta*(1 - pi_r - pi_d));
V_s_init = ((c_s_init^(1-sigma))/(1-sigma) - xi*(n_s_init^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_init*(w_init*n_r_init - c_r_init) + beta*delta_v*V_r_init)/(1-beta*(1-delta_v));

pi_s1 = (share/(1-share))*(1+share/(1-share))*pi_s3/((1 - (share/(1-share))^2)*(c_s_init*c_i_init));
pi_s2 = (share/(1-share))*(1+share/(1-share))*pi_s3/((1 - (share/(1-share))^2)*(n_s_init*n_i_init));

// FINAL STEADY STATE //

// Enter a guess for number of deaths (as a fraction of the population), number of susceptibles & number of recovered
s_end = 0.5868;
r_end = 1-s_end;
d_end = 0.003289;

w_end = w_init;
n_s_end = n_s_init;
n_i_end = n_i_init;
n_r_end = n_r_init;
y_end = A*(s_end*n_s_end + r_end*n_r_end);
c_s_end = c_s_init;
c_i_end = c_i_init;
c_r_end = c_r_init;
ffr_end = ffr_init;
lambda_end = c_s_end^(-sigma);
V_r_end = ((c_r_end^(1-sigma))/(1-sigma) - xi*(n_r_end^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_end*(w_end*n_s_end - c_s_end))/(1-beta);
V_i_end = ((c_i_end^(1-sigma))/(1-sigma) - xi*(n_i_end^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_end*(w_end*phi*n_i_end - c_i_end) + beta*pi_r*V_r_end)/(1-beta*(1 - pi_r - pi_d));
V_s_end = ((c_s_end^(1-sigma))/(1-sigma) - xi*(n_s_end^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda_end*(w_end*n_r_end - c_r_end) +  beta*delta_v*V_r_end)/(1-beta*(1-delta_v));

model;

// SIR BLOCK //

[name = 'Number of newly-infected']
t = pi_s1*s(-1)*c_s*i(-1)*c_i + pi_s2*s(-1)*n_s*i(-1)*n_i + pi_s3*s(-1)*i(-1);

[name = 'Law of motion of susceptibles']
s = s(-1) - t - covid_19;

[name = 'Law of motion of infected']
i = (1 - pi_r - pi_d)*i(-1) + t + covid_19;
   
[name = 'Law of motion of recovered']
r = r(-1) + pi_r*i(-1);

[name = 'Law of motion of deaths']
d = d(-1) + pi_d*i(-1);

[name = 'Probability of getting infected']
tau = t/s(-1);

[name = 'R0']
R0 = (pi_s1*c_s*c_i + pi_s2*n_s*n_i + pi_s3)/(pi_r + pi_d);

// NK BLOCK //

// Aggregate equations

[name = 'Price Phillips curve']
(1 - teta)*(1+subsidy) + teta*w/A - phip*infl*(infl-infl_target) + beta*phip*infl(+1)*(infl(+1)-infl_target)*y(+1)/y = 0;

[name = 'Resource constraint']
y*(1 - (phip/2)*((infl-infl_target)^2)) = s(-1)*c_s + i(-1)*c_i + r(-1)*c_r;

[name = 'Production function']
y = A*(s(-1)*n_s + i(-1)*phi*n_i + r(-1)*n_r);

[name = 'Aggregate consumption']
c = s(-1)*c_s + i(-1)*c_i + r(-1)*c_r;

[name = 'Shadow rate', mcp = 'ffr_shadow > 0']
1+ffr_shadow = (1+ffr_init)*((infl/infl_target)^delta_infl)*((y/y_end)^delta_gap) + mp_shock;

[name = 'Monetary policy']
1+ffr = 1+ffr_shadow;

// Households

[name = 'Marginal utility of consumption - susceptibles']                      
lambda*(1+mu) = c_s^(-sigma) +  pi_s1*i(-1)*c_i*beta*(V_i(+1) - (1-delta_v)*V_s(+1) - delta_v*V_r(+1));

[name = 'Marginal utility of consumption - infected']                      
lambda*(1+mu) = c_i^(-sigma);

[name = 'Marginal utility of consumption - recovered']                      
lambda*(1+mu) = c_r^(-sigma);

[name = 'Labor supply - susceptibles']                      
lambda*w = xi*(n_s^(1/varphi)) - pi_s2*i(-1)*n_i*beta*(V_i(+1) - (1-delta_v)*V_s(+1) - delta_v*V_r(+1));

[name = 'Labor supply - infected']                      
lambda*w*phi = xi*(n_i^(1/varphi));

[name = 'Labor supply - recovered']                      
lambda*w = xi*(n_r^(1/varphi));

[name = 'Euler equation']
lambda = beta*lambda(+1)*(1+ffr)/infl(+1);

[name = 'Lifetime utility of a susceptible']
V_s = (c_s^(1-sigma))/(1-sigma) - xi*(n_s^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda*(w*n_s - (1+mu)*c_s) + (1-tau)*(1-delta_v)*beta*V_s(+1) + tau*beta*V_i(+1) + (1-tau)*delta_v*beta*V_r(+1);

[name = 'Lifetime utility of an infected']
V_i = (c_i^(1-sigma))/(1-sigma) - xi*(n_i^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda*(w*phi*n_i - (1+mu)*c_i) + (1-pi_r-pi_d)*beta*V_i(+1) + pi_r*beta*V_r(+1);

[name = 'Lifetime utility of a recovered']
V_r = (c_r^(1-sigma))/(1-sigma) - xi*(n_r^(1+1/varphi))/(1+1/varphi) + flow_value_life + lambda*(w*n_r - (1+mu)*c_r) + beta*V_r(+1);

[name = 'Containment rate process']
mu = 0*i;

@# include "commonVarEq.mod"  // added by epi-mmb team
end;

initval;
s = 1;
i = 0;
r = 0;
end;

resid;

endval;
t = 0;
s = s_end;
i = 0;
r = r_end;
d = d_end;
tau = 0;
R0 = (pi_s1*c_s_end*c_i_end + pi_s2*n_s_end*n_i_end + pi_s3)/(pi_r + pi_d);
w = w_end;
infl = infl_target;
ffr = ffr_end;
ffr_shadow = ffr_end;
y = y_end;
lambda = lambda_end;
c = y_end;
c_s = c_s_end;
c_i = c_i_end; 
c_r = c_r_end;
n_s = n_s_end;
n_i = n_i_end;
n_r = n_r_end;
V_s = V_s_end; 
V_i = V_i_end; 
V_r = V_r_end;
mu = 0;
@# include "commonVarSS.mod" //added by epi-mmb team
end;

resid;

check;

shocks;
var covid_19;
periods 1;
values (i_ini);
var mp_shock;
periods 1:500;
values 0;
end;

perfect_foresight_setup(periods=500);
perfect_foresight_solver(lmmcp,noprint);
//simul(periods=500,lmmcp);
@# include "saveResults.mod"   //added by epi-mmb team