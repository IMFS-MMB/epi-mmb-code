
function [ys,params,check] = BKM_22_steadystate(ys,exe,M_,options)

%This function calculates the NO-INFECTION steady state


npar = size(M_.param_names,1);                            
for jj = 1:npar                             
  paramname = M_.param_names{jj};                
  eval([ paramname ' = M_.params(' int2str(jj) ');']);        
end                                                       
check = 0;
thet_new = thet; %(thet)^(1/2)-1;
pir = pir_mt;
pid = pid_mt;

Pop = 1;
I = 0;
S = 1;
R = 0; 
D = 0; 
T = 0;
IA= 0;
RA= 0;

chi=0;
psi=0;
pi_SS = 1; 
pi_SI = 0; 
pi_SR = 0;

lI = 0;
lS = 1;
lR = 0; 
lIA= 0;
lRA= 0;

t = 0;
muc = 0;
Gamma = 0;

%w = (varepsilon-1)/varepsilon*A;
w = A;
p_star = 1;
pic = 1;
RR = pic/bet;

n_s = 1/(1+thet_new); %(1/thet)^(1/2);
c_s = w*n_s;
lambb_s = c_s^(-sigmac);
B_s = 0;

n_r = 1/(1+thet_new);
c_r = w*n_r;
lambb_r = c_r^(-sigmac);
B_r = 0;
                      

f = @(c_i) 1 - exp(-c_i)/(1/RR-1)*((1-pir-pid) + pir*lambb_r*c_i^(sigmac));

c_i = fzero(f,2);
n_i = 0; %1-thet_new*c_i/w;
n_i =  max((1-thet_new*c_i/w/phi_Z),0);
fprintf('n_i   = %12.6f ; \n',1-thet_new*c_i/w/phi_Z);
lambb_i = c_i^(-sigmac);
B_i = (phi_Z*w*n_i - c_i)/(1/RR-1);

c = S*c_s + I*c_i + R*c_r;
n = S*n_s + I*phi_Z*n_i + R*n_r;
y = c;
Consumption = c;
Labour = n;
Output = y;
Inflation = pic; 
Interest = RR;
Susceptibles = 1;
Infected = 0;
Recovered = 0;
Deaths = 0;


%U_s = (c_s^(1-sigmac)/(1-sigmac) - thet*n_s^(1+sigman)/(1+sigman))/(1-bet);
%U_r = (c_r^(1-sigmac)/(1-sigmac) - thet*n_r^(1+sigman)/(1+sigman))/(1-bet);
%U_i = (c_i^(1-sigmac)/(1-sigmac) - thet*n_i^(1+sigman)/(1+sigman) + bet*pir*U_r)/(1-bet*(1-pir-pid));


U_s = (log(c_s) + thet*log(1-n_s))/(1-bet);
pi_dd = 1/1000000/52; %1 micromort translated into weeks
x_dd = 10/50413; %percent of income agent willing to pay for certain life instead of 1 micromort probability of death
U_dd = U_s-log(1+x_dd)/pi_dd/bet;
fprintf('U_dd = %15.6f; \n', U_dd);
U_r = (log(c_r) + thet*log(1-n_r))/(1-bet);
U_i = (log(c_i) + thet*log(1-n_i)+ bet*pir*U_r + bet*pid*U_d)/(1-bet*(1-pir-pid)); 

lambt = bet*rho*(U_i-U_s);

Delta = 1;
Omega = w/A*y*c^(-sigmac);
Upsilon = y*c^(-sigmac);

feas = 0;
lock = 0;
Util = log(c_s) + thet*log(1-n_s);
n_i_con = n_i + lock;
n_s_con = n_s + lock;

tauc_s = 0;
tauc_i = 0;
tauc_r = 0;
taun_s = 0;
taun_i = 0;
taun_r = 0;

params=NaN(npar,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

nvar = size(M_.endo_names,1);                   
ys = zeros(nvar,1);                   
for jj = 1:nvar                   
  varname = deblank(M_.endo_names{jj});   
  eval(['ys(' int2str(jj) ') = ' varname ';']);   
end                                                         

%  for jj=1:nvar 
%      disp([M_.endo_names(jj,:),'= ',num2str(ys(jj,1),100),';']); 
%  end
