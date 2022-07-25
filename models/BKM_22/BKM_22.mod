@#define capacity = 1
//0 - no limit, 1 - hard limit on health care system capacity, 2 - smooth approximation

@#define lockdown = 0
//0 - no lockdown , 1 - lockdown rule on everybody, 2 - isolation of I, 3 - mix: isolation of I + lockdown on everybody 4: strict lockdown

@#define mp = 0
//0 - react to inf & GAP, 1 - react to inf & GDP

var S I R D T Pop 
    IA RA chi psi 
    pid pir pi_SS pi_SI pi_SR
    U_s c_s n_s lambb_s t lambt 
    U_i c_i n_i lambb_i 
    U_r c_r n_r lambb_r
    w pic 
    p_star Omega Upsilon Delta
    Gamma 
    B_s B_r B_i
    RR
    y c n feas 
    lS lI lR lIA lRA
    tauc_s tauc_i tauc_r taun_s taun_i taun_r
    Util 
    ;
@#include "commonVar.mod"    // added by epi-mmb team

varexo eps_tauc_s eps_tauc_i eps_tauc_r eps_taun_s eps_taun_i eps_taun_r  eps_RR GAP ;%eps_lock

parameters pis1 pis2 pis3 pid_mt pir_mt lim 
            kappa rho piT 
            sigmac thet bet zeta_c phi_Z zeta_n
            A varepsilon delt
            U_d 
            Phi_pic Phi_y Phi_I Phi_I_d Phi_lock_c Phi_lock_n Phi_lock_3
            rho_tauc_s rho_tauc_i rho_tauc_r rho_taun_s rho_taun_i rho_taun_r
	    i_ini;

load inf_ini
i_ini=helper; 
//i_ini=0.001;

//SIR-Macro rescaled such that c=2, n=1
pis3 = 0.5696; 
pis1 = 0.2118; 
pis2 = 1.1853; 

pid_mt = 7*0.006/18; // probability of I dying with medical treatment
pir_mt = 7/18 - pid_mt; // probability of I recovering with medical treatment
lim = 0.01; // max capacity of healhcare system 

U_d = -3863;

rho = 0.6; //rho = 1 no asymptomatic
kappa = 0.5;
piT = 0*0.1;//0

sigmac = 1;
thet = 1.4474; //weight on labor disutility
bet = 0.96^(1/52); //discount factor

zeta_n = 0;     // relative infectiousness of infected via work 
tc=7.65;    // lockdown tax on c
tn=3.8;     // lockdown tax on n
Phi_lock_3 = 3.8; // lockdown impact via 3rd channel 

@#if (lockdown == 0)
    zeta_c = 1;    // relative infectiousness of infected via cons 
    Phi_lock_c = 0; 
    Phi_lock_n = 0; 
@#endif

@#if lockdown == 1
    zeta_c = 1;              
    Phi_lock_c = tc;        
    Phi_lock_n = tn;        
    
@#endif

@#if lockdown == 2
    zeta_c = 0;  
    Phi_lock_c = 0; 
    Phi_lock_n = 0; 
@#endif

@#if lockdown == 3
    zeta_c = 0;  
    Phi_lock_c = tc; 
    Phi_lock_n = tn; 
@#endif

@#if lockdown == 4
    zeta_c = 1;  
    Phi_lock_c = 3*tc; 
    Phi_lock_n = 3*tn; 
@#endif


phi_Z = .8; //productivity of infected

A = 2; //productivity
varepsilon = 6; //EOS between product varieties
delt = 0.75^(1/13); //Calvo probability

Phi_pic = 1.5;
Phi_y = 0.5/52;
Phi_I = 0; 
Phi_I_d = 0;

rho_tauc_s = 0;
rho_tauc_i = 0.99;
rho_tauc_r = 0.99;
rho_taun_s = 0;
rho_taun_i = 0.99;
rho_taun_r = 0.99;

model;

//Death rates
@#if capacity == 1
pid = min(3*pid_mt,pid_mt*(1+I(-1)/(lim))); //increases linearly until 2*capacity constraint, then constant
pir = 7/18-pid;
@#endif

@#if capacity == 2
pid = pid_mt+(2*pid_mt)/(1+exp(-400*(I(-1)-0.01)));  //Logistic
pir = 7/18-pid;
@#endif

@#if capacity == 0
pir = pir_mt;
pid = pid_mt;
@#endif

//Epidemics block

T = S(-1)*t;
S = S(-1) - T;
I = I(-1) + rho*T - (pir+pid)*I(-1);
IA = IA(-1) +(1-rho)*T - (pir+pid)*IA(-1);
R = R(-1) + pir*I(-1) + piT*RA(-1);
RA = (1-piT)*RA(-1)+(pir+pid)*IA(-1);
D = D(-1) + pid*I(-1);
Pop = Pop(-1) - pid*I(-1);

//Households (assuming standard CRRA separable utility)
t =  pis1*c_s*(zeta_c*I(-1)*c_i+kappa*IA(-1)*c_s) + pis2*n_s*(zeta_n*I(-1)*n_i+kappa*IA(-1)*n_s)  + (1-tauc_s)^Phi_lock_3*pis3*(I(-1)+kappa*IA(-1));
psi = ((1-piT)*psi(-1) + (pir+pid)*chi(-1))/(1-t);
chi = (1-pir-pid)*chi(-1)/(1-t) + (1-rho)*t/(1-t);
pi_SS = 1-pi_SI-pi_SR; 
pi_SI = rho*t/(1+chi(-1)+psi(-1)); 
pi_SR = piT*psi(-1)/(1+chi(-1)+psi(-1));

//Susceptible
U_s = log(c_s) + thet*log(1-n_s) + bet*((1-pi_SI-pi_SR)*U_s(1) + pi_SI*U_i(1) + pi_SR*U_r(1));
(S(-1)+IA(-1)+RA(-1))*(c_s + B_s - w*n_s - Gamma) = B_s(-1)*(lS(-1) - rho*(lS(-1)-S(-1)) + lIA(-1) + lRA(-1)*(1-piT)) *RR(-1)/pic;                                                       
c_s^(-1) = (1+tauc_s)*lambb_s - lambt*pis1*(zeta_c*I(-1)*c_i+kappa*IA(-1)*c_s); 
bet*rho*(U_i(1)-U_s(1)) = lambt*(1+chi(-1)+psi(-1));
thet/(1-n_s) = (1-taun_s)*w*lambb_s + lambt*pis2*(zeta_n*I(-1)*n_i+kappa*IA(-1)*n_s);
lambb_s = bet*(1-pi_SI-pi_SR)*lambb_s(1)*RR/pic(1) + bet*pi_SI*lambb_i(1)*RR/pic(1) + bet*pi_SR*lambb_r(1)*RR/pic(1);

//Infected
U_i = log(c_i) + thet*log(1-n_i) + bet*((1-pir-pid)*U_i(+1) + pir*U_r(+1) + pid*U_d);
I(-1)*(c_i + B_i - phi_Z*w*n_i - Gamma) = (B_i(-1)*lI(-1)*(1-pir-0*pid) + B_s(-1)*rho*(lS(-1)-S(-1))) * RR(-1)/pic;    
c_i^(-sigmac) = (1+tauc_i)*lambb_i;
n_i =  max((1-thet/(phi_Z*w*lambb_i)),0);
lambb_i = bet*(1-pir-pid)*lambb_i(1)*RR/pic(1)/(1-pid) + bet*pir*lambb_r(1)*RR/pic(1)/(1-pid);  

//Recovered
U_r = log(c_r) + thet*log(1-n_r) + bet*U_r(+1);
(R(-1))*(c_r + B_r - w*n_r - Gamma) = (B_r(-1)*(lR(-1)) + B_i(-1)*pir*lI(-1) + B_s(-1)*piT*lRA(-1)) * RR(-1)/pic;
c_r^(-sigmac) = (1+tauc_r)*lambb_r;
thet/(1-n_r) = (1-taun_r)*w*lambb_r;
lambb_r = bet*lambb_r(1)*RR/pic(1);

//Firms
p_star = Omega/Upsilon;
Omega/(1-delt*bet) = w/A*y* ((S(-1)+IA(-1)+RA(-1))*c_s^(-sigmac) + I(-1)*c_i^(-sigmac) + R(-1)*c_r^(-sigmac)) + delt*bet*pic(+1)^varepsilon*Omega(+1)/(1-delt*bet); //Here we use weighted marginal utilities
Upsilon/(1-delt*bet) = y*((S(-1)+IA(-1)+RA(-1))*c_s^(-sigmac) + I(-1)*c_i^(-sigmac) + R(-1)*c_r^(-sigmac)) + delt*bet*pic(+1)^(varepsilon-1)*Upsilon(+1)/(1-delt*bet); //Here we use weighted marginal utilities
1 = delt*(pic^(varepsilon-1)) + (1-delt)*p_star^(1-varepsilon);
Delta = delt*pic^varepsilon*Delta(-1) + (1-delt)*p_star^(-varepsilon); //37

//Government
y - w*n = Gamma*(S(-1) + I(-1) + R(-1) + IA(-1) + RA(-1));

//Monetary policy 
RR/(steady_state(RR)) = ((pic/steady_state(pic))^(Phi_pic)*((y/Pop)/(steady_state(y)/steady_state(Pop)+GAP))^(Phi_y)*(exp(I)^Phi_I)*exp(I-I(-1))^Phi_I_d)*exp(eps_RR);

//Market clearing
B_s*(S(-1)+IA(-1)+RA(-1)) + B_r*R(-1) + B_i*I(-1) = 0;
c = (S(-1)+IA(-1)+RA(-1))*c_s + I(-1)*c_i + R(-1)*c_r;
n = (S(-1)+IA(-1)+RA(-1))*n_s + I(-1)*n_i*phi_Z + R(-1)*n_r;

y*Delta = A*n;
feas = y - c;
lS = S(-1);
lI = I(-1);
lR = R(-1);
lIA = IA(-1);
lRA = RA(-1);
@#include "commonVarEq.mod"  // added by epi-mmb team

//Lockdown
@#if (lockdown == 0)
tauc_s = rho_tauc_s*tauc_s(-1) + eps_tauc_s;
tauc_i = rho_tauc_i*tauc_i(-1) + eps_tauc_i;
tauc_r = rho_tauc_r*tauc_r(-1) + eps_tauc_r;
taun_s = rho_taun_s*taun_s(-1) + eps_taun_s;
taun_i = rho_taun_i*taun_i(-1) + eps_taun_i;
taun_r = rho_taun_r*taun_r(-1) + eps_taun_r;
@#endif

@#if (lockdown != 0)
tauc_s = Phi_lock_c*I(-1);
taun_s = Phi_lock_n*I(-1);
tauc_i = tauc_s; 
tauc_r = tauc_s;
taun_i = taun_s; 
taun_r = taun_s;
@#endif

//Aggregate period utility (not used, can be removed)
Util = (S(-1)+IA(-1)+RA(-1))*(log(c_s) + thet*log(1-n_s)) + I(-1)*(log(c_i) + thet*log(1-n_i)) + R(-1)*(log(c_r) + thet*log(1-n_r)) - D(-1)*(log(c_r) - thet*log(1-n_r));

end;

model_diagnostics;

//Simulation

initval;
I = i_ini;
S = 1-i_ini;
R = 0;
D = 0;
Pop = 1;
T = 0;
Delta = 1;
pic = 1;
RR = pic/bet;
@#include "commonVarSS.mod" //added by epi-mmb team
end;
model_diagnostics;

shocks;
var GAP;
periods 1:250;
values 0;
end;


//This section solves the flex price model and saves the GAP
    set_param_value('delt',0);
    set_param_value('Phi_y',0);  //so that stupid guess of gdp_flexprice above does not matter
  perfect_foresight_setup(periods=250);

    load init_covid;

    @#if lockdown == 0 
    oo_.endo_simul=[xx0;oo_.endo_simul(57:65,:)]; // edited by epi-mmb team
    @#else
    oo_.endo_simul=xx_lockdown; 
    @#endif
    perfect_foresight_solver;

    GAP_series=y(2:251)./Pop(2:251)-y(1)/Pop(1);
    yf = y;

@#if mp==0    
    shocks;
    var GAP;
    periods 1:250;
    values (GAP_series);
    end;
@#endif

//This section solves the sticky price model

set_param_value('delt', 0.75^(1/13));
set_param_value('Phi_y',0.5/52);

perfect_foresight_setup(periods=250);

load init_covid;

@#if lockdown == 0 
oo_.endo_simul=[xx0;oo_.endo_simul(57:65,:)]; // edited by epi-mmb team
@#else
oo_.endo_simul=xx_lockdown; //xx0; // xx_lockdown; //xx_test; 
@#endif

perfect_foresight_solver(stack_solve_algo=0,tolf=1e-8, maxit=10,no_homotopy); 
  
OG=(y-yf)./(yf);


// save results

hor = 250;
a=repelem(1,51);
b=repelem(0,51);
PIC=vertcat(a',pic(1:hor));
nomRweek=vertcat(RR(1:hor).^1-RR(1)^1,b');
realRweek=vertcat((1+RR(1:hor)-RR(1)-pic(2:hor+1)+pic(1)).^1-1,b');
inf=(repelem(NaN,hor))';
nR=(repelem(NaN,hor))';
rR=(repelem(NaN,hor))';
for i=1:1:hor
 temp_inf = PIC(i:i+51);
 temp_nr = nomRweek(i:i+51);
 temp_rr = realRweek(i:i+51);
 clear sum;
 inf(i)=prod(temp_inf);
 nR(i)=sum(temp_nr);
 rR(i)=sum(temp_rr);
end;


//save results_no_contaiment_baselineMP             // capacity = 1, lockdown = 0, mp = 0
//save results_lockdown_baselineMP                  // capacity = 1, lockdown = 1, mp = 0
//save results_isolation_baselineMP                 // capacity = 1, lockdown = 2, mp = 0
//save results_isolation_and_lockdown_baselineMP    // capacity = 1, lockdown = 3, mp = 0
//save results_strict_lockdown_baselineMP           // capacity = 1, lockdown = 4, mp = 0

//save results_no_contaiment_standardMP             // capacity = 1, lockdown = 0, mp = 1
//save results_lockdown_standardMP                  // capacity = 1, lockdown = 1, mp = 1
//save results_isolation_standardMP                 // capacity = 1, lockdown = 2, mp = 1
//save results_isolation_and_lockdown_standardMP    // capacity = 1, lockdown = 3, mp = 1
//save results_strict_lockdown_standardMP           // capacity = 1, lockdown = 4, mp = 1

//save results_no_contaiment_optimalMP              // capacity = 1, lockdown = 0, mp = 0, Phi_I = 0.0625
//save results_lockdown_optimalMP                   // capacity = 1, lockdown = 1, mp = 0, Phi_I = 0.0629
//save results_isolation_optimalMP                  // capacity = 1, lockdown = 2, mp = 0, Phi_I = −0.0096
//save results_isolation_and_lockdown_optimalMP     // capacity = 1, lockdown = 3, mp = 0, Phi_I = −0.0272
//save results_strict_lockdown_optimalMP            // capacity = 1, lockdown = 4, mp = 0, Phi_I = −0.0093

 @#include "saveResults.mod"   //added by epi-mmb team
