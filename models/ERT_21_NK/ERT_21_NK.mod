//This file solves and simulates the model developed in Eichenbaum,
//Rebelo and Trabandt (2020),'Epidemics in the New Keynesian Model'.

//This code was written and run with Matlab R2021a and DYNARE 4.6.1 The nonlinear
//model is solved and simulated with Dynare's nonlinear deterministic solver.

////////////////////////////////////////////////////////////
//IMPORTANT: You must use DYNARE 4.6.1 to run this code.//// 
////////////////////////////////////////////////////////////

//We solve the fully nonlinear model. Given the strong nonlinearities of the 
//model due to the epidemic dynamics we use homotopy for a number of 
//parameters, i.e. we increase parameter values stepwise until we impose
//their final desired value. Homotopy calculations can take a while. 

//This code only produces plots for the specific model below. The plots 
//provided in the paper are produced using the figs05_...through figs10_... 
//files in the main directory of this code kit.

//Mathias Trabandt, mathias.trabandt@gmail.com


//Endogenous variables of the actual (sticky price) economy
var y k n w rk x c s i r ns ni nr cs ci cr tau 
lambtilde lamtau lami lams lamr dd pop Rb pie mc F Kf rr 
dcs dns dci dni dw dlams dlamtau dlambtilde dlami dlamr 
dcr dnr drk dF dKf pbreve dpie;

//Endogenous variables of the flexible price economy
var yF kF nF wF rkF xF cF sF iF rF nsF niF nrF csF ciF crF 
tauF lambtildeF lamtauF lamiF lamsF lamrF ddF popF RbF pieF 
mcF FF KfF rrF dcsF dnsF dciF dniF dwF dlamsF dlamtauF 
dlambtildeF dlamiF dlamrF dcrF dnrF drkF dFF dKfF pbreveF dpieF;  
@#include "commonVar.mod"    // added by epi-mmb team 
var Interest Inflation;


parameters 
xi rpi rx gam pi1 pi2 pi3 pir pid betta i_ini A theta 
alfa inc_target n_target delta g_ss eta xi_flex;


//Initialize parameters
betta=0.98^(1/52); //Weekly household discount factor
pid=7*0.002/14;    //Weekly probability of dying
pir=7*1/14-pid;    //Weekly probability of recovering

delta=0.06/52;     //capital depreciation rate (weekly)
alfa=2/3;          //labor share

gam=1.35;          //Steady state gross price markup; Keep in mind that
                   //gam must be larger than unity. 

xi=0.98;           //Calvo price stickiness in actual economy (weekly)
xi_flex=0;         //Calvo price stickiness in flexible price economy (weekly)

rpi=1.5;           //Taylor rule coefficient inflation
rx=0.5/52;         //Taylor rule coefficient output gap     

n_target=28;        //Calibration target for weekly hours
inc_target=58000/52;//Calibration target for weekly income
eta=0.19;           //Calibration target for gov. cons to GDP ratio 

load inf_ini
i_ini=helper; 
//i_ini=0.001;        //Initial seed of infection


format long;
//Calibation targets for shares of pi-terms in T-function in SIR model
pi3_shr_target=2/3;                   //share of T_0 jump due general infections
pi1_shr_target=(1-pi3_shr_target)/2;  //share of T_0 jump due to consumption-based infections
pi2_shr_target=(1-pi3_shr_target)/2;  //share of T_0 jump due to work-based infections
RplusD_target=0.60;                   //total share of people infected and then either recovered or dead after epidemic


//pre-infection steady states
y_ss=inc_target;
n_ss=n_target;
mc_ss=1/gam;
w_ss=mc_ss*alfa*y_ss/n_ss;
rk_ss=1/betta-1+delta;
kn_ss=(1-alfa)*w_ss/alfa/rk_ss;
yk_ss=(y_ss/n_ss)/kn_ss;
A=(y_ss/n_ss)^alfa*yk_ss^(1-alfa);
k_ss=(y_ss/A/n_ss^alfa)^(1/(1-alfa));
x_ss=delta*k_ss;
g_ss=eta*y_ss;
c_ss=(1-eta)*y_ss-x_ss;
g_ss=y_ss-c_ss-x_ss;
ns_ss=n_ss;
cs_ss=c_ss;
tau_ss=0;
s_ss=1;
i_ss=0;
r_ss=0;
lambtilde_ss=1/cs_ss;
ci_ss=cs_ss;
cr_ss=cs_ss;
theta=lambtilde_ss*w_ss/ns_ss;
ni_ss=lambtilde_ss*w_ss/theta;
nr_ss=ns_ss;
lams_ss=(log(cs_ss)-theta/2*ns_ss^2 + lambtilde_ss*(w_ss*ns_ss-cs_ss) ) / ( 1/betta-1 );
lamr_ss=(log(cr_ss)-theta/2*nr_ss^2 + lambtilde_ss*(w_ss*nr_ss-cr_ss) ) / ( 1/betta-1 );
lami_ss=(log(ci_ss)-theta/2*ni_ss^2 + lambtilde_ss*(w_ss*ni_ss-ci_ss) + pir*lamr_ss) / ( 1/betta-1+pir+pid );
lamtau_ss=lami_ss-lams_ss;
rr_ss=1/betta;
dcs_ss=1;
dns_ss=1;
dci_ss=1;
dni_ss=1;
dw_ss=1;
dlams_ss=1;
dlamtau_ss=1; 
dlambtilde_ss=1;
dlami_ss=1;
dlamr_ss=1;
dcr_ss=1;
dnr_ss=1;
drk_ss=1;
Kf_ss=1/(1-betta*xi)*gam*mc_ss*lambtilde_ss*y_ss;
F_ss=1/(1-betta*xi)*lambtilde_ss*y_ss;
pie_ss=1;
Rb_ss=rr_ss;
KfF_ss=1/(1-betta*xi_flex)*gam*mc_ss*lambtilde_ss*y_ss;
FF_ss=1/(1-betta*xi_flex)*lambtilde_ss*y_ss;

Consumption_ss = c_ss;            //added by epi-mmb team
Labour_ss = n_ss;
Output_ss = y_ss;
Susceptibles_ss = s_ss; 
Infected_ss = i_ss;
Recovered_ss = r_ss;
//Deaths_ss = d_ss; //not available

//Some useful command window output
cons_share=c_ss/y_ss
inv_share=x_ss/y_ss
gov_share=g_ss/y_ss
value_of_life=1/(1-betta)*(log(c_ss)-theta/2*n_ss^2)*c_ss
ann_capoutputratio=k_ss/(52*y_ss)
thetaval=theta
Aval=A
 

//calibrate the pi's in the transmission (tau) - function
go_calibrate_pi;
 

//final numbers for pi1,pi2,pi3 and xi will be imposed below using homotopy
pi1_final=pi1;
pi2_final=pi2;
pi3_final=pi3;
xi_final=xi;


//put scaled down values of pi1,pi2,pi3 into Dynare M_. structure
//if you put the final numbers from the get go, no solution will be found
//use smaller numbers first then compute a solution. Then use the solution 
//as an initial guess for slightly larger values. Proceed using this homotopy 
//until final values are imposed (all below, after model block). 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final/3;
//Homotopy setup for xi is done after the model block.


model;
//Keep in mind that Dynare wants state variables denoted with t-1 which differs 
//from the timing convention used in the manuscript where state variables 
//are denoted with time subscript t. All model equations have been implemented
//using the Dynare timing convention. 
 
//Note that we use Dynare's time stacking (two point boundary value) algorithm to 
//solve the nonlinear model. The two points are the initial pre-infection steady state
//and the terminal steady state. Dynare solves for the model transition dynamics between 
//these two points that satisfy all nonlinear model equations. By default, Dynare 
//assumes that the terminal steady state is the same as the initial steady state which 
//is not true in our model. The economy will not return to the initial (pre-infection) 
//steady state after the epidemic has ended. Since the terminal steady state depends on 
//epidemic dynamics and policies, we dont know the terminal value of the new steady state.
//In other words the terminal steady state is an endogenous function of the transition 
//dynamics. Now, the model has many level variables dated t+1 which implies that Dynare will
//replace the value of these variables at the end of the simulation by their pre-infection steady
//state variables. We dont want this for the reason described above. So, we dont write 
//level variables dated t+1 but use a transformation of that variable. Specifically, 
//we replace model variables, say, X(+1) by X*dX(+1) and introduce the new variable dX=X/X(-1). 
//Note that when substituting dX(+1) back into X*dX(+1) we get the original model formulation, 
//i.e. X(+1). With this transformation, Dynare now replaces dX(+1) at the end of the simulation 
//with its pre-infection steady state which is correct since the growth rate of a variable 
//in the pre- and terminal steady states is the same.


/////////////////////////////////////////////////////// 
//equilibrium equations: actual (sticky price) economy  
//////////////////////////////////////////////////////
//Production 
y=pbreve*A*k(-1)^(1-alfa)*n^alfa;

//Marginal cost
mc=1/(A*alfa^alfa*(1-alfa)^(1-alfa))*w^alfa*rk^(1-alfa);

//Cost mininizing inputs
w=mc*alfa*A*n^(alfa-1)*k(-1)^(1-alfa);

//Law of motion capital
k=x+(1-delta)*k(-1);

//Aggregate resources
y=c+x+g_ss;

//Aggregate hours
n=s(-1)*ns+i(-1)*ni+r(-1)*nr;

//Aggregate consumption
c=s(-1)*cs+i(-1)*ci+r(-1)*cr;

//New infections
tau=pi1*s(-1)*cs*i(-1)*ci+pi2*s(-1)*ns*i(-1)*ni+pi3*s(-1)*i(-1);

//Total suseptibles
s=s(-1)-tau;

//Total infected
i=i(-1)+tau-(pir+pid)*i(-1);

//Total recovered
r=r(-1)+pir*i(-1);

//Total deaths 
dd=dd(-1)+pid*i(-1);

//Total population
pop=pop(-1)-pid*i(-1);

//First order condition, consumption susceptibles
1/cs=lambtilde-lamtau*pi1*i(-1)*ci;

//First order condition, consumption infected
1/ci=lambtilde;

//First order condition, consumption recovered
1/cr=lambtilde;

//First order condition, hours susceptibles
theta*ns=lambtilde*w+lamtau*pi2*i(-1)*ni;

//First order condition, hours infected
theta*ni=lambtilde*w;

//First order condition, hours recovered
theta*nr=lambtilde*w;

//First order condition, capital
lambtilde=betta*(rk*drk(+1)+(1-delta))*lambtilde*dlambtilde(+1);

//First order condition, new infected
lami=lamtau+lams;

//First order condition, susceptibles
log(cs*dcs(+1))-theta/2*(ns*dns(+1))^2
+lamtau*dlamtau(+1)*(pi1*cs*dcs(+1)*i*ci*dci(+1)+pi2*ns*dns(+1)*i*ni*dni(+1)+pi3*i)
+lambtilde*dlambtilde(+1)*( w*dw(+1)*ns*dns(+1)-cs*dcs(+1) )
-lams/betta+lams*dlams(+1);

//First order condition, infected
log(ci*dci(+1))-theta/2*(ni*dni(+1))^2
+lambtilde*dlambtilde(+1)*( w*dw(+1)*ni*dni(+1)-ci*dci(+1) )
-lami/betta+lami*dlami(+1)*(1-pir-pid)+lamr*dlamr(+1)*pir;

//First order condition, recovered
log(cr*dcr(+1))-theta/2*(nr*dnr(+1))^2
+lambtilde*dlambtilde(+1)*( w*dw(+1)*nr*dnr(+1)-cr*dcr(+1) )
-lamr/betta+lamr*dlamr(+1);

//First order condition, bonds
lambtilde=betta*Rb/(pie*dpie(+1))*lambtilde*dlambtilde(+1);

//Real interest rate
rr=Rb/(pie*dpie(+1));

//Nonlinear price setting 1
Kf=gam*mc*lambtilde*y+betta*xi*(pie*dpie(+1))^(gam/(gam-1))*Kf*dKf(+1);

//Nonlinear price setting 2
F=lambtilde*y+betta*xi*(pie*dpie(+1))^(1/(gam-1))*F*dF(+1);

//Nonlinear price setting 3
Kf=F*( (1-xi*pie^(1/(gam-1)) ) / (1-xi) )^(-(gam-1));

//Inverse price dispersion
pbreve^(-1)=(1-xi)*( (1-xi*pie^(1/(gam-1)))/(1-xi) )^gam
 + xi*pie^(gam/(gam-1))/pbreve(-1);

//Taylor rule
Rb=STEADY_STATE(Rb)+rpi*log(pie/STEADY_STATE(pie))+rx*log(y/yF);

//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dF=F/F(-1);
dKf=Kf/Kf(-1);  
dpie=pie/pie(-1);
dcs=cs/cs(-1);
dns=ns/ns(-1);
dci=ci/ci(-1);
dni=ni/ni(-1);
dw=w/w(-1);
dlams=lams/lams(-1);
dlamtau=lamtau/lamtau(-1);
dlambtilde=lambtilde/lambtilde(-1);
dlami=lami/lami(-1);
dlamr=lamr/lamr(-1);
dcr=cr/cr(-1);
dnr=nr/nr(-1);
drk=rk/rk(-1);


//////////////////////////////////////////////// 
//equilibrium equations: flexible price economy
//////////////////////////////////////////////// 
//Production
yF=pbreveF*A*kF(-1)^(1-alfa)*nF^alfa;

//Marginal cost
mcF=1/(A*alfa^alfa*(1-alfa)^(1-alfa))*wF^alfa*rkF^(1-alfa);

//Cost minimizing inputs
wF=mcF*alfa*A*nF^(alfa-1)*kF(-1)^(1-alfa);

//Law of motion for capital
kF=xF+(1-delta)*kF(-1);

//Aggregate Resources
yF=cF+xF+g_ss;

//Aggregate hours
nF=sF(-1)*nsF+iF(-1)*niF+rF(-1)*nrF;

//Aggregate consumption
cF=sF(-1)*csF+iF(-1)*ciF+rF(-1)*crF;

//New infected
tauF=pi1*sF(-1)*csF*iF(-1)*ciF+pi2*sF(-1)*nsF*iF(-1)*niF+pi3*sF(-1)*iF(-1);

//Total suseptibles
sF=sF(-1)-tauF;

//Total infected
iF=iF(-1)+tauF-(pir+pid)*iF(-1);

//Total recovered
rF=rF(-1)+pir*iF(-1);

//Total deaths 
ddF=ddF(-1)+pid*iF(-1);

//Total population
popF=popF(-1)-pid*iF(-1);

//First order condition, consumption susceptibles
1/csF=lambtildeF-lamtauF*pi1*iF(-1)*ciF;

//First order condition, consumption infected
1/ciF=lambtildeF;

//First order condition, consumption recovered
1/crF=lambtildeF;

//First order condition, hours susceptibles
theta*nsF=lambtildeF*wF+lamtauF*pi2*iF(-1)*niF;

//First order condition, hours infected
theta*niF=lambtildeF*wF;

//First order condition, hours recovered
theta*nrF=lambtildeF*wF;

//First order condition, capital
lambtildeF=betta*(rkF*drkF(+1)+(1-delta))*lambtildeF*dlambtildeF(+1);

//First order condition, new infected
lamiF=lamtauF+lamsF;

//First order condition, susceptibles
log(csF*dcsF(+1))-theta/2*(nsF*dnsF(+1))^2
+lamtauF*dlamtauF(+1)*(pi1*csF*dcsF(+1)*iF*ciF*dciF(+1)+pi2*nsF*dnsF(+1)*iF*niF*dniF(+1)+pi3*iF)
+lambtildeF*dlambtildeF(+1)*( wF*dwF(+1)*nsF*dnsF(+1)-csF*dcsF(+1) )
-lamsF/betta+lamsF*dlamsF(+1);

//First order condition, infected
log(ciF*dciF(+1))-theta/2*(niF*dniF(+1))^2
+lambtildeF*dlambtildeF(+1)*( wF*dwF(+1)*niF*dniF(+1)-ciF*dciF(+1) )
-lamiF/betta+lamiF*dlamiF(+1)*(1-pir-pid)+lamrF*dlamrF(+1)*pir;

//First order condition, recovered
log(crF*dcrF(+1))-theta/2*(nrF*dnrF(+1))^2
+lambtildeF*dlambtildeF(+1)*( wF*dwF(+1)*nrF*dnrF(+1)-crF*dcrF(+1) )
-lamrF/betta+lamrF*dlamrF(+1);

//First order condition, bonds
lambtildeF=betta*RbF/(pieF*dpieF(+1))*lambtildeF*dlambtildeF(+1);

//Real interest rate
rrF=RbF/(pieF*dpieF(+1));

//Nonlinear price setting 1
KfF=gam*mcF*lambtildeF*yF+betta*xi_flex*(pieF*dpieF(+1))^(gam/(gam-1))*KfF*dKfF(+1);

//Nonlinear price setting 2
FF=lambtildeF*yF+betta*xi_flex*(pieF*dpieF(+1))^(1/(gam-1))*FF*dFF(+1);

//Nonlinear price setting 3
KfF=FF*( (1-xi_flex*pieF^(1/(gam-1)) ) / (1-xi_flex) )^(-(gam-1));

//Inverse price dispersion
pbreveF^(-1)=(1-xi_flex)*( (1-xi_flex*pieF^(1/(gam-1)))/(1-xi_flex) )^gam
 + xi_flex*pieF^(gam/(gam-1))/pbreveF(-1);

//Taylor rule
RbF=STEADY_STATE(Rb)+rpi*log(pieF/STEADY_STATE(pieF))+rx*log(yF/yF);

//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dFF=FF/FF(-1);
dKfF=KfF/KfF(-1);  
dpieF=pieF/pieF(-1);
dcsF=csF/csF(-1);
dnsF=nsF/nsF(-1);
dciF=ciF/ciF(-1);
dniF=niF/niF(-1);
dwF=wF/wF(-1);
dlamsF=lamsF/lamsF(-1);
dlamtauF=lamtauF/lamtauF(-1);
dlambtildeF=lambtildeF/lambtildeF(-1);
dlamiF=lamiF/lamiF(-1);
dlamrF=lamrF/lamrF(-1);
dcrF=crF/crF(-1);
dnrF=nrF/nrF(-1);
drkF=rkF/rkF(-1);
@# include "commonVarEq.mod"  // added by epi-mmb team
Interest=Rb;
Inflation=pie;
end; 



//initial (pre-infection) steady state used in simulations 
initval;
y=y_ss; 
k=k_ss;
n=n_ss;
w=w_ss;
rk=rk_ss;
x=x_ss;
c=c_ss;
s=1;
i=0; 
r=0; 
ns =ns_ss;
ni=ni_ss; 
nr=nr_ss; 
cs=cs_ss;
ci=ci_ss;
cr=cr_ss;
tau=tau_ss;
lambtilde=lambtilde_ss;
lamtau=lamtau_ss;
lami=lami_ss;
lams=lams_ss;
lamr=lamr_ss;
dd=0;
pop=1;
dcs=dcs_ss;
dns=dns_ss;
dci=dci_ss;
dni=dni_ss;
dw=dw_ss;
dlams=dlams_ss;
dlamtau=dlamtau_ss;
dlambtilde=dlambtilde_ss;
dlami=dlami_ss;
dlamr=dlamr_ss; 
dcr=dcr_ss;
dnr=dnr_ss;
drk=drk_ss;
rr=rr_ss;
Rb=Rb_ss; 
pie=pie_ss; 
mc=mc_ss; 
F=F_ss; 
Kf=Kf_ss; 
dF=1; 
dKf=1; 
pbreve=1; 
dpie=1;
yF=y_ss; 
kF=k_ss;
nF=n_ss;
wF=w_ss;
rkF=rk_ss;
xF=x_ss;
cF=c_ss;
sF=1;
iF=0; 
rF=0; 
nsF=ns_ss;
niF=ni_ss; 
nrF=nr_ss; 
csF=cs_ss;
ciF=ci_ss;
crF=cr_ss;
tauF=tau_ss;
lambtildeF=lambtilde_ss;
lamtauF=lamtau_ss;
lamiF=lami_ss;
lamsF=lams_ss;
lamrF=lamr_ss;
ddF=0;
popF=1;
dcsF=dcs_ss;
dnsF=dns_ss;
dciF=dci_ss;
dniF=dni_ss;
dwF=dw_ss;
dlamsF=dlams_ss;
dlamtauF=dlamtau_ss;
dlambtildeF=dlambtilde_ss;
dlamiF=dlami_ss;
dlamrF=dlamr_ss; 
dcrF=dcr_ss;
dnrF=dnr_ss;
drkF=drk_ss;
rrF=rr_ss;
RbF=Rb_ss; 
pieF=pie_ss; 
mcF=mc_ss;
FF=FF_ss; 
KfF=KfF_ss;
dFF=1; 
dKfF=1; 
pbreveF=1; 
dpieF=1;
@# include "commonVarSS.mod" //added by epi-mmb team
Interest=Rb_ss;
Inflation=pie_ss;
end;

//calculate residuals of dynamic equations with provided steady state
//and then invoke steady state computation followed by another check.
resid;
steady;
resid;
 
 
//set unknown type state variables to zero. Initialize initial seed 
M_.endo_histval=oo_.steady_state;
M_.endo_histval(strmatch('i',M_.endo_names,'exact'))=i_ini;
M_.endo_histval(strmatch('s',M_.endo_names,'exact'))=1-i_ini;
M_.endo_histval(strmatch('iF',M_.endo_names,'exact'))=i_ini;
M_.endo_histval(strmatch('sF',M_.endo_names,'exact'))=1-i_ini;


//Back up on price stickiness parameter (if nonzero); use homotopy below
//to impose final desired value.
M_.params(strmatch('xi',M_.param_names,'exact'))=xi_final/1.1;


//solve and simulate model
options_.slowc=1;
options_.simul.maxit=100;
simul(periods=500,stack_solve_algo=0);
 

//Homotopy for pi3 and xi parameters (otherwise no solution if you set them to final values right away)
pi3_final_steps_vec=[pi3_final/3:0.02:pi3_final,pi3_final];
for pi3_final_step=pi3_final_steps_vec
    pi3_final_step
    M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final_step;
    [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
end

xi_final_steps_vec=[M_.params(strmatch('xi',M_.param_names,'exact')):0.005:xi_final];
for xi_final_step=xi_final_steps_vec
    xi_final_step
    M_.params(strmatch('xi',M_.param_names,'exact'))=xi_final_step;
    [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
end


//Simulate model with final (desired, calibrated) values for pi1,pi2,pi3 and xi 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final;
M_.params(strmatch('xi',M_.param_names,'exact'))=xi_final;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above    
   


//plot results
//plot_agg_results;
//plot_by_type_results;


//Save results
//save all_results   
@# include "saveResults.mod"   //added by epi-mmb team