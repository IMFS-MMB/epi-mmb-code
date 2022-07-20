//This file solves and simulates the model developed in Ascari, Colciago 
//and Silvestrini (2021). Based on Eichenbaum, Rebelo and Trabandt (2020),
//'Epidemics in the Neoclassical and New Keynesian Models'.
// 2 Sectors, Nominal Frictions, Roundabout

//This code was written and run with Matlab 2018b and DYNARE 4.4.3. The nonlinear
//model is solved and simulated with Dynare's nonlinear deterministic solver.

////////////////////////////////////////////////////////////
//IMPORTANT: You must use DYNARE 4.4.3 to run this code.////
// addpath c:\dynare\4.4.3\matlab ////////////////////////// 
////////////////////////////////////////////////////////////

// Care, since it is 4.4.3, no preprocessor can be used, so change by hand!!

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var fe1 fe2 C1 C2 ls Ld w mutilde lambda lambdat lambdai lambdas 
lambdad C zc1 zc2 N1 N2 rho1 rho2 e1 e2 R pi wtilde f1 f2 X Ne1 
Ne2 stilde Y tt ss ii dd rr no1 no2 GDP dC1 dls dLd dw dmutilde 
dlambda dlambdai dlambdad dC de1 de2 dpi df1 df2 dzc1 dzc2 
dlambdat dfe1 dfe2 Ztilde harm_av omega;
var Consumption Labour Output Susceptibles Infected Recovered Deaths;     // added by epi-mmb team
var Interest Inflation Investment;


%----------------------------------------------------------------
% 2. Defining parameters
%----------------------------------------------------------------

parameters pibar Ybar Rbar zmin Zbar nu phi betta eta chi kappa 
theta delta thetaw alphatilde alpha alphaw tr_pi tr_y tr_smut tau 
fx pi1 pi2 pi3 pir pid ud i_ini psi0 psi1 gamma;

%----------------------------------------------------------------
% 3. calibration
%----------------------------------------------------------------

chi=0.5;

eta=1.5;

theta=3.8;   %elasticity of sub

kappa=6;     %tail Pareto (to have HHI)

zmin=1;   %min Pareto

Zbar=1;    %aggregate productivity

fx=0.47;     %fixed costs production

thetaw=4;    %elasticity labor inputs

delta=0.00211;   %exit rate

betta=0.98^(1/52);    %Weekly household discount factor

phi=4;   %Frisch elasticity labor supply

tau=(kappa/(kappa-(theta-1)))^(1/(theta-1));   %ztilde=tau*zc

alphaw=0;   %fraction wage bill in advance

nu=1;   %disutility labor multiplier

alpha=1/3;   %cobb-douglas

alphatilde=0.98;  %calvo union wage resetting, change price once per year

tr_pi=1.5;  %Taylor rule

tr_smut=0.8;  %Taylor rule

tr_y=0.5/52;  %Taylor rule

pid=7*0.002/14;    %Weekly probability of dying

pir=7*1/14-pid;    %Weekly probability of recovering 

ud=10;    %%%%%%% Psycological cost of death

//i_ini=1e-3;        %Initial seed of infection
load inf_ini
i_ini=helper; 

psi0=1;

psi1=1000;

gamma=1.5;

format long;

pi3_shr_target=2/3;                   %share of T_0 jump due general infections

pi1_shr_target=(1-pi3_shr_target)/2;  %share of T_0 jump due to consumption-based infections

pi2_shr_target=(1-pi3_shr_target)/2;  %Share of T_0 jump due to work-based infections

RplusD_target=0.6;                   %total share of people infected and then either recovered or dead after epidemic


%--------------------------------------------------------------------------------------------------------------------------------
% Steady State: to denote ss values attach "bar" at the end of the name of the variable. Es. Z=>Zbar  
%--------------------------------------------------------------------------------------------------------------------------------

%Load Numerical part: 11 variables

load SS_for_dynare y;

Cbar=y(1);

C1bar=y(2);

lambdabar=y(3);

zc1bar=y(4);

N1bar=y(5);

zc2bar=y(6);

N2bar=y(7);

wbar=y(8);

Xbar=y(9);

Ybar=y(10);

C2bar=y(11);


%recursive steady state 

rho1bar=(theta/(theta-1))*(1/(Zbar*tau*zc1bar))*(((alphaw/betta+1-alphaw)*wbar)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

rho2bar=(theta/(theta-1))*(1/(Zbar*tau*zc2bar))*(((alphaw/betta+1-alphaw)*wbar)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

e1bar=(1/theta)*chi*rho1bar^(1-eta)*((zmin/zc1bar)^kappa*N1bar)^((theta-eta)/(1-theta))*(lambdabar*Cbar)^(-eta)*Ybar-fx;

e2bar=(1/theta)*(1-chi)*rho2bar^(1-eta)*((zmin/zc2bar)^kappa*N2bar)^((theta-eta)/(1-theta))*(lambdabar*Cbar)^(-eta)*Ybar-fx;

mutildebar=(thetaw)/(thetaw-1);

lsbar=((lambdabar*wbar)/(nu*mutildebar))^(1/phi);

fe1bar=psi0+psi1*(delta/(1-delta)*N1bar)^gamma;

fe2bar=psi0+psi1*(delta/(1-delta)*N2bar)^gamma;

ddbar=0;

iibar=0;

ssbar=1;

ttbar=0;

rrbar=0;

lambdasbar=0;

lambdadbar=(betta*log(Cbar)-betta*nu*(lsbar^(1+phi))/(1+phi)+betta*ud-betta*lambdabar*Cbar+betta*lambdabar*wbar/mutildebar*lsbar)/(1-betta);

lambdaibar=betta*pid*lambdadbar/(1-betta*(1-pir)+betta*pid);

lambdatbar=lambdaibar;

Ldbar=lsbar;

Rbar=1/betta;

Ne1bar=N1bar*(delta/(1-delta));

Ne2bar=N2bar*(delta/(1-delta));

stildebar=1;

pibar=1;

wtildebar=wbar;

f1bar=Ldbar*lambdabar*wbar^thetaw/(1-alphatilde*betta);

f2bar=Ldbar*wbar^thetaw*nu*lsbar^phi/(1-alphatilde*betta);

GDPbar=Ybar-Xbar;

no1bar=(zmin/zc1bar)^kappa*N1bar;

no2bar=(zmin/zc2bar)^kappa*N2bar;

Ztildebar=(chi*(no1bar/(no1bar+no2bar))^(1/(1-theta))*(1/(zc1bar*tau))+(1-chi)*(no2bar/(no1bar+no2bar))^(1/(1-theta))*(1/(zc2bar*tau)))^(-1);

harm_avbar=2*(1/(tau*zc1bar)+1/(tau*zc2bar))^(-1);

omegabar=no1bar/(no1bar+no2bar);

Consumption_ss = Cbar;            //added by epi-mmb team
Labour_ss = Ldbar;
Output_ss = Ybar;
Susceptibles_ss = ssbar; 
Infected_ss = iibar;
Recovered_ss = rrbar;
Deaths_ss = ddbar;

//calibrate the pi's in the transmission (tau) - function
go_calibrate_pi;
 

//final numbers for pi1,pi2,pi3 and will be imposed below using homotopy
pi1_final=pi1;
pi2_final=pi2;
pi3_final=pi3;

//put scaled down values of pi1,pi2,pi3 into Dynare M_. structure
//if you put the final numbers from the get go, no solution will be found
//use smaller numbers first then compute a solution. Then use the solution 
//as an initial guess for slightly larger values. Proceed using this homotopy 
//until final values are imposed (all below, after model block). 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final/3;


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
//equilibrium equations
//////////////////////////////////////////////////////

tt=pi1*ss(-1)*ii(-1)*C1^2/(1-dd(-1))+ss(-1)*ii(-1)*pi2*ls*Ld+pi3*ss(-1)*ii(-1);

ii=ii(-1)+tt-(pir+pid)*ii(-1);

ss=ss(-1)-tt;

rr=rr(-1)+pir*ii(-1);

dd=dd(-1)+pid*ii(-1);

nu*ls^phi=lambda*w/mutilde-lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi2*Ld;

lambdat=lambdai-lambdas;

lambdai=-betta*(-lambdai*dlambdai(+1)*(1-pir)+pid*(lambdai*dlambdai(+1)-lambdad*dlambdad(+1)));

lambdas=-betta*(lambdat*dlambdat(+1)*(-ii/(1-dd)*pi1*(C1*dC1(+1))^2-ii*pi2*ls*dls(+1)*Ld*dLd(+1)-pi3*ii)-lambdas(+1));

lambdad=-betta*(-log(C*dC(+1)/(1-dd))+nu*(ls*dls(+1))^(1+phi)/(1+phi)-ud+lambda*dlambda(+1)*C*dC(+1)/(1-dd)-lambda*dlambda(+1)*w*dw(+1)/(mutilde*dmutilde(+1))*ls*dls(+1)-lambdad*dlambdad(+1));

C=C1^(1/(1-eta))*(1-dd(-1))^(eta/(eta-1))*chi^(1/(eta-1))*(lambda*rho1*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi1*C1)^(eta/(1-eta));

C2=(1-chi)*(1-dd(-1))^eta*(lambda*rho2*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta)))^(-eta)*C^(1-eta);

C=(chi^(1/eta)*C1^((eta-1)/eta)+(1-chi)^(1/eta)*C2^((eta-1)/eta))^(eta/(eta-1));

fe1=psi0+psi1*(Ne1)^gamma;

fe2=psi0+psi1*(Ne2)^gamma;

fe1=betta*(1-delta)*(dlambda(+1))*(fe1*dfe1(+1)+((zmin/(zc1*dzc1(+1)))^kappa)*e1*de1(+1));

fe2=betta*(1-delta)*(dlambda(+1))*(fe2*dfe2(+1)+((zmin/(zc2*dzc2(+1)))^kappa)*e2*de2(+1));

1=betta*(dlambda(+1))*R/(dpi(+1)*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unions

wtilde=thetaw/(thetaw-1)*f2/f1;

f1=Ld*w^(thetaw)*lambda+alphatilde*betta*(pi*dpi(+1))^(thetaw-1)*f1*df1(+1);

f2=Ld*w^(thetaw)*(nu*ls^(phi)+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi2*Ld)+alphatilde*betta*(pi*dpi(+1))^(thetaw)*f2*df2(+1);

%%%%%%%%%%%%%%%%%%%%%%%%% Firms 

rho1=theta/(theta-1)*(1/(Zbar*tau*zc1))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

rho2=theta/(theta-1)*(1/(Zbar*tau*zc2))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

e1=(1/theta)*chi*rho1*((zmin/zc1)^kappa*N1(-1))^(theta/(1-theta))*(lambda*rho1*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi1*C1)^(-eta)*(C/(1-dd(-1)))^(-eta)*Y-fx;

e2=(1/theta)*(1-chi)*rho2*((zmin/zc2)^kappa*N2(-1))^(theta/(1-theta))*(lambda*rho2*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta)))^(-eta)*(C/(1-dd(-1)))^(-eta)*Y-fx;

%%%%%%%%%%%%%%%%%%% Entry and Exit

N1=(1-delta)*(N1(-1)+Ne1);

N2=(1-delta)*(N2(-1)+Ne2);

no1=(zmin/zc1)^kappa*N1(-1);

no2=(zmin/zc2)^kappa*N2(-1);

zc1=(theta^(theta/(theta-1)))/(theta-1)*(1/Zbar)*(fx/(rho1^theta*((zmin/zc1)^kappa*N1(-1))^(theta/(1-theta))*(lambda*rho1*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi1*C1)^(-eta)*(C/(1-dd(-1)))^(-eta)*Y))^(1/(theta-1))*chi^(1/(1-theta))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

zc2=(theta^(theta/(theta-1)))/(theta-1)*(1/Zbar)*(fx/(rho2^theta*((zmin/zc2)^kappa*N2(-1))^(theta/(1-theta))*(lambda*rho2*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta)))^(-eta)*(C/(1-dd(-1)))^(-eta)*Y))^(1/(theta-1))*(1-chi)^(1/(1-theta))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

%%%%%%%%%%%%%%%%%%%%% Taylor rule

R/Rbar=((pi/pibar)^tr_pi*(Y/Ybar)^tr_y)^(1-tr_smut)*(R(-1)/Rbar)^tr_smut;

%%%%%%%%%%%%%%%%% Aggregation

C+Ne1*fe1+Ne2*fe2=(alphaw*R+1-alphaw)*w*Ld+e1*N1(-1)*(zmin/zc1)^kappa+e2*N2(-1)*(zmin/zc2)^kappa;

C^eta=(1-dd(-1))^eta*(chi*rho1*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))*(lambda*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))*rho1+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi1*C1)^(-eta)+(1-chi)*lambda^(-eta)*(rho2*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta)))^(1-eta));

Y=C+X+fe1*Ne1+fe2*Ne2+fx*(N1(-1)*(zmin/zc1)^kappa+N2(-1)*(zmin/zc2)^kappa);

X=chi*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))/(Zbar*tau*zc1)*(alpha*w*(alphaw*R+1-alphaw)/(1-alpha))^(1-alpha)*(lambda*rho1*((zmin/zc1)^kappa*N1(-1))^(1/(1-theta))+lambdat*ss(-1)*ii(-1)/(1-dd(-1))*pi1*C1)^(-eta)*(C/(1-dd(-1)))^(-eta)*Y+(1-chi)*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta))/(Zbar*tau*zc2)*(alpha*w*(alphaw*R+1-alphaw)/(1-alpha))^(1-alpha)*(lambda*rho2*((zmin/zc2)^kappa*N2(-1))^(1/(1-theta)))^(-eta)*(C/(1-dd(-1)))^(-eta)*Y;

w^(1-thetaw)=(1-alphatilde)*wtilde^(1-thetaw)+alphatilde*(w(-1)/pi)^(1-thetaw);

(1-dd(-1))*ls=stilde*Ld;

stilde=(1-alphatilde)*(wtilde/w)^(-thetaw)+alphatilde*(w(-1)/w)^(-thetaw)*pi^(thetaw)*stilde(-1);

GDP=Y-X;

%%%%%%%%%%%% Productivity

Ztilde=(chi*(no1/(no1+no2))^(1/(1-theta))*(no1^(1/(1-theta))*rho1)^(-eta)*(1/(zc1*tau))+(1-chi)*(no2/(no1+no2))^(1/(1-theta))*(no2^(1/(1-theta))*rho2)^(-eta)*(1/(zc2*tau)))^(-1);

harm_av=2*(1/(tau*zc1)+1/(tau*zc2))^(-1);

omega=no1/(no1+no2);



//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.

dC1=C1/C1(-1); 
dls=ls/ls(-1); 
dLd=Ld/Ld(-1);  
dw=w/w(-1); 
dmutilde=mutilde/mutilde(-1); 
dlambda=lambda/lambda(-1); 
dlambdat=lambdat/lambdat(-1); 
dlambdai=lambdai/lambdai(-1); 
dlambdad=lambdad/lambdad(-1); 
dC=C/C(-1);  
de1=e1/e1(-1); 
de2=e2/e2(-1);  
dfe1=fe1/fe1(-1); 
dfe2=fe2/fe2(-1);  
dpi=pi/pi(-1); 
df1=f1/f1(-1);
df2=f2/f2(-1); 
dzc1=zc1/zc1(-1);
dzc2=zc2/zc2(-1); 

// ===============
// epi-mmb varialbes
Consumption = C;
Labour = Ld;
Output = Y;
Susceptibles = ss; 
Infected = ii;
Recovered = rr;
Deaths = dd; 
Interest=R;
Inflation=pi;
Investment=0;
end;


//initial (pre-infection) steady state used in simulations 
initval;

C1=C1bar; 
C2=C2bar; 
ls=lsbar; 
Ld=Ldbar;  
w=wbar; 
mutilde=mutildebar; 
lambda=lambdabar;
lambdat=lambdatbar; 
lambdai=lambdaibar; 
lambdas=lambdasbar; 
lambdad=lambdadbar; 
C=Cbar;  
zc1=zc1bar; 
zc2=zc2bar; 
N1=N1bar; 
N2=N2bar; 
no1=N1bar*(zmin/zc1bar)^kappa; 
no2=N2bar*(zmin/zc2bar)^kappa; 
wtilde=wtildebar; 
rho1=rho1bar; 
rho2=rho2bar; 
e1=e1bar; 
e2=e2bar;
fe1=fe1bar; 
fe2=fe2bar; 
R=Rbar; 
pi=pibar; 
f1=f1bar;
f2=f2bar; 
X=Xbar; 
Ne1=Ne1bar; 
Ne2=Ne2bar; 
stilde=stildebar; 
Y=Ybar; 
ss=ssbar;
ii=iibar; 
dd=ddbar; 
tt=ttbar;
rr=rrbar;
GDP=GDPbar;
dC1=1; 
dls=1; 
dLd=1;  
dw=1; 
dmutilde=1; 
dlambda=1; 
dlambdat=1; 
dlambdai=1; 
dlambdad=1; 
dC=1;   
de1=1; 
de2=1; 
dfe1=1; 
dfe2=1; 
dpi=1; 
df1=1;
df2=1; 
dzc1=1;
dzc2=1; 
Ztilde=Ztildebar;
harm_av=harm_avbar;
omega=omegabar;
// steady state for common variables
Consumption = Consumption_ss;
Labour = Labour_ss;
Output = Output_ss;
Susceptibles = Susceptibles_ss; 
Infected = Infected_ss;
Recovered = Recovered_ss;
Deaths = Deaths_ss;
Interest=Rbar;
Inflation=pibar;
Investment=0;
end;

//calculate residuals of dynamic equations with provided steady state
//and then invoke steady state computation followed by another check.
resid;
steady;
resid;
 
 
//set unknown type state variables to zero. Initialize initial seed 
M_.endo_histval=oo_.steady_state;
M_.endo_histval(strmatch('ii',M_.endo_names,'exact'))=i_ini;
M_.endo_histval(strmatch('ss',M_.endo_names,'exact'))=1-i_ini;



//solve and simulate model
options_.slowc=1;
options_.simul.maxit=100;
simul(periods=520,stack_solve_algo=0);
 

//Homotopy for pi3 parameters (otherwise no solution if you set them to final values right away)
pi3_final_steps_vec=[pi3_final/3:0.02:pi3_final,pi3_final];
for pi3_final_step=pi3_final_steps_vec
    pi3_final_step
    M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final_step;
    [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
end



//Simulate model with final (desired, calibrated) values for pi1,pi2,pi3
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above       


//plot results
//plot_agg_results;


//Save results
//save all_results   
// save results
results.oo_ = oo_ ;
results.M_ = M_;
save simulated_results results;
