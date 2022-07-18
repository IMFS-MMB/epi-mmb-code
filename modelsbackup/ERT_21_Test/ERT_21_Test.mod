//This file solves and simulates the model developed in Eichenbaum,
//Rebelo and Trabandt (2020),'The Macroeconomics of Testing and Quarantining'

//This code was written and run with Matlab 2016b and DYNARE 4.6.1. The nonlinear
//model is solved and simulated with Dynare's nonlinear deterministic solver

////////////////////////////////////////////////////////////
//IMPORTANT: You must use DYNARE 4.6.1 to run this code.//// 
////////////////////////////////////////////////////////////

//We solve the fully nonlinear model. Given the strong nonlinearities of the 
//model due to the epidemic dynamics we use homotopy for a number of 
//parameters, i.e. we increase parameter values stepwise until we impose
//their final desired value. Homotopy calculations can take a while. 

//This code only produces plots for disaggregated model variables. Aggregate 
//variables are plotted using the fig1_...through fig8_... files in the 
//main directory of this code kit.

//Mathias Trabandt, mathias.trabandt@gmail.com


//Endogenous model variables 
var 

p_s_a p_i_a p_r_a Uu Us Ui Ur cu nu cs ns ci ni cr nr
tauu taus Iu Ik Su Sk Ru Rk Tu Tk Du Dk pid_endo Gamma Gammai
aggC aggN Pop Tests lambu lamtauu lamiu lamsu lambs lamtaus 
lambi lambr lamru masku masks maski maskr Mstar dUi dUu dpid_endo 
dUs dUr dlamiu dp_i_a dlamsu dp_s_a dlamru dtauu dp_r_a;
@#include "commonVar.mod"    // added by epi-mmb team
///////////////////////////////////////////////////////////////////////////
//Important: in several places in the code below, the word 'containment' 
//as used in parameters/switches/variables refers to the word 'quarantines' 
//as used in the manuscript.
//Also, the word 'mask' in the code below refers to 'NPI' as used in the 
//manuscript.
///////////////////////////////////////////////////////////////////////////

//Shocks (aka exogenous 0/1 switches for containment simulations) 
varexo cont_switch ni_cont_val Gammai_cont_val;

//Define parameters
parameters pi1 pi2 pi3 pir pid betta i_ini A theta 
alfa inc_target kappa pis b do_strict gam;

//Initialize parameters

stepps=100; //homotopy steps for pi1, pi2 and pi3

betta=0.96^(1/52);   //Weekly household discount factor
days=14;
pid=7*0.0035/days;      //Weekly probability of dying
pir=7*1/days-pid;      //Weekly probability of recovering

sumpidpir=pid+pir   //print sum of pid and pir

gam=1;%1e7;     //utility parameter masks

pis=1/104;        //Weekly probability of re-infection

kappa=0;           //slope endog. mortality rate

alfa=0;              //share of people tested in pre-infection steady state
                     //the final value of alfa will be imposed below using homotopy 
                     //DO NOT CHANGE alfa here as it must be kept at zero for the 
                     //pre-infection steady state (no testing).

b=0;%-4.05           //constant in utility function (for recalibration of value of life)

load inf_ini
i_ini=helper; 
//i_ini=0.00066;         //Initial seed for infections

n_target=28;         //Calibration target for weekly hours
inc_target=58000/52; //Calibration target for weekly income

//Calibration targets for shares of pi-terms in T-function in SIR model
pi3_shr_target=1/6;                   //share of T_0 jump due general infections
pi1_shr_target=(1-pi3_shr_target)/2;  //share of T_0 jump due to consumption-based infections
pi2_shr_target=(1-pi3_shr_target)/2;  //share of T_0 jump due to work-based infections
RplusD_target=0.83; 

//Per switches below, testing, containment and sensitivity can be activated
//Recall that smart containment is containment and testing, 
//i.e. set do_containment =  1; and do_set_alfa_nonzero=1;

do_containment=1;      //if =0, no containment. if =1 containment

do_strict=0;           //only affects results if do_containment = 1;
                       //if do_strict = 0, containment is smart. 
                       //if do_strict = 1, containment is strict.

do_set_alfa_nonzero=1; //if =1 then alfa=0.02/4. Otherwise alfa is kept at zero.

do_sensitivity=0;      //if =1, sensitivity analysis wrt alfa



//calibrate theta and A using pre-infection steady states
theta=n_target^(-2);
A=inc_target/n_target;

//calibrate the pi's in T-function
go_calibrate;
 
//final numbers for pi1,pi2,pi3 will be imposed below using homotopy
pi1_final=pi1;
pi2_final=pi2;
pi3_final=pi3;
pis_final=pis;

//put scaled down values of pi1,pi2,pi3 into Dynare M_. structure
//if you put the final numbers from the get go, no solution will be found
//use smaller numbers first then compute a solution. Then use the solution 
//as an initial guess for slightly larger values. Proceed using this homotopy 
//until final values are imposed (all below, after model block). 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final/1000;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final/1000;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final/100;
M_.params(strmatch('pis',M_.param_names,'exact'))=pis_final*0; 


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
//in the pre- and terminal steady states is the same. Please keep in mind that for variables 
//which have a pre-infection steady state that is non-zero we replace X(+1) by X*dX(+1) 
//and define dX=X/X(-1). For variables which have a zero pre-infection steady state we 
//replace X(+1) by ( X+dX(+1) ) and define dX=X-X(-1). 


///////////////////////////////////////////////////////////// 
//equilibrium equations for people with unknown health status
/////////////////////////////////////////////////////////////  
//Auxilliary variable: derivative of Uu wrt p_i_a
#dUu_dp_i_a=alfa*betta*(1-pid_endo*p_i_a(-1))^2*Ui*dUi(+1)
-pid_endo*dpid_endo(+1)*(1-alfa)^2*betta^2*(1-pid_endo*p_i_a(-1))*Uu*dUu(+1)*dUu(+2)
-pid_endo*dpid_endo(+1)*alfa*(1-alfa)*betta^2*(1-pid_endo*p_i_a(-1))*
( (p_s_a+dp_s_a(+1))*Us*dUs(+1)*dUs(+2) + (p_i_a+dp_i_a(+1))*Ui*dUi(+1)*dUi(+2) 
+ (p_r_a+dp_r_a(+1))*Ur*dUr(+1)*dUr(+2) );

//Auxilliary variable: derivative of Uu wrt p_s_a
#dUu_dp_s_a=alfa*betta*(1-pid_endo*p_i_a(-1))*Us*dUs(+1);

//Auxilliary variable: derivative of Uu wrt p_r_a
#dUu_dp_r_a=alfa*betta*(1-pid_endo*p_i_a(-1))*Ur*dUr(+1);

//Utility of people with unknown health status
Uu=b+log(cu)-theta/2*nu^2-gam/2*masku^2+(1-alfa)*betta*(1-pid_endo*p_i_a(-1))*Uu*dUu(+1)
+alfa*betta*(1-pid_endo*p_i_a(-1))*( p_s_a*Us*dUs(+1) + p_i_a*Ui*dUi(+1) + p_r_a*Ur*dUr(+1) );

//Transmission function for people with unknown health status
//Per switch containment is imposed -- either smart or strict for known types. 
tauu=(1-masku)*(pi1*cu*(Iu(-1)*cu*(1-masku)+Ik(-1)*(1-cont_switch)*ci*(1-maski))+pi2*nu*(Iu(-1)*nu*(1-masku)+Ik(-1)*(1-cont_switch)*ni*(1-maski))
+pi3*(Iu(-1)*(1-masku)+(1-cont_switch*do_strict)*Ik(-1)*(1-maski)));

//Budget for people with unknown health status
cu=A*nu+Gamma;

//Transition function probability of being infected
p_i_a*(1-pid_endo*p_i_a(-1))=tauu*p_s_a(-1)+(1-pir-pid_endo)*p_i_a(-1);

//Transition function probability of being susceptible
p_s_a*(1-pid_endo*p_i_a(-1))=p_s_a(-1)*(1-tauu)+pis*p_r_a(-1);

//Transition function probability of being recovered
p_r_a*(1-pid_endo*p_i_a(-1))=p_r_a(-1)*(1-pis)+pir*p_i_a(-1);

//FOC consumption, unknown health status  
1/cu-lambu+lamtauu*(1-masku)*pi1*(Iu(-1)*cu*(1-masku)+Ik(-1)*(1-cont_switch)*ci*(1-maski))=0;

//FOC labor, unknown health status
-theta*nu+lambu*A+lamtauu*(1-masku)*pi2*(Iu(-1)*nu*(1-masku)+Ik(-1)*(1-cont_switch)*ni*(1-maski))=0;

//FOC masks, unknown health status
-gam*masku-lamtauu*tauu/(1-masku)=0;

//FOC tauu
-lamtauu+lamiu*p_s_a(-1)-lamsu*p_s_a(-1)=0;

//FOC p_i_a
dUu_dp_i_a*1/(1-pid_endo*p_i_a(-1))-lamiu+betta*(lamiu+dlamiu(+1))*pid_endo*dpid_endo(+1)*(p_i_a+dp_i_a(+1))
+betta*(lamiu+dlamiu(+1))*(1-pir-pid_endo*dpid_endo(+1))
+betta*(lamsu+dlamsu(+1))*pid_endo*dpid_endo(+1)*(p_s_a+dp_s_a(+1))
+(lamru+dlamru(+1))*betta*pid_endo*dpid_endo(+1)*(p_r_a+dp_r_a(+1))+(lamru+dlamru(+1))*betta*pir;

//FOC p_s_a
dUu_dp_s_a*1/(1-pid_endo*p_i_a(-1))+betta*(lamiu+dlamiu(+1))*(tauu+dtauu(+1))-lamsu
+betta*(lamsu+dlamsu(+1))*(1-(tauu+dtauu(+1)))=0; 

//FOC p_r_a
dUu_dp_r_a*1/(1-pid_endo*p_i_a(-1))+betta*(lamsu+dlamsu(+1))*pis-lamru
+betta*(lamru+dlamru(+1))*(1-pis)=0; 



/////////////////////////////////////////////////////////////////////////// 
//equilibrium equations for people with known health status (after testing)
/////////////////////////////////////////////////////////////////////////// 
//Budget susceptibles
cs=A*ns+Gamma;

//Budget infected
//Note that Gammai=Gamma without containment. 
//With containment, Gammi is set to desired value, financed by Gamma of all other agents. 
ci=A*ni+Gammai;  

//Budget recovered
cr=A*nr+Gamma;

//Utility susceptibles
Us=b+log(cs)-theta/2*ns^2-gam/2*masks^2+betta*( (1-taus)*Us*dUs(+1) + taus*Ui*dUi(+1) );

//Probability of infection, susceptibles
//Per switch containment is imposed -- either smart or strict for known types. 
taus=(1-masks)*(pi1*cs*(Iu(-1)*cu*(1-masku)+Ik(-1)*(1-cont_switch)*ci*(1-maski))
+pi2*ns*(Iu(-1)*nu*(1-masku)+Ik(-1)*(1-cont_switch)*ni*(1-maski))
+pi3*(Iu(-1)*(1-masku)+(1-cont_switch*do_strict)*Ik(-1)*(1-maski)));

//FOC consumption, susceptibles
1/cs-lambs+lamtaus*(1-masks)*pi1*(Iu(-1)*cu*(1-masku)+Ik(-1)*(1-cont_switch)*ci*(1-maski))=0;

//FOC labor, susceptibles
-theta*ns+lambs*A+lamtaus*(1-masks)*pi2*(Iu(-1)*nu*(1-masku)+Ik(-1)*(1-cont_switch)*ni*(1-maski))=0;

//FOC masks, susceptibles
-gam*masks-lamtaus*taus/(1-masks)=0;

//FOC taus, susceptibles
betta*(Ui*dUi(+1)-Us*dUs(+1))-lamtaus=0;

//Utility infected
Ui=b+log(ci)-theta/2*ni^2-gam/2*maski^2+betta*( (1-pir-pid_endo)*Ui*dUi(+1) + pir*Ur*dUr(+1) );

//FOC consumption, infected
//note that with containment, equation pins down lambi only. 
1/ci-lambi=0;

//FOC labor, infected
//per switch set ni to desired value under containment
(1-cont_switch)*(-theta*ni+lambi*A)+cont_switch*(ni-ni_cont_val)=0;

//FOC masks, infected
(1-cont_switch)*(-gam*maski)+cont_switch*(maski-masks)=0;

//Utility recovered
Ur=b+log(cr)-theta/2*nr^2-gam/2*maskr^2+betta*(1-pis)*Ur*dUr(+1)+betta*pis*Us*dUs(+1);

//FOC consumption, recovered
1/cr-lambr=0;

//FOC labor, recovered
-theta*nr+lambr*A=0;

//FOC masks, recovered
(1-cont_switch)*(-gam*maskr)+cont_switch*(maskr-masks)=0;

//////////////////////
//Population Dynamics
/////////////////////
//Number of new infections, unknown health status
Tu=p_s_a(-1)*Mstar(-1)*tauu;

//Number of susceptibles, unknown health status
Su=p_s_a*Mstar;

//Number of infected, unknown health status
Iu=p_i_a*Mstar;

//Number of recovered, unknown health status
Ru=p_r_a*Mstar;

//Number of deaths, unknown health status
Du=Du(-1)+pid_endo*Iu(-1);

//Auxilliary variable to keep track of population of unknowns over time 
//Mstar is denoted Popu in manuscript
Mstar=Mstar(-1)*(1-pid_endo*p_i_a(-1))*(1-alfa);

//Number of new infections, known health status
//Per switch containment is imposed -- either smart or strict for known types.
Tk=pi1*Sk(-1)*cs*(1-masks)*(Iu(-1)*cu*(1-masku)+Ik(-1)*(1-cont_switch)*ci*(1-maski))
+pi2*Sk(-1)*ns*(1-masks)*(Iu(-1)*nu*(1-masku)+Ik(-1)*(1-cont_switch)*ni*(1-maski))
+pi3*Sk(-1)*(1-masks)*(Iu(-1)*(1-masku)+(1-cont_switch*do_strict)*Ik(-1)*(1-maski));

//Number of susceptibles, known health status
Sk=Sk(-1)-Tk+pis*Rk(-1)+alfa*(Su(-1)-Tu+pis*Ru(-1));

//Number of infected, known health status
Ik=Tk+(1-pir-pid_endo)*Ik(-1)+alfa*(Tu+(1-pir-pid_endo)*Iu(-1));

//Number of recovered, known health status
Rk=Rk(-1)+pir*Ik(-1)-pis*Rk(-1)+alfa*(Ru(-1)+pir*Iu(-1)-pis*Ru(-1));

//Number of deaths, known health status
Dk=Dk(-1)+pid_endo*Ik(-1);

//Total population
Pop=Pop(-1)-pid_endo*(Iu(-1)+Ik(-1));

//Medical preparedness -- endogenous mortality rate
pid_endo=pid+kappa*(Ik(-1)+Iu(-1))^2;

//Number of tests
Tests=Sk(-1)+Ik(-1)+alfa*(Su(-1)+(1-pid_endo)*Iu(-1)+Ru(-1));

//Government budget contstraint
0=Gammai*Ik(-1)+Gamma*(Sk(-1)+Rk(-1)+Su(-1)+Iu(-1)+Ru(-1));

//Transfer rule for containment
//When containment, then Gammai for known infected is set to desired value. 
//This controls how much consumption infected get under containment.
//Without containment, Gammai=Gamma
(1-cont_switch)*(Gammai-Gamma)+cont_switch*(Gammai-Gammai_cont_val); 

//Aggregate consumption
aggC=Sk(-1)*cs+Ik(-1)*ci+Rk(-1)*cr+(Su(-1)+Iu(-1)+Ru(-1))*cu;

//Aggregate labor
aggN=Sk(-1)*ns+Ik(-1)*ni+Rk(-1)*nr+(Su(-1)+Iu(-1)+Ru(-1))*nu;

@# include "commonVarEq.mod"  // added by epi-mmb team

//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dUi=Ui/Ui(-1);
dUu=Uu/Uu(-1);
dpid_endo=pid_endo/pid_endo(-1);
dUs=Us/Us(-1);
dUr=Ur/Ur(-1);

//Auxilliary variables: absolute rates of changes of variables with zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dlamiu=lamiu-lamiu(-1); 
dp_i_a=p_i_a-p_i_a(-1); 
dlamsu=lamsu-lamsu(-1); 
dp_s_a=p_s_a-p_s_a(-1);
dlamru=lamru-lamru(-1); 
dtauu=tauu-tauu(-1);
dp_r_a=p_r_a-p_r_a(-1);
end; 


//Initial values for endogenous model variables
cs_SS=inc_target;
ci_SS=inc_target;
cr_SS=inc_target;
cu_SS=inc_target;
ns_SS=n_target;
ni_SS=n_target;
nr_SS=n_target;
nu_SS=n_target;
Us_SS=1/(1-betta)*(b+log(cs_SS)-theta/2*ns_SS^2);
Ur_SS=1/(1-betta*(1-pis))*(b+log(cs_SS)-theta/2*ns_SS^2+betta*pis*Us_SS);
Ui_SS=1/(1-betta*(1-pir-pid))*(b+log(ci_SS)-theta/2*ni_SS^2 + betta*pir*Ur_SS);
Us_SSConsUnits=Us_SS*inc_target;
lamtaus_SS=betta*(Ui_SS-Us_SS);
Uu_SS=1/(1-(1-alfa)*betta)*(b+log(cu_SS)-theta/2*nu_SS^2+alfa*betta*Us_SS);
dUu_dp_i_a_SS=alfa*betta*(Ui_SS-Ur_SS)-pid*(1-alfa)^2*betta^2*Uu_SS-2*pid*alfa*(1-alfa)*betta^2*Us_SS;
lamsu_SS=1/(1-betta)*alfa*betta*(Us_SS-Ur_SS);
lamiu_SS=1/(1-betta*(1-pir-pid))*(dUu_dp_i_a_SS+betta*lamsu_SS*pid);
lamtauu_SS=lamiu_SS-lamsu_SS;

//initial (pre-infection, no testing) steady state used in simulations 
initval;
Ui=Ui_SS; 
Ur=Ur_SS; 
Us=Us_SS; 
Uu=Uu_SS;
lamtaus=lamtaus_SS;
p_i_a=0; 
p_s_a=1; 
tauu=0; 
taus=0; 
Iu=0; 
Ik=0; 
lambu=1/inc_target;
lambi=1/inc_target;
lambs=1/inc_target; 
lambr=1/inc_target;
lamtauu=lamtauu_SS;
lamiu=lamiu_SS;  
lamsu=lamsu_SS; 
Su=1; 
Mstar=1; 
Ru=0;  
Tu=0; 
Tk=0; 
Sk=0; 
Rk=0; 
Du=0;
Dk=0;
Pop=1;
cu=inc_target;  
nu=n_target; 
ci=inc_target; 
ni=n_target; 
cs=inc_target;  
ns=n_target;  
cr=inc_target;  
nr=n_target; 
aggC=inc_target; 
aggN=n_target;
Gamma=0; 
Gammai=0;
Tests=0;
pid_endo=pid;
p_r_a=0; 
lamru=0;
dUi=1;
dUu=1;
dpid_endo=1;
dUs=1;
dUr=1;
dlamiu=0; 
dp_i_a=0; 
dlamsu=0; 
dp_s_a=0;
dlamru=0; 
dtauu=0;
dp_r_a=0;

masku=0;
maski=0;
maskr=0;
masks=0;


@# include "commonVarSS.mod"  // added by epi-mmb team
end;


//check residuals of dynamic equations given initial steady state
resid;
steady;
resid;
 
 
 
//set initial seed of infection
M_.endo_histval=oo_.steady_state;
M_.endo_histval(strmatch('Iu',M_.endo_names,'exact'))=(i_ini);
M_.endo_histval(strmatch('Su',M_.endo_names,'exact'))=(1-i_ini);
M_.endo_histval(strmatch('p_i_a',M_.endo_names,'exact'))=(i_ini);
M_.endo_histval(strmatch('p_s_a',M_.endo_names,'exact'))=(1-i_ini);
%M_.endo_histval(strmatch('Ik',M_.endo_names,'exact'))=i_ini*0;
%M_.endo_histval(strmatch('Sk',M_.endo_names,'exact'))=1-i_ini*0;


//solve and simulate model
options_.slowc=1;
options_.simul.maxit=20;
options_.dynatol.f=1e-8;
options_.dynatol.x=1e-8;
simul(periods=500,stack_solve_algo=0);


 
//Containment; per switch do_strict either strict or smart containment. do_strict
//switch affects model equations directly. 
if do_containment==1
    //containment is implemented as follows:
    //-set ni=0
    //-set ci=inc_target, financed by Gammai transfers
    //-set ci and ni in ALL transmission functions and FOCS to zero, i.e. quarantine
    //-finance Gammai transfers by general lump sum taxes for everybody else 
    shocks;
    var cont_switch; //containment switch
    periods 1:500;    
    values 1;

    var ni_cont_val; //value of ni 
    periods 1:500;
    values 0;

    var Gammai_cont_val; //value of Gammai
    periods 1:500;
    values(inc_target);

    end;

    oo_=make_ex_(M_,options_,oo_);

    [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 
end
 

%[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
%dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 

 
//Homotopy for pi1 and pi2 parameters (otherwise no solution if you set them to final values right away)
pi3_final_steps_vec=linspace(pi3_final/1000,pi3_final/1,stepps);
pi2_final_steps_vec=linspace(pi2_final/1000,pi2_final/1,stepps);
pi1_final_steps_vec=linspace(pi1_final/100,pi1_final/1,stepps);
pis_final_steps_vec=linspace(pis_final/1000,pis_final/1,stepps);
for pi3_final_step=pi3_final_steps_vec
    M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
end
for pi2_final_step=pi2_final_steps_vec
    M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
end
for pi1_final_step=pi1_final_steps_vec
    M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess

end
for pis_final_step=pis_final_steps_vec
    M_.params(strmatch('pis',M_.param_names,'exact'))=pis_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    
end
     dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 
oo_endo_simul=oo_.endo_simul;
save oo_endo_simul_old oo_endo_simul; //save simulation results for later use as initial guess
 


//Simulate model with final (desired, calibrated) values for pi1,pi2,pi3 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pi3',M_.param_names,'exact'))=pi3_final;
M_.params(strmatch('pis',M_.param_names,'exact'))=pis_final;
load oo_endo_simul_old;
oo_.endo_simul=oo_endo_simul;
 [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 
  



//activate testing if desired using homotopy
if do_set_alfa_nonzero==1
    //Homotopy for alfa (otherwise no solution) 
    alfa_steps_vec=[0.0001:0.001/10:0.02/4,0.02/4];
    for alfa_step=alfa_steps_vec;
        M_.params(strmatch('alfa',M_.param_names,'exact'))=alfa_step;
        [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
    end
 dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 

end
 
@# include "saveResults.mod"  // added by epi-mmb team
