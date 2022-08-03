//This file solves and simulates the model developed in Eichenbaum,
//Rebelo and Trabandt (2021),'Inequality in Life and Death', IMF Economic Review

//This code was written and run with Matlab 2021a and DYNARE 4.6.1. The nonlinear
//model is solved and simulated with Dynare's nonlinear deterministic solver

////////////////////////////////////////////////////////////
//IMPORTANT: You must use DYNARE 4.6.1 to run this code.//// 
////////////////////////////////////////////////////////////

//We solve the fully nonlinear model. Given the strong nonlinearities of the 
//model due to the epidemic dynamics we use homotopy for a number of 
//parameters, i.e. we increase parameter values stepwise until we impose
//their final desired value. Homotopy calculations can take a while. 

//Mathias Trabandt, mathias.trabandt@gmail.com

//Endogenous model variables 
var consshift taushift Gaml Gamh inch incl incshrh Cl Ch Nl Nh NINCl NINCh C 
N P1 P2 bgov bhstar wl wh csl1 csh1 csl2 csh2 cil1 cih1 cil2 cih2 crl1 crh1 
crl2 crh2 nsl nsh nil nih nrl nrh taul tauh sl sh il ih rl rh dl dh Ul Uh U
residEuler diffbhstar profits lamtaul lamtauh lambl lambh 
lamsl lamsh lamil lamih lamrl lamrh dcsh1 dcsh2 dcih1 dcih2 dcil1 dcil2 dnsh 
dnih dwh dlamsh  dlamtauh dlambh dP1 dP2 dcsl1 dcsl2 dnsl dnil dwl dlamsl 
dlamtaul dlambl dlamih dlamrh dlamil dlamrl dcrh1 dcrh2 
dnrh dcrl1 dcrl2 dnrl dUh dUl taulSIRmod tauhSIRmod slSIRmod shSIRmod ilSIRmod 
ihSIRmod rlSIRmod rhSIRmod dlSIRmod dhSIRmod;
@#include "commonVar.mod"    // added by epi-mmb team

//Exogenous variables
varexo  bhstareps transhk taushk contain gov_switch;

//Define parameters
parameters pi1 pi2 Al Ah pil3 pih3 pi4 pidl pidh eta
cbar betta i_ini theta pirl pirh ml mh rstar 
do_sticky_wages P2_ss_numeraire_normalization;

//Initialize parameters
betta=0.98^(1/52);     //Weekly household discount factor
pidh=7*0.00357/14; //Weekly probability of dying
pidl=7*0.005/14;   //Weekly probability of dying

pirh=7*1/14-pidh;   //Weekly probability of recovering
pirl=7*1/14-pidl;  //Weekly probability of recovering

P2_ss_numeraire_normalization=1; //Numeraire

sh_ss=0.18;         //Share of high skilled people of total population
bhstar_ss=380000;   //Asset holdings high skilled
              
eta=0.5;            //Utility weight good c2
ml=6.5;             //constant in utility function high skilled
mh=6.5;             //constant in utility function low skilled

do_sticky_wages=1; //if =1, sticky wages; if =0 flexible wages

load inf_ini
i_ini=helper; 
//i_ini=0.001;       //Initial seed for infections

n_target=28;           //Calibration target for weekly hours
inc_target=58000/52;   //Calibration target for weekly income (labor and interest income)
winc_share_target=0.38;//Calibration target for wage income share (high skilled/total)

pop=330000000; //US population (total)


//Options
horzz=150;//periods to find evolution of assets before convergence; check residEuler -- must be small.
load_ini_guess_diffbhstar=1; //load existing path of assets as initial guess

wakeupWeek=13; //week number when people take epi into account; Standard SIR model before;
SIRstartWeek=1;//week number when epi starts in standard SIR model


//Solver options 
homotopy_steps=50;//homotopy steps for pi4 when solving the model
scale_pi4_homotopy=20; //scaling pi4 for homotopy

options_.slowc=1;
options_.simul.maxit=20;
options_.periods=175;


//Calibration targets for shares of pi-terms in T-function in SIR model
pi12_shr_target=1/12; %share of T_0 jump due to consumption-based infections
pi3_shr_target=5/12;  %share of T_0 jump due to work-based infections
pi4_shr_target=6/12;  %share of T_0 jump due general infections
RplusD_target=0.6;    %total share of people infected and then either recovered or dead after epidemic

pi1vs2factor=1.05;    %pi1 vs pi2 
pi3lvshfactor=20;     %pi3l vs pi3h

if (pi12_shr_target+pi3_shr_target+pi4_shr_target)~=1
    error('pi shares dont add up to one');
    return;
end

//Calibrate Pre-epidemic steady state
calib_steady;

theta=theta_calib;
Ah=Ah_calib;
Al=Al_calib;
cbar=cbar_calib;
csl1_ss=csl1_calib;
P1_ss=P1_calib;

sl_ss=1-sh_ss;
rstar=1/betta-1;
P2_ss=P2_ss_numeraire_normalization;

nsl_ss=Al/theta*(1-eta)/(csl1_ss-cbar);
csh1_ss=1/sh_ss*(Al*sl_ss*nsl_ss-sl_ss*csl1_ss);
nsh_ss=Ah/theta*(1-eta)/(csh1_ss-cbar)*P2_ss/P1_ss; 
    
wl_ss=P1_ss*Al;    
wh_ss=P2_ss*Ah;

csl2_ss=eta/(P2_ss/P1_ss*(1-eta)/(csl1_ss-cbar));
csh2_ss=eta/(P2_ss/P1_ss*(1-eta)/(csh1_ss-cbar));

C_ss=P1_ss*(sh_ss*csh1_ss+sl_ss*csl1_ss)+P2_ss*(sh_ss*csh2_ss+sl_ss*csl2_ss);
N_ss=sl_ss*nsl_ss+sh_ss*nsh_ss;
lambh_ss=1/P1_ss*(1-eta)/(csh1_ss-cbar);
lambl_ss=1/P1_ss*(1-eta)/(csl1_ss-cbar);

cih1_ss=csh1_ss; 
crh1_ss=cih1_ss;
cil1_ss=csl1_ss; 
crl1_ss=csl1_ss;

cih2_ss=csh2_ss; 
cil2_ss=csl2_ss; 
crh2_ss=csh2_ss;
crl2_ss=csl2_ss;

nih_ss=nsh_ss;
nil_ss=nsl_ss;
nrh_ss=nsh_ss;
nrl_ss=nsl_ss;

taul_ss=0;
tauh_ss=0;
ih_ss=0;
il_ss=0;
rh_ss=0;
rl_ss=0;
dh_ss=0;
dl_ss=0;

lamsh_ss=1/(1/betta-1)*(mh+(1-eta)*log(csh1_ss-cbar)+eta*log(csh2_ss)-theta/2*(nsh_ss)^2+lambh_ss*( wh_ss*nsh_ss-P1_ss*csh1_ss-P2_ss*csh2_ss));
lamsl_ss=1/(1/betta-1)*(ml+(1-eta)*log(csl1_ss-cbar)+eta*log(csl2_ss)-theta/2*(nsl_ss)^2+lambl_ss*( wl_ss*nsl_ss-P1_ss*csl1_ss-P2_ss*csl2_ss));

lamrh_ss=1/(1/betta-1)*(mh+(1-eta)*log(crh1_ss-cbar)+eta*log(crh2_ss)-theta/2*(nrh_ss)^2+lambh_ss*( wh_ss*nrh_ss-P1_ss*crh1_ss-P2_ss*crh2_ss ));
lamrl_ss=1/(1/betta-1)*(ml+(1-eta)*log(crl1_ss-cbar)+eta*log(crl2_ss)-theta/2*(nrl_ss)^2+lambl_ss*( wl_ss*nrl_ss-P1_ss*crl1_ss-P2_ss*crl2_ss ));

lamih_ss=1/((1/betta-(1-pirh-pidh)))*(mh+(1-eta)*log(cih1_ss-cbar)+eta*log(cih2_ss)-theta/2*(nih_ss)^2+lambh_ss*( wh_ss*nih_ss-P1_ss*cih1_ss-P2_ss*cih2_ss )+lamrh_ss*pirh);
lamil_ss=1/((1/betta-(1-pirl-pidl)))*(ml+(1-eta)*log(cil1_ss-cbar)+eta*log(cil2_ss)-theta/2*(nil_ss)^2+lambl_ss*( wl_ss*nil_ss-P1_ss*cil1_ss-P2_ss*cil2_ss )+lamrl_ss*pirl);
 
lamtauh_ss=lamih_ss-lamsh_ss;
lamtaul_ss=lamil_ss-lamsl_ss;

Uh_ss=1/(1-betta)*sh_ss*(mh+(1-eta)*log(csh1_ss-cbar)+eta*log(csh2_ss)-theta/2*(nsh_ss)^2);
Ul_ss=1/(1-betta)*sl_ss*(ml+(1-eta)*log(csl1_ss-cbar)+eta*log(csl2_ss)-theta/2*(nsl_ss)^2);
U_ss=Uh_ss+Ul_ss;

VoL_h_mill_pc=Uh_ss/lambh_ss/10000000/sh_ss;
VoL_l_mill_pc=Ul_ss/lambl_ss/10000000/sl_ss;
 
Cl_ss=(P1_ss*sl_ss*csl1_ss+P2_ss*sl_ss*csl2_ss); 
Ch_ss=(P1_ss*sh_ss*csh1_ss+P2_ss*sh_ss*csh2_ss); 

Nl_ss=sl_ss*nsl_ss; 
Nh_ss=sh_ss*nsh_ss; 

NINCl_ss=sl_ss*wl_ss*nsl_ss;  
NINCh_ss=sh_ss*wh_ss*nsh_ss; 

profits_ss=(P1_ss*Al*( sl_ss*nsl_ss )-wl_ss*( sl_ss*nsl_ss))+(P2_ss*Ah*( sh_ss*nsh_ss)-wh_ss*( sh_ss*nsh_ss));
int_ss=1+rstar;

Gamh_ss=0;
Gaml_ss=0;

inch_ss=wh_ss*(sh_ss*nsh_ss)+rstar*P2_ss*bhstar_ss+profits_ss+Gamh_ss;
incl_ss=wl_ss*(sl_ss*nsl_ss)+Gaml_ss;

incshrh_ss=inch_ss/(inch_ss+incl_ss);

Score1=(1-incshrh_ss)*((1-sh_ss)+2*sh_ss);
Score2=incshrh_ss*(sh_ss);
Gini_ss=1-Score1-Score2; 

//calibrate the pi's in Transmisson function
go_calib_pis;  
 
//final numbers for pi's will be imposed below using homotopy
pi1_final=pi1;
pi2_final=pi2;
pil3_final=pil3;
pih3_final=pih3;
pi4_final=pi4;

//put scaled value of pi4 into Dynare M_. structure
//if you put the final number from the get go, no solution will be found
//use smaller numbers first then compute a solution. Then use the solution 
//as an initial guess for slightly larger values. Proceed using this homotopy 
//until final values are imposed (all below, after model block). 
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pil3',M_.param_names,'exact'))=pil3_final;
M_.params(strmatch('pih3',M_.param_names,'exact'))=pih3_final;
M_.params(strmatch('pi4',M_.param_names,'exact'))=pi4_final/scale_pi4_homotopy;



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

//Equilibrium equations
(sh(-1)*csh1+ih(-1)*cih1+rh(-1)*crh1 ) + ( sl(-1)*csl1+il(-1)*cil1+rl(-1)*crl1 )=Al*( sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl );

P2*diffbhstar+P1*(sh(-1)*csh1+ih(-1)*cih1+rh(-1)*crh1)+P2*(sh(-1)*csh2+ih(-1)*cih2+rh(-1)*crh2)=inch;
diffbhstar=bhstar-bhstar(-1);
diffbhstar=bhstareps;
residEuler=dlambh(+1)*dP2(+1)*(1+rstar)*betta-1;
//note that the path of assets is nonstationary in the model. We solve for the path of assets as follows. bhstareps is a
//vector of exogenous variables. In the script comp_bhstar.m (called below) we are searching for the realizations of bhstareps
//over time such that the first order condition for bonds holds, i.e. the Euler equation residuals
//are zero. Check out the file comp_bhstar.m for the details.

inch=wh*(sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh)+rstar*P2*bhstar(-1)+profits+Gamh;
incl=wl*(sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl)+Gaml;

incshrh=inch/(inch+incl);

profits=(P1*Al*( sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl )-wl*( sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl ))
+(P2*Ah*( sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh )-wh*( sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh ));

P1*(sl(-1)*csl1+il(-1)*cil1+rl(-1)*crl1)+P2*(sl(-1)*csl2+il(-1)*cil2+rl(-1)*crl2)=incl;
P2=P2_ss_numeraire_normalization;

wl=P1*Al;
wh=P2*Ah;

C=P1*( sh(-1)*csh1+ih(-1)*cih1+rh(-1)*crh1 + sl(-1)*csl1+il(-1)*cil1+rl(-1)*crl1   )
+P2*(  sh(-1)*csh2+ih(-1)*cih2+rh(-1)*crh2 + sl(-1)*csl2+il(-1)*cil2+rl(-1)*crl2    );

N=(sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl ) + ( sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh );

taushift=(1-taushk);

tauh=taushift*(pi1*sh(-1)*csh1*( ih(-1)*cih1 + il(-1)*cil1 ) + pi2*sh(-1)*csh2*( ih(-1)*cih2 + il(-1)*cil2 )+pih3*sh(-1)*nsh*ih(-1)*nih + pi4*sh(-1)*( ih(-1)+il(-1) ));
taul=taushift*(pi1*sl(-1)*csl1*( ih(-1)*cih1 + il(-1)*cil1 ) + pi2*sl(-1)*csl2*( ih(-1)*cih2 + il(-1)*cil2 )+pil3*sl(-1)*nsl*il(-1)*nil + pi4*sl(-1)*( ih(-1)+il(-1) ));

sh=sh(-1)-tauh;
sl=sl(-1)-taul;

ih=(1-pirh-pidh)*ih(-1)+tauh;
il=(1-pirl-pidl)*il(-1)+taul;

rh=rh(-1)+pirh*ih(-1);
rl=rl(-1)+pirl*il(-1);

dh=dh(-1)+pidh*ih(-1);
dl=dl(-1)+pidl*il(-1);

consshift=(1-contain); //Containment variable

mh+(consshift(+1))*(1-eta)*log(csh1*dcsh1(+1)-cbar)+(consshift(+1))*eta*log(csh2*dcsh2(+1))-theta/2*(nsh*dnsh(+1))^2
+lamtauh*dlamtauh(+1)*taushift(+1)*(pi1*csh1*dcsh1(+1)*( ih*cih1*dcih1(+1) + il*cil1*dcil1(+1) ) + pi2*csh2*dcsh2(+1)*( ih*cih2*dcih2(+1) + il*cil2*dcil2(+1) )+pih3*nsh*dnsh(+1)*ih*nih*dnih(+1) + pi4*( ih+il))
+lambh*dlambh(+1)*( wh*dwh(+1)*nsh*dnsh(+1)-P1*dP1(+1)*csh1*dcsh1(+1)-P2*dP2(+1)*csh2*dcsh2(+1) )
-lamsh/betta+lamsh*dlamsh(+1);
%level consshift(+1)) abive is ok since it is a temporary shock

ml+(consshift(+1))*(1-eta)*log(csl1*dcsl1(+1)-cbar)+(consshift(+1))*eta*log(csl2*dcsl2(+1))-theta/2*(nsl*dnsl(+1))^2
+lamtaul*dlamtaul(+1)*taushift(+1)*(pi1*csl1*dcsl1(+1)*( ih*cih1*dcih1(+1) + il*cil1*dcil1(+1) ) + pi2*csl2*dcsl2(+1)*( ih*cih2*dcih2(+1) + il*cil2*dcil2(+1) )+pil3*nsl*dnsl(+1)*il*nil*dnil(+1) + pi4*( ih+il))
+lambl*dlambl(+1)*( wl*dwl(+1)*nsl*dnsl(+1)-P1*dP1(+1)*csl1*dcsl1(+1)-P2*dP2(+1)*csl2*dcsl2(+1) )
-lamsl/betta+lamsl*dlamsl(+1);

mh+(consshift(+1))*(1-eta)*log(cih1*dcih1(+1)-cbar)+(consshift(+1))*eta*log(cih2*dcih2(+1))-theta/2*(nih*dnih(+1))^2
+lambh*dlambh(+1)*( wh*dwh(+1)*nih*dnih(+1)-P1*dP1(+1)*cih1*dcih1(+1)-P2*dP2(+1)*cih2*dcih2(+1) )
-lamih/betta+lamih*dlamih(+1)*(1-pirh-pidh)+lamrh*dlamrh(+1)*pirh;

ml+(consshift(+1))*(1-eta)*log(cil1*dcil1(+1)-cbar)+(consshift(+1))*eta*log(cil2*dcil2(+1))-theta/2*(nil*dnil(+1))^2
+lambl*dlambl(+1)*( wl*dwl(+1)*nil*dnil(+1)-P1*dP1(+1)*cil1*dcil1(+1)-P2*dP2(+1)*cil2*dcil2(+1) )
-lamil/betta+lamil*dlamil(+1)*(1-pirl-pidl)+lamrl*dlamrl(+1)*pirl;

mh+(consshift(+1))*(1-eta)*log(crh1*dcrh1(+1)-cbar)+(consshift(+1))*eta*log(crh2*dcrh2(+1))-theta/2*(nrh*dnrh(+1))^2
+lambh*dlambh(+1)*( wh*dwh(+1)*nrh*dnrh(+1)-P1*dP1(+1)*crh1*dcrh1(+1)-P2*dP2(+1)*crh2*dcrh2(+1) )
-lamrh/betta+lamrh*dlamrh(+1);

ml+(consshift(+1))*(1-eta)*log(crl1*dcrl1(+1)-cbar)+(consshift(+1))*eta*log(crl2*dcrl2(+1))-theta/2*(nrl*dnrl(+1))^2
+lambl*dlambl(+1)*( wl*dwl(+1)*nrl*dnrl(+1)-P1*dP1(+1)*crl1*dcrl1(+1)-P2*dP2(+1)*crl2*dcrl2(+1) )
-lamrl/betta+lamrl*dlamrl(+1);

lamih=lamtauh+lamsh;
lamil=lamtaul+lamsl;

(consshift)*(1-eta)/(csh1-cbar)=P1*lambh-lamtauh*taushift*pi1*(ih(-1)*cih1+il(-1)*cil1);
(consshift)*(1-eta)/(csl1-cbar)=P1*lambl-lamtaul*taushift*pi1*(ih(-1)*cih1+il(-1)*cil1);

(consshift)*eta/csh2=P2*lambh-lamtauh*taushift*pi2*(ih(-1)*cih2+il(-1)*cil2);
(consshift)*eta/csl2=P2*lambl-lamtaul*taushift*pi2*(ih(-1)*cih2+il(-1)*cil2);

theta*nsh-wh*lambh-lamtauh*taushift*pih3*ih(-1)*nih;
do_sticky_wages*(wl-STEADY_STATE(wl))=(1-do_sticky_wages)*(theta*nsl-wl*lambl-lamtaul*taushift*pil3*il(-1)*nil);

(consshift)*(1-eta)/(cih1-cbar)=P1*lambh;
(consshift)*(1-eta)/(cil1-cbar)=P1*lambl;

(consshift)*eta/cih2=P2*lambh;
(consshift)*eta/cil2=P2*lambl;

theta*nih=wh*lambh;
do_sticky_wages*(nil-nsl)=(1-do_sticky_wages)*(theta*nil-wl*lambl);

(consshift)*(1-eta)/(crh1-cbar)=P1*lambh;
(consshift)*(1-eta)/(crl1-cbar)=P1*lambl;

(consshift)*eta/crh2=P2*lambh;
(consshift)*eta/crl2=P2*lambl;

theta*nrh=wh*lambh;
do_sticky_wages*(nrl-nsl)=(1-do_sticky_wages)*(theta*nrl-wl*lambl);

Gaml=transhk;
bgov-Gamh*(sh(-1)+ih(-1)+rh(-1))=Gaml*(sl(-1)+il(-1)+rl(-1))+(1+rstar)*bgov(-1);
gov_switch*(Gamh-0)=(1-gov_switch)*((sh(-1)+ih(-1)+rh(-1))*Gamh+rstar*bgov(-1)); 

Cl=(P1*(sl(-1)*csl1+il(-1)*cil1+rl(-1)*crl1)+P2*(sl(-1)*csl2+il(-1)*cil2+rl(-1)*crl2)); 
Ch=(P1*(sh(-1)*csh1+ih(-1)*cih1+rh(-1)*crh1)+P2*(sh(-1)*csh2+ih(-1)*cih2+rh(-1)*crh2));  

Nl=(sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl);
Nh=(sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh);

NINCl=wl*(sl(-1)*nsl+il(-1)*nil+rl(-1)*nrl);
NINCh=wh*(sh(-1)*nsh+ih(-1)*nih+rh(-1)*nrh);

Ul=sl(-1)*(ml+(consshift)*(1-eta)*log(csl1-cbar)+(consshift)*eta*log(csl2)-theta/2*(nsl)^2)
+il(-1)*(ml+(consshift)*(1-eta)*log(cil1-cbar)+(consshift)*eta*log(cil2)-theta/2*(nil)^2)
+rl(-1)*(ml+(consshift)*(1-eta)*log(crl1-cbar)+(consshift)*eta*log(crl2)-theta/2*(nrl)^2)
+betta*Ul*dUl(+1);

Uh=sh(-1)*(mh+(consshift)*(1-eta)*log(csh1-cbar)+(consshift)*eta*log(csh2)-theta/2*(nsh)^2)
+ih(-1)*(mh+(consshift)*(1-eta)*log(cih1-cbar)+(consshift)*eta*log(cih2)-theta/2*(nih)^2)
+rh(-1)*(mh+(consshift)*(1-eta)*log(crh1-cbar)+(consshift)*eta*log(crh2)-theta/2*(nrh)^2)
+betta*Uh*dUh(+1);

U=Uh+Ul;

%Standard SIR model equations
tauhSIRmod=taushift*(pi1*shSIRmod(-1)*STEADY_STATE(csh1)*( ihSIRmod(-1)*STEADY_STATE(cih1) + ilSIRmod(-1)*STEADY_STATE(cil1) ) + pi2*shSIRmod(-1)*STEADY_STATE(csh2)*( ihSIRmod(-1)*STEADY_STATE(cih2) + ilSIRmod(-1)*STEADY_STATE(cil2) )+pih3*shSIRmod(-1)*STEADY_STATE(nsh)*ihSIRmod(-1)*STEADY_STATE(nih) + pi4*shSIRmod(-1)*( ihSIRmod(-1)+ilSIRmod(-1) ));
taulSIRmod=taushift*(pi1*slSIRmod(-1)*STEADY_STATE(csl1)*( ihSIRmod(-1)*STEADY_STATE(cih1) + ilSIRmod(-1)*STEADY_STATE(cil1) ) + pi2*slSIRmod(-1)*STEADY_STATE(csl2)*( ihSIRmod(-1)*STEADY_STATE(cih2) + ilSIRmod(-1)*STEADY_STATE(cil2) )+pil3*slSIRmod(-1)*STEADY_STATE(nsl)*ilSIRmod(-1)*STEADY_STATE(nil) + pi4*slSIRmod(-1)*( ihSIRmod(-1)+ilSIRmod(-1) ));

shSIRmod=shSIRmod(-1)-tauhSIRmod;
slSIRmod=slSIRmod(-1)-taulSIRmod;

ihSIRmod=(1-pirh-pidh)*ihSIRmod(-1)+tauhSIRmod;
ilSIRmod=(1-pirl-pidl)*ilSIRmod(-1)+taulSIRmod;

rhSIRmod=rhSIRmod(-1)+pirh*ihSIRmod(-1);
rlSIRmod=rlSIRmod(-1)+pirl*ilSIRmod(-1);

dhSIRmod=dhSIRmod(-1)+pidh*ihSIRmod(-1);
dlSIRmod=dlSIRmod(-1)+pidl*ilSIRmod(-1);
 
//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dcsh1	=	csh1	/	csh1	(-1);
dcsh2	=	csh2	/	csh2	(-1);
dcih1	=	cih1	/	cih1	(-1);
dcih2	=	cih2	/	cih2	(-1);
dcil1	=	cil1	/	cil1	(-1);
dcil2	=	cil2	/	cil2	(-1);
dnsh	=	nsh	/	nsh	(-1);
dnih	=	nih	/	nih	(-1);
dwh	=	wh	/	wh	(-1);
dlamsh	=	lamsh	/	lamsh	(-1);
dlamtauh	=	lamtauh	/	lamtauh	(-1);
dlambh	=	lambh	/	lambh	(-1);
dP1	=	P1	/	P1	(-1);
dP2	=	P2	/	P2	(-1);
dcsl1	=	csl1	/	csl1	(-1);
dcsl2	=	csl2	/	csl2	(-1);
dnsl	=	nsl	/	nsl	(-1);
dnil	=	nil	/	nil	(-1);
dwl	=	wl	/	wl	(-1);
dlamsl	=	lamsl	/	lamsl	(-1);
dlamtaul	=	lamtaul	/	lamtaul	(-1);
dlambl	=	lambl	/	lambl	(-1);
dlamih	=	lamih	/	lamih	(-1);
dlamrh=lamrh/lamrh(-1);
dlamil =	lamil	/	lamil	(-1);
dlamrl	=	lamrl	/	lamrl	(-1);
dcrh1	=	crh1	/	crh1	(-1);
dcrh2	=	crh2	/	crh2	(-1);
dnrh	=	nrh	/	nrh	(-1);
dcrl1	=	crl1	/	crl1	(-1);
dcrl2	=	crl2	/	crl2	(-1);
dnrl	=	nrl	/	nrl	(-1);
dUh=Uh/Uh(-1);
dUl=Ul/Ul(-1);

@# include "commonVarEq.mod"  // added by epi-mmb team
end; 
 
 
initval;
P1=P1_ss;
P2=P2_ss;
C=C_ss;
N=N_ss;
wl=wl_ss;
wh=wh_ss;
taul=taul_ss;
tauh=tauh_ss;
sh=sh_ss;
ih=ih_ss;
rh=rh_ss;
dh=dh_ss;
sl=sl_ss;
il=il_ss;
rl=rl_ss;
dl=dl_ss;
nsh=nsh_ss;
nih=nih_ss;
nrh=nrh_ss;
nsl=nsl_ss;
nil=nil_ss;
nrl=nrl_ss;
csh1=csh1_ss;
cih1=cih1_ss;
crh1=crh1_ss;
csh2=csh2_ss;
cih2=cih2_ss;
crh2=crh2_ss;
csl1=csl1_ss;
cil1=cil1_ss;
crl1=crl1_ss;
csl2=csl2_ss;
cil2=cil2_ss;
crl2=crl2_ss;
lambh=lambh_ss;
lambl=lambl_ss;
lamtaul=lamtaul_ss;
lamtauh=lamtauh_ss;
lamil=lamil_ss;
lamih=lamih_ss;
lamsl=lamsl_ss; 
lamsh=lamsh_ss;
lamrl=lamrl_ss;
lamrh=lamrh_ss;
Cl=Cl_ss; 
Ch=Ch_ss; 
Nl=Nl_ss; 
Nh=Nh_ss;
dcsh1	=	1;
dcsh2	=	1;
dcih1	=	1;
dcih2	=	1;
dcil1	=	1;
dcil2	=	1;
dnsh	=	1;
dnih	=	1;
dwh	=	1;
dlamsh	=	1;
dlamtauh	=	1;
dlambh	=	1;
dP1	=	1;
dP2	=	1;
dcsl1	=	1;
dcsl2	=	1;
dnsl	=	1;
dnil	=	1;
dwl	=	1;
dlamsl	=	1;
dlamtaul	=	1;
dlambl	=	1;
dlamih	=	1;
dlamrh=1;
dlamil	=	1;
dlamrl	=	1;
dcrh1	=	1;
dcrh2	=	1;
dnrh	=	1;
dcrl1	=	1;
dcrl2	=	1;
dnrl	=	1;
bhstar=bhstar_ss;
NINCl=NINCl_ss; 
NINCh=NINCh_ss;
Uh=Uh_ss;
Ul=Ul_ss;
U=U_ss;
dUh=1;
dUl=1;
Gamh=0;
Gaml=0;
profits=profits_ss;
diffbhstar=0;
residEuler=0;
inch=inch_ss;
incl=incl_ss;
incshrh=incshrh_ss;
consshift=1;
taushift=1;
bgov=0;
taulSIRmod=0; 
tauhSIRmod=0; 
slSIRmod=sl_ss; 
shSIRmod=sh_ss; 
ilSIRmod=0; 
ihSIRmod=0; 
rlSIRmod=0; 
rhSIRmod=0; 
dlSIRmod=0; 
dhSIRmod=0;
@# include "commonVarSS.mod" //added by epi-mmb team
end;

//check residuals of dynamic equations given initial steady state
resid;
steady;
resid;


//set initial seed of infection
M_.endo_histval=oo_.steady_state;

M_.endo_histval(strmatch('ih',M_.endo_names,'exact'))=ihSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('il',M_.endo_names,'exact'))=ilSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('sh',M_.endo_names,'exact'))=shSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('sl',M_.endo_names,'exact'))=slSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('dh',M_.endo_names,'exact'))=dhSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('dl',M_.endo_names,'exact'))=dlSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('rh',M_.endo_names,'exact'))=rhSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('rl',M_.endo_names,'exact'))=rlSIR(wakeupWeek-SIRstartWeek);

M_.endo_histval(strmatch('ihSIRmod',M_.endo_names,'exact'))=ihSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('ilSIRmod',M_.endo_names,'exact'))=ilSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('shSIRmod',M_.endo_names,'exact'))=shSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('slSIRmod',M_.endo_names,'exact'))=slSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('dhSIRmod',M_.endo_names,'exact'))=dhSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('dlSIRmod',M_.endo_names,'exact'))=dlSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('rhSIRmod',M_.endo_names,'exact'))=rhSIR(wakeupWeek-SIRstartWeek);
M_.endo_histval(strmatch('rlSIRmod',M_.endo_names,'exact'))=rlSIR(wakeupWeek-SIRstartWeek);


//initialize exogenous variables
shocks;
var bhstareps;
periods 1:1;
values 0;

var  transhk;
periods 1:1;
values 0;

var  taushk;
periods 1:1;
values 0;

var  contain;
periods 1:1;
values 0;

var gov_switch;
periods 1:1;
values 0;
end;

//solve and simulate model
simul(stack_solve_algo=0,no_homotopy);

%get Chetty et al data
getdata;
 
//set exogenous variables
%transfers to low skilled
oo_.exo_simul(5,2)=15; 
oo_.exo_simul(6:18,2)=40;  
oo_.exo_simul(19:19+26-1,2)=linspace(35,0,26);
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);

//reduction in virus spread due to summer/reorga
oo_.exo_simul(6,3)=0.25;
oo_.exo_simul(7,3)=0.5;
oo_.exo_simul(8:18,3)=0.7;
oo_.exo_simul(19,3)=0.5;
oo_.exo_simul(20,3)=0.25;
oo_.exo_simul(21,3)=0.0;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);

//containment
oo_.exo_simul(2,4)=0.2;
oo_.exo_simul(3:11,4)=0.3; 
oo_.exo_simul(12,4)=0.25;
oo_.exo_simul(13,4)=0.225;
oo_.exo_simul(14:18,4)=0.2;
oo_.exo_simul(19:23,4)=0.15;
oo_.exo_simul(24:28,4)=0.1;
oo_.exo_simul(29:33,4)=0.05;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
 
//fiscal rule deactived in first two years
oo_.exo_simul(2:104,5)=1;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
 

//Homotopy for pi4 parameter (otherwise no solution if you set it to final values right away)
disp('pi4 homotopy');
pi4_final_steps_vec=linspace(pi4_final/scale_pi4_homotopy,pi4_final,homotopy_steps);
for pi4_final_step=pi4_final_steps_vec
    M_.params(strmatch('pi4',M_.param_names,'exact'))=pi4_final_step;
    //simulate using previous solution as initial guess
    [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
end    

 
//Simulate model with final (desired, calibrated) values for pi's
M_.params(strmatch('pi1',M_.param_names,'exact'))=pi1_final;
M_.params(strmatch('pi2',M_.param_names,'exact'))=pi2_final;
M_.params(strmatch('pih3',M_.param_names,'exact'))=pih3_final;
M_.params(strmatch('pil3',M_.param_names,'exact'))=pil3_final;
M_.params(strmatch('pi4',M_.param_names,'exact'))=pi4_final;
//simulate using previous solution as initial guess
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
//convert simulation results to model variable names defined above 
dyn2vec(M_, oo_, options_);   


 
///////////////////////////////////////////////////////////////////////////////////
//compute assetholdings bhstar such that asset FOC holds and welfare is maximized//
comp_bhstar;
//////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////
//Report Stats in Command Window
/////////////////////////////////
//Annual agg. consumption per capita
Cann_pc=C_ss*52

//Total Hours per capita
Ntot_pc=N_ss

//target 18%
sh_ss

//target 38%
wageincomeshare_high_skilled=(sh_ss*wh_ss*nsh_ss)/(sh_ss*wh_ss*nsh_ss+sl_ss*wl_ss*nsl_ss)

//skill premium
skill_premium=wageincomeshare_high_skilled/sh_ss

//wage premium
wh_wl_ratio=wh_ss/wl_ss

//P1 and P2
P1ss=P1_ss
P2ss=P2_ss

cbar_ss=cbar
theta_ss=theta
Ah_ss=Ah
Al_ss=Al

//value of life
VoL_h_mill_pc=Uh_ss/lambh_ss/10000000/sh_ss
VoL_l_mill_pc=Ul_ss/lambl_ss/10000000/sl_ss
VoL_wavg_mill_pc=sh_ss*VoL_h_mill_pc+sl_ss*VoL_l_mill_pc

//Some checking
den=pi1*(sl_ss*csl1_ss+sh_ss*csh1_ss)*(sh_ss*cih1_ss+sl_ss*cil1_ss)+pi2*(sl_ss*csl2_ss+sh_ss*csh2_ss)*(sh_ss*cih2_ss+sl_ss*cil2_ss)+pil3*sl_ss*nsl_ss*nil_ss*sl_ss+pih3*sh_ss*nsh_ss*nih_ss*sh_ss+pi4;
pi1_share=pi1*(sl_ss*csl1_ss+sh_ss*csh1_ss)*(sh_ss*cih1_ss+sl_ss*cil1_ss)/den
pi2_share=pi2*(sl_ss*csl2_ss+sh_ss*csh2_ss)*(sh_ss*cih2_ss+sl_ss*cil2_ss)/den
pil3_share=pil3*sl_ss*nsl_ss*nil_ss*sl_ss/den
pih3_share=pih3*sh_ss*nsh_ss*nih_ss*sh_ss/den
pi4_share=pi4/den
maxabs_residEuler=max(abs(residEuler))


//merge SIR and SIR-macro sims
oo_.endo_simul=[repmat(oo_.steady_state,1,wakeupWeek-SIRstartWeek) oo_.endo_simul(:,2:end)];
oo_.endo_simul(strmatch('ih',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=ihSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('il',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=ilSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('sh',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=shSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('sl',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=slSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('rh',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=rhSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('rl',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=rlSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('dh',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=dhSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('dl',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=dlSIR(1:wakeupWeek-SIRstartWeek);


oo_.endo_simul(strmatch('ihSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=ihSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('ilSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=ilSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('shSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=shSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('slSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=slSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('rhSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=rhSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('rlSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=rlSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('dhSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=dhSIR(1:wakeupWeek-SIRstartWeek);
oo_.endo_simul(strmatch('dlSIRmod',M_.endo_names,'exact'),SIRstartWeek:wakeupWeek-1)=dlSIR(1:wakeupWeek-SIRstartWeek);
dyn2vec(M_, oo_, options_); 

//Gini and deaths
incshrh_ss=incshrh(1)
Score1=(1-incshrh(1))*((1-sh_ss)+2*sh_ss)
Score2=incshrh(1)*(sh_ss)
Gini=1-Score1-Score2
dl_end=dl(end-1)
dh_end=dh(end-1)
D_end=dl(end-1)+dh(end-1)


//Plot results
% plot_results;

//Save results 
% save all_results_baseline 
results.oo_ = oo_ ;
results.M_ = M_;
save simulated_results_base;
@# include "saveResults.mod"   //added by epi-mmb team
 


 
