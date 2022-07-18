// This file solves and simulates the model developed in 
// Giagheddu and Papetti (2020), 
// 'The macroeconomics of age-varying epidemics'

//This code was written and run with Matlab 2017a and DYNARE 4.6.1. The nonlinear
//model is solved and simulated with Dynare's nonlinear deterministic solver
////////////////////////////////////////////////////////////
//IMPORTANT: You must use DYNARE 4.6.1 to run this code.//// 
////////////////////////////////////////////////////////////
// Andrea Papetti, papetti.andrea@gmail.com
///////////////////////////////////////////////////////////
// This code follows closely the publicly available code for
// Eichenbaum, Rebelo and Trabandt (2020),'The Macroeconomics of Testing and Quarantining'
// see https://sites.google.com/site/mathiastrabandt/home/research


// endogenous variables
var lambry lambiy lamsy lamtauy
cry cro ciy cio csy cso
nry niy nsy
Gamma
ury uro uiy uio usy uso
Ury Uro Uiy Uio Usy Uso
Ty Sy Iy Ry Dy Ny
To So Io Ro Do No
tauy tauo N 
dUry dUro dUiy dUio dUsy dUso
piy_d_endo pio_d_endo
aggC aggH ;

@#include "commonVar.mod"    // added by epi-mmb team

parameters 
do_SIR_canonical betta gamma
pio_d piy_d pio_r piy_r fy fo
Zyy Zyo Zoy Zoo
R0 eta eps phii deltav deltac kappa
n_target inc_target
A theta Pbar
pisy1 pisy2 pisy3 pisy4 piso1 piso2 piso3 other_share;

varexo muc;

//Initialize parameters
do_SIR_canonical=0; // if=1, SIRage_canonical (no macroeconomic interactions); if =0 SIRage_macro

stepps=100; //homotopy steps for pis

betta=0.96^(1/52);  %Weekly household discount factor
gamma=7/18;         % weekly ``removal'' rate (the rate at which an infected individual becomes either recovered or dead)

pio_d=0.1273*gamma;     %Weekly probability of dying for old
piy_d=0.0045*gamma;     %Weekly probability of dying for young
pio_r=gamma-pio_d;      %Weekly probability of recovering for old
piy_r=gamma-piy_d;      %Weekly probability of recovering for young

fy=.825; % share of young [0,70) in the economy for Italy, see demo_covidMortality.xlsx 
fo=1-fy; %share of old (70+) in the economy, 

% Contact matrix (see http://sherrytowers.com/2012/12/11/sir-model-with-age-classes/#code and references therein)
% From R code: load_contact_matrix.R (YOUNG: [0,70); OLD: 70+)
Zyy=19.123321*7;
Zyo=1.337484*7;
Zoy=fy*Zyo*7/(1-fy);
Zoo=1.391304*7;
% zyy=19.123321;
% zyo=1.337484;
% zoo=1.391304;
% Z=[zyy*7 zyo*7;
%    fy*zyo*7/(1-fy) zoo*7];

R0=1.590695809; % We change it to have Ry(end)+Ro(end)+Dy(end)+Do(end)=0.6 in "SIRage" model
% similar to see Eichenbaum et al. (2020 April 3, p 16)

% Auxilary MM-matrix to target R0 (see http://sherrytowers.com/2012/12/11/sir-model-with-age-classes/#code and references therein)
% M=[Z(1,1)*fy/fy Z(1,2)*fy/fo;
%     Z(2,1)*fo/fy Z(2,2)*fo/fo];
% eta=R0*gamma/(max(eig(M))); % transmission rate on contact (equal across age classes)
eta=R0*gamma/1.371075936665521e+02; % transmission rate on contact (equal across age classes)

load inf_ini
eps=helper; 
//eps=0.001; % share of initial aggregate infected (as in Eichenbaum et al. (2020))

phii=0.8;           %Productivity of infected (young) people
 
deltav=0/52;        %Weekly probability of discovering a vaccine
deltac=0/52;        %Weekly probability of discovering a treatment
kappa=0;            %Slope of pid-function in endog. pid scenario 

%Calibration targets for hours and income
n_target=28;         %Weekly hours
inc_target=58000/52; %weekly income

theta=1/n_target^2;     %calculate disutility of labor parameter theta
                        %so that pre-infection steady state labor is equal to
                        %desired target (using n=(1/theta)^(1/2) pre-infection 
                        %steady state equation)
A=inc_target/n_target;  %calculate parameter A such that income is equal to
                        %desired target, c=inc_target=A*n_target
                        
% muc=0;

//pre-infection steady state
nry_ss=(1/theta)^(1/2);           %labor young recovered (same as post-infection steady state)
cry_ss=A*nry_ss;                    %consumption young recovered
alp=.8;
Pbar=alp*cry_ss;
cro_ss=Pbar;                      % condumption old recovered
ury_ss=log(cry_ss)-theta/2*nry_ss^2;  %utility young recovered
Ury_ss=ury_ss/(1-betta);          %PV utility young recovered
uro_ss=log(cro_ss);                %utility old recovered
Uro_ss=uro_ss/(1-betta);          %PV utility old recovered
% UrssConsUnits_y=Ury_ss*cry_ss;        %PV utility in cons. units (Urss*Marg.Util.Cons); value of life
UrssConsUnits_o=Uro_ss*cro_ss;
niy_ss=(1/theta)^(1/2);           %labor young infected
ciy_ss=phii*A*niy_ss;               %consumption infected
cio_ss=Pbar; 
uiy_ss=log(ciy_ss)-theta/2*niy_ss^2;  %utility young infected
Uiy_ss=1/(1-(1-deltac)*betta*(1-piy_r-piy_d))*(uiy_ss...
+(1-deltac)*betta*piy_r*Ury_ss+deltac*betta*Ury_ss);  %PV utility infected
uio_ss=log(cio_ss);  %utility old infected
Uio_ss=1/(1-(1-deltac)*betta*(1-pio_r-pio_d))*(uio_ss...
+(1-deltac)*betta*pio_r*Uro_ss+deltac*betta*Uro_ss);  %PV utility infected
Usy_ss=Ury_ss; %PV utility of susceptibles same as recovered in steady state
Uso_ss=Uro_ss; %PV utility of susceptibles same as recovered in steady state

% complete the set of steady state variables (following baseline dyn. recursion)
muc_ss=0;
lambry_ss=theta*nry_ss/A;
% cry_ss=((1+muc).*lambry_ss).^(-1);
Gamma_ss=(1+muc_ss).*cry_ss-A.*nry_ss;
% cro_ss=(Pbar+Gamma_ss)./(1+muc_ss)
aggC_ss=(cry_ss*fy+cro_ss*fo);
aggH_ss=nry_ss*fy;
lambiy_ss=(theta*niy_ss)./(phii*A);
lamtauy_ss=(1-deltav)*betta*(Uiy_ss-Usy_ss);
lamsy_ss=(cry_ss^(-1))/(1+muc_ss);


% if do_SIR_canonical==1
% pisy1_final=0;
% pisy2_final=0;
% pisy3_final=0;
% pisy4_final=1;
% piso1_final=0;
% piso2_final=0;
% piso3_final=1;
% else

// Set the weight on contacts other than consumption/work to calibrate pis
// other_share=1 ==> SIR_age_canonical with no consumption/work impact on infection
// corresponding to case: pisy1=pisy2=pisy3=piso1=piso2=0; pisy4=piso3=1;
// But still need to go through the go_calibrate_pis due to the numerical
// procedure to solve the model (otherwise not even the canonical SIR transition is solved)
other_share=1*do_SIR_canonical +2/3*(1-do_SIR_canonical);

//calibrate the pi's in T-function
go_calibrate_pis;

//final numbers for pis will be imposed below using homotopy
pisy1_final=pisy1;
pisy2_final=pisy2;
pisy3_final=pisy3;
pisy4_final=pisy4;
piso1_final=piso1;
piso2_final=piso2;
piso3_final=piso3;
% end

//put scaled down values of pis into Dynare M_. structure
//if you put the final numbers from the get go, no solution will be found
//use smaller numbers first then compute a solution. Then use the solution 
//as an initial guess for slightly larger values. Proceed using this homotopy 
//until final values are imposed (all below, after model block). 
M_.params(strmatch('pisy1',M_.param_names,'exact'))=pisy1_final/1000;
M_.params(strmatch('pisy2',M_.param_names,'exact'))=pisy2_final/1000;
M_.params(strmatch('pisy3',M_.param_names,'exact'))=pisy3_final/1000;
M_.params(strmatch('pisy4',M_.param_names,'exact'))=pisy4_final/10000;
M_.params(strmatch('piso1',M_.param_names,'exact'))=piso1_final/1000;
M_.params(strmatch('piso2',M_.param_names,'exact'))=piso2_final/1000;
M_.params(strmatch('piso3',M_.param_names,'exact'))=piso3_final/10000;

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

/////////////////////////////////////////////////////////////////////////// 
// economic equilibrium equations for diffeerent types (S, I, R), age (y,o)
/////////////////////////////////////////////////////////////////////////// 
% RECOVERED
// (1) FOC labor, young recovered
lambry=(theta*nry)/A;

// (2) FOC consumption, young recovered
cry=((1+muc)*lambry)^(-1);

// (3) Budget, young recovered 
Gamma=(1+muc)*cry-A*nry;

// (4) Consumption, old recoverd
cro=(Pbar+Gamma)/(1+muc);

// (5) Flow utility, young recovered
ury=log(cry)-theta/2*nry^2;

// (6) Flow utility, old recovered
uro=log(cro);

// (7) Utility, young recovered
Ury=ury+betta*Ury*dUry(+1);

// (8) Utility, old recovered
Uro=uro+betta*Uro*dUro(+1);

% INFECTED
// (9) FOC labor, young infected
lambiy=(theta*niy)/(phii*A);

// (10) FOC consumptioon, young infected
ciy=((1+muc)*lambiy)^(-1);

// (11) Budget, young infected
(1+muc)*ciy-phii*A*niy-Gamma=0;

// (12) Consumption, old infected
cio=(Pbar+Gamma)/(1+muc);

// (13) Flow utility, young infected
uiy=log(ciy)-theta/2*niy^2;

// (14) Flow utility, old infected
uio=log(cio);

// (15) Utility, young infected
Uiy=uiy+(1-deltac)*betta*((1-piy_r-piy_d_endo)*Uiy*dUiy(+1)
+piy_r*Ury*dUry(+1))
+deltac*betta*Ury*dUry(+1);

// (16) Utility, old infected
Uio=uio+(1-deltac)*betta*((1-pio_r-pio_d_endo)*Uio*dUio(+1)
+pio_r*Uro*dUro(+1))
+deltac*betta*Uro*dUro(+1);

% SUSCEPTIBLE
// (17) FOC labor, young susceptible
-theta*nsy+A*lamsy+eta*Zyy*pisy3*(Iy(-1)/fy)*niy*lamtauy=0;

// (18) FOC consumption, young susceptible
lamsy=(csy^(-1)+eta*lamtauy*(pisy1*Zyy*(Iy(-1)/fy)*ciy
+pisy2*Zyo*(Io(-1)/fo)*cio))/(1+muc);

// (19) FOC taus, susceptible
lamtauy=(1-deltav)*betta*(Uiy*dUiy(+1)-Usy*dUsy(+1));

// (20) Budget, young susceptible
csy=1/(1+muc)*(A*nsy+Gamma);

// (21) Consumption, old susceptible
cso=(Pbar+Gamma)/(1+muc);

// (22) Flow utility, young susceptible
usy=log(csy)-theta/2*nsy^2;

// (23) Flow utility, old susceptible
uso=log(cso);

// (24) Utility, young suscepitble
Usy=usy+(1-deltav)*betta*(1-tauy)*Usy*dUsy(+1)
       +(1-deltav)*betta*tauy*Uiy*dUiy(+1)
       +deltav*betta*Ury*dUry(+1);

// (25) Utility, old suscepitble
Uso=uso+(1-deltav)*betta*(1-tauo)*Uso*dUso(+1)
                 +(1-deltav)*betta*tauo*Uio*dUio(+1)
                 +deltav*betta*Uro*dUro(+1);

/////////////////////////
// Population Dynamics
////////////////////////
// (26) New infections, young
Ty=eta*Sy(-1)*(
    pisy1*Zyy*(Iy(-1)/fy)*csy*ciy
    +pisy2*Zyo*(Io(-1)/fo)*csy*cio
    +pisy3*Zyy*(Iy(-1)/fy)*nsy*niy
    +pisy4*(Zyy*(Iy(-1)/fy)+Zyo*(Io(-1)/fo)));

// (27) New infections, old
To=eta*So(-1)*(
    piso1*Zoo*(Io(-1)/fo)*cso*cio
    +piso2*Zoy*(Iy(-1)/fy)*cso*ciy
    +piso3*(Zoy*(Iy(-1)/fy)+Zoo*(Io(-1)/fo))); 

// (28) Total susceptibles, young
Sy=Sy(-1)-Ty;

// (29) Total susceptibles, old
So=So(-1)-To;

// (30) Endogenous death probability, young ["medical preparedness model": kappa \neq 0]
piy_d_endo=piy_d+kappa*Iy(-1)^2;

// (31) Endogenous death probability, old ["medical preparedness model": kappa \neq 0]
pio_d_endo=pio_d+kappa*Io(-1)^2;

// (32) Total infected, young
Iy=Iy(-1)+Ty-(piy_d_endo+piy_r)*Iy(-1);

// (33) Total infected, old
Io=Io(-1)+To-(pio_d_endo+pio_r)*Io(-1);

// (34) Total recovered, young
Ry=Ry(-1)+piy_r*Iy(-1);

// (35) Total recovered, old
Ro=Ro(-1)+pio_r*Io(-1);

// (36) Total deaths, young
Dy=Dy(-1)+piy_d_endo*Iy(-1);

// (37) Total deaths, old
Do=Do(-1)+pio_d_endo*Io(-1);

// (38) Total population, young
Ny=Sy+Iy+Ry;

// (39) Total population, old
No=So+Io+Ro;

// (40) Total population, young + old
N=Ny+No;

// (41) Infection probability, young
tauy=Ty/Sy(-1);

// (42) Infection probability, old
tauo=To/So(-1);

% CLEARING
// (43) Aggregate consumption
aggC=Sy(-1)*csy+Iy(-1)*ciy+Ry(-1)*cry+So(-1)*cso+Io(-1)*cio+Ro(-1)*cro;

// (44) Aggregate hours
aggH=Sy(-1)*nsy+Iy(-1)*niy*phii+Ry(-1)*nry;

// (45) Budget, government
muc*aggC-Gamma*N(-1)=0;

@# include "commonVarEq.mod"  // added by epi-mmb team


//Auxilliary variables: gross growth rates of variables with non-zero 
//pre-infection steady states. These variables are needed to calculate
//numerically accurate simulations when the terminal steady state differs
//from the pre-infection steady state and if you do not know the terminal 
//steady state a priori since it depends on the epidemic dynamics.
dUry=Ury/Ury(-1);
dUro=Uro/Uro(-1); 
dUiy=Uiy/Uiy(-1);
dUio=Uio/Uio(-1);
dUsy=Usy/Usy(-1);
dUso=Uso/Uso(-1);
end;

//initial (pre-infection) steady state used in simulations 
initval;
lambry=lambry_ss;
lambiy=lambiy_ss;
lamsy=lamsy_ss; 
lamtauy=lamtauy_ss;
cry=cry_ss;
cro=cro_ss;
ciy=ciy_ss;
cio=cio_ss;
csy=cry_ss;
cso=cro_ss;
nry=nry_ss;
niy=niy_ss;
nsy=nry_ss;
Gamma=Gamma_ss;
ury=ury_ss;
uro=uro_ss;
uiy=uiy_ss;
uio=uio_ss;
usy=ury_ss;
uso=uro_ss;
Ury=Ury_ss;
Uro=Uro_ss;
Uiy=Uiy_ss;
Uio=Uio_ss;
Usy=Usy_ss;
Uso=Uso_ss;
Ty=0;
Sy=fy;
Iy=0;
Ry=0;
Dy=0;
Ny=fy;
To=0;
So=fo;
Io=0;
Ro=0; 
Do=0; 
No=fo;
tauy=0;
tauo=0;
N=1;
dUry=1;
dUro=1;
dUiy=1;
dUio=1;
dUsy=1;
dUso=1;
piy_d_endo=piy_d;
pio_d_endo=pio_d;
aggC=aggC_ss;
aggH=aggH_ss;
@# include "commonVarSS.mod" //added by epi-mmb team
muc=muc_ss;
end;

shocks;
    var     muc;
    periods 1:1;
    values  0.000;
end;

//check residuals of dynamic equations given initial steady state
resid;
steady;
resid;

//set unknown type state variables to zero. Initialize initial seed 
M_.endo_histval=oo_.steady_state;
M_.endo_histval(strmatch('Iy',M_.endo_names,'exact'))=eps*fy;
M_.endo_histval(strmatch('Io',M_.endo_names,'exact'))=eps*fo;
M_.endo_histval(strmatch('Sy',M_.endo_names,'exact'))=(1-eps)*fy;
M_.endo_histval(strmatch('So',M_.endo_names,'exact'))=(1-eps)*fo;

//solve and simulate model
options_.slowc=1;
options_.simul.maxit=20;
options_.dynatol.f=1e-8;
options_.dynatol.x=1e-8;
options_.nograph=1;
//simul(periods=500,stack_solve_algo=0);
simul(periods=500,stack_solve_algo=0,noprint);


//Homotopy for pi1 and pi2 parameters (otherwise no solution if you set them to final values right away)
pisy1_final_steps_vec=linspace(pisy1_final/1000,pisy1_final/1,stepps);
pisy2_final_steps_vec=linspace(pisy2_final/1000,pisy2_final/1,stepps);
pisy3_final_steps_vec=linspace(pisy3_final/1000,pisy3_final/1,stepps);
pisy4_final_steps_vec=linspace(pisy4_final/10000,pisy4_final/1,stepps);

piso1_final_steps_vec=linspace(piso1_final/1000,piso1_final/1,stepps);
piso2_final_steps_vec=linspace(piso2_final/1000,piso2_final/1,stepps);
piso3_final_steps_vec=linspace(piso3_final/10000,piso3_final/1,stepps);


for pi4_final_step=pisy4_final_steps_vec
    M_.params(strmatch('pisy4',M_.param_names,'exact'))=pi4_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
end
for pi3_final_step=pisy3_final_steps_vec
    M_.params(strmatch('pisy3',M_.param_names,'exact'))=pi3_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
end
for pisy2_final_step=pisy2_final_steps_vec
    M_.params(strmatch('pisy2',M_.param_names,'exact'))=pisy2_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
  
end
for pisy1_final_step=pisy1_final_steps_vec
    M_.params(strmatch('pisy1',M_.param_names,'exact'))=pisy1_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess

end
for piso3_final_step=piso3_final_steps_vec
    M_.params(strmatch('piso3',M_.param_names,'exact'))=piso3_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    
end
for piso2_final_step=piso2_final_steps_vec
    M_.params(strmatch('piso2',M_.param_names,'exact'))=piso2_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    
end
for piso1_final_step=piso1_final_steps_vec
    M_.params(strmatch('piso1',M_.param_names,'exact'))=piso1_final_step;
     [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess
    
end
dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 
oo_endo_simul=oo_.endo_simul;
save oo_endo_simul_old oo_endo_simul; //save simulation results for later use as initial guess

//Simulate model with final (desired, calibrated) values for pi1,pi2,pi3 
M_.params(strmatch('pisy1',M_.param_names,'exact'))=pisy1_final;
M_.params(strmatch('pisy2',M_.param_names,'exact'))=pisy2_final;
M_.params(strmatch('pisy3',M_.param_names,'exact'))=pisy3_final;
M_.params(strmatch('pisy4',M_.param_names,'exact'))=pisy4_final;
M_.params(strmatch('piso1',M_.param_names,'exact'))=piso1_final;
M_.params(strmatch('piso2',M_.param_names,'exact'))=piso2_final;
M_.params(strmatch('piso3',M_.param_names,'exact'))=piso3_final;

load oo_endo_simul_old;
oo_.endo_simul=oo_endo_simul;
 [oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_); //simulate using previous solution as initial guess

//dyn2vec(M_, oo_, options_);   //convert simulation results to model variable names defined above 
@# include "saveResults.mod"   //added by epi-mmb team

% //plot all endogenous variables
% ia=7; ib=7;
% numvar=size(oo_.endo_simul,1)-9;
% horz=500;
 
%First, plot all periods
% figure; 
% for ii=1:1:numvar
%     subplot(ia,ib,ii)
%     plot(oo_.endo_simul(ii,1:horz),'b-'); hold on;
%     %plot(oo_.steady_state(ii)*ones(numel(oo_.endo_simul(1,1:horz))),'k--'); hold off;
%     title(M_.endo_names(ii,:),'Interpreter','none');
% end
% orient landscape
% print -dpdf -fillpage sim_fig3  

//Plot results
// //plot_results;

//Save results
// //save base_SIRage_macro_GP
