// Dynare model file ADS Model 2020 
// Created by Veronica Acurio Vasconez 
// June 2020

// Variables 

var L C Z K Ym In Y R Rk Q Pm vp Pi Pistar W b i lambda_c Lambda phi eta nu Omega Gamma f fprime g1 g2 psi G V Q_Z V1 V2 ;

predetermined_variables K;

@#include "commonVar.mod"    // added by epi-mmb team

// Exogenous Variables 

varexo Itilde; 

// Parameters 

parameters alphav betav gammav epsilonv betta h theta lambda epsilon alfa delta omega_g kappa epsilon_p theta_p chi zetta phi_i phi_pi phi_y tau omega psi_bar ;


alphav=0.40000000;
betav=0.20000000;
gammav=0.10000000;
epsilonv=0.05000000;
betta=0.99000000;
h=0.71000000;
theta=0.97500000;
lambda=0.41126405;
epsilon=0.00135038;
alfa=0.33000000;
delta=0.02500000;
omega_g=0.18000000;
kappa=5.74000000;
epsilon_p=4.16700000;
theta_p=0.66000000;
chi=0.24000000;
zetta=0.50000000;
phi_i=0.81000000;
phi_pi=2.04000000;
phi_y=0.08000000;
tau=0.00100000;
omega=10.00000000;
psi_bar=0.00000000;


model;

//
//Necessary steady state values from baseline
//

        
//
//Epidemiological block (SIR)
//

//1. Subceptible dynamics
    //[mcp='S>0']
    //S = S(-1) - alphav*Itilde(-1)*S(-1);

//2. Infected dynamics
    //[mcp='Itilde<1']
    //Itilde = Itilde(-1) + alphav*Itilde(-1)*S(-1) - gammav*Itilde(-1);

//3. Recovered dynamics
    //[mcp='Rtilde<1']
    //Rtilde = Rtilde(-1) + gammav*Itilde(-1);

//4. Labor supply
L = 1 - Itilde;


//
//Households
//

//5. Marginal utility of consumption
lambda_c = 1/(C-h*C(-1)) - betta*h/(C(+1)-h*C);

//6. Euler equation
    //R(+1)*Lambda(+1) = 1;
(1+i)/Pi(+1)*Lambda(+1) = 1;

//7. Stochastic discount rate
    //Lambda(+1)  =   betta*lambda_c(+1)/lambda_c;
Lambda  =   betta*lambda_c/lambda_c(-1);


//
//Private banks
//

//8. Value of banks' capital
nu = Lambda(+1)*Gamma(+1)*(Rk(+1)-R(+1));

//9. Value of banks' net wealth
eta = Lambda(+1)*Gamma(+1)*R(+1);

//10. Auxiliar variable
Gamma = 1 - theta + theta*lambda*phi;

//11. Optimal leverage
phi = eta/(lambda - nu);

//12. Banks' net worth
Omega = theta*((Rk-R)*phi(-1) + R)*Omega(-1) + epsilon*Q*Z(-1);

//12b. Expected discounted terminal wealth
V = nu*Q*Z + eta*Omega;

//12c. Expected discounted value of banks' capital
Q_Z = Q*Z;

//12d. Expected discounted value of banks' capital
V1 = nu*Q*Z;

//12e. Expected discounted value of banks' net wealth
V2 = eta*Omega;


//
//Corona Bonds
//

//13. Aggregate capital
Q*Z = phi/(1-psi)*Omega;


//
//Intermediate physical goods
//

//14. Assets and capital
Q*K(+1) = Q*Z;

//15. Production function
Ym = (K)^alfa*L^(1-alfa);

//16. Wage equilibrium
W = (1-alfa)*Pm*Ym/L;

//17. Return on capital
Rk = (alfa*Pm*Ym/K + (1-delta)*Q)/Q(-1);


//
//Capital Producers
//

//18. Capital accumulation
K(+1) = (1-delta)*K + In;

//19. Optimal investment decision
Q =   1+ f + fprime*In/In(-1) - Lambda(+1)*fprime(+1)*(In(+1)/In)^2;

//20-21. Auxiliars variables
f = kappa/2*(In/In(-1)-1)^2;
fprime = kappa*(In/In(-1) - 1);


//
//Retailers
//

//22. Aggregation final goods
Ym = vp*Y;

//23. Price dispersion
vp = theta_p*(Pi(-1)^chi)^(-epsilon_p)*vp(-1) + (1-theta_p)*Pistar^(-epsilon_p);

//24. Price dynamics
1 = theta_p*(Pi(-1)^chi/Pi)^(1-epsilon_p) + (1-theta_p)*Pistar^(1-epsilon_p);

//25. Optimal price
g1 = epsilon_p/(epsilon_p-1)*g2;

//26-27. Auxiliars variables 
g1 = lambda_c*Ym*Pistar + betta*theta_p*Pistar/Pistar(+1)*(Pi^chi/Pi(+1))^(1-epsilon_p)*g1(+1);
g2 = lambda_c*Ym*Pm + betta*theta_p*(Pi^chi/Pi(+1))^(-epsilon_p)*g2(+1);


//
//Government and Central Bank
//

//28. Unemployment benefits
b = zetta*W;

//29. Unconventional Monetary Policy
psi = psi_bar + omega*(log(Rk(+1)) - log(R(+1)) - (log(steady_state(Rk))-log(steady_state(R))));

//30. Interest Rate Rule
1+i = (1+i(-1))^phi_i*(1/betta* (Pi/steady_state(Pi))^phi_pi * (Y/steady_state(Y))^phi_y)^(1-phi_i);

//31. Fisher equation
    //1+i = R(+1)*Pi(+1);
1+i(-1) = R*Pi;

//
//Equilibrium 
//

//32. Aggregate resource constraint
Y = C + In + f*In + G + tau*psi*Q*K(+1);

//33. Government consumption
G = omega_g*steady_state(Y);

@# include "commonVarEq.mod"  // added by epi-mmb team
end;


initval;
L=0.90000000;
C=1.49705506;
Z=15.28624109;
K=15.28624109;
Ym=2.29172084;
In=0.38215603;
Y=2.29172084;
R=1.01010101;
Rk=1.01260101;
Q=1.00000000;
Pm=0.76001920;
vp=1.00000000;
Pi=1.00000000;
Pistar=1.00000000;
W=1.29663748;
b=0.64831874;
i=0.01010101;
lambda_c=0.68433205;
Lambda=0.99000000;
phi=4.00000000;
eta=1.62892979;
nu=0.00403160;
Omega=3.82156027;
Gamma=1.62892979;
f=0.00000000;
fprime=0.00000000;
g1=4.52480675;
g2=3.43894000;
psi=0.00000000;
G=0.41250975;
V=6.28668141;
Q_Z=15.28624109;
V1=0.06162803;
V2=6.22505338;
Itilde=0.10000000;
@# include "commonVarSS.mod" // added by epi-mmb team
end;


endval;
L=0.99999912;
C=1.66339304;
Z=16.98469732;
K=16.98469732;
Ym=2.54635424;
In=0.42461743;
Y=2.54635424;
R=1.01010101;
Rk=1.01260101;
Q=1.00000000;
Pm=0.76001920;
vp=1.00000000;
Pi=1.00000000;
Pistar=1.00000000;
W=1.29663748;
b=0.64831874;
i=0.01010101;
lambda_c=0.61589939;
Lambda=0.99000000;
phi=4.00000000;
eta=1.62892979;
nu=0.00403160;
Omega=4.24617433;
Gamma=1.62892979;
f=0.00000000;
fprime=0.00000000;
g1=4.52480675;
g2=3.43894000;
psi=0.00000000;
G=0.45834376;
V=6.98519540;
Q_Z=16.98469732;
V1=0.06847553;
V2=6.91671987;
Itilde=0.00000088;
end;

check;
steady;

xx=[0.1260000000000000 
0.1569456000000000 
0.1927577232322560 
0.2327701953765473 
0.2755682054461491 
0.3189519488719189 
0.3601151064359844 
0.3960669137792167 
0.4242070257207728 
0.4428509932168035 
0.4514973461882993 
0.4507531490007773 
0.4420038041316804 
0.4270018016610735 
0.4075218765104307 
0.3851455254833873 
0.3611668743153598 
0.3365811478371578 
0.3121166048305015 
0.2882824909099549 
0.2654176908895255 
0.2437332646937143 
0.2233468814312737 
0.2043094713137938 
0.1866252861872388 
0.1702667250856095 
0.1551851511872727 
0.1413187051304198 
0.1285978970547742 
0.1169495690883054 
0.1062996683980133 
0.0965751550796213 
0.0877052826082385 
0.0796224247197424 
0.0722625758035439 
0.0655656177042633 
0.0594754208715694 
0.0539398295583133 
0.0489105674145362 
0.0443430900279510 
0.0401964037564295 
0.0364328648874300 
0.0330179692378328 
0.0299201394082004 
0.0271105147605684 
0.0245627476028305 
0.0222528078919588 
0.0201587979067714 
0.0182607777103221 
0.0165406017639844 
0.0149817667264581 
0.0135692702387448 
0.0122894803361713 
0.0111300150224435 
0.0100796314745927 
0.0091281243110936 
0.0082662323404398 
0.0074855532080948 
0.0067784653714497 
0.0061380568518110 
0.0055580602369400 
0.0050327934353423 
0.0045571057129402 
0.0041263285728871 
0.0037362310693449 
0.0033829791754804 
0.0030630988543791 
0.0027734425087522 
0.0025111585110745 
0.0022736635400464 
0.0020586174719800 
0.0018639005968780 
0.0016875929486262 
0.0015279555569171 
0.0013834134453091 
0.0012525402152994 
0.0011340440704977 
0.0010267551480346 
0.0009296140362869 
0.0008416613689367 
0.0007620283953707 
0.0006899284365529 
0.0006246491438224 
0.0005655454856535 
0.0005120333943258 
0.0004635840107363 
0.0004197184713122 
0.0003800031861791 
0.0003440455624725 
0.0003114901309741 
0.0002820150381569 
0.0002553288692639 
0.0002311677712650 
0.0002092928474497 
0.0001894877980640 
0.0001715567838008 
0.0001553224911293 
0.0001406243804261 
0.0001273170996595 
0.0001152690480022 
0.0001043610752176 
0.0000944853040009 
0.0000855440636600 
0.0000774489246211 
0.0000701198242310 
0.0000634842752308 
0.0000574766490861 
0.0000520375271000 
0.0000471131129007 
0.0000426547005010 
0.0000386181926754 
0.0000349636648983 
0.0000316549705326 
0.0000286593833698 
0.0000259472739877 
0.0000234918167271 
0.0000212687243926 
0.0000192560080535 
0.0000174337595716 
0.0000157839547053 
0.0000142902748439 
0.0000129379456089 
0.0000117135907276 
0.0000106050997332 
0.0000096015081833 
0.0000086928892132 
0.0000078702553505 
0.0000071254696207 
0.0000064511650649 
0.0000058406718723 
0.0000052879514096 
0.0000047875364914 
0.0000043344773039 
0.0000039242924453 
0.0000035529245993 
0.0000032167004031 
0.0000029122941135 
0.0000026366947107 
0.0000023871761156 
0.0000021612702253 
0.0000019567424999 
0.0000017715698604 
0.0000016039206773 
0.0000014521366534 
0.0000013147164208 
0.0000011903006904 
0.0000010776588067 
0.0000009756765746 
0.0000008833452387 
];

shocks;
var Itilde;
periods 1:149;
values (xx);
end;





perfect_foresight_setup(periods=149);
//perfect_foresight_solver(lmmcp);
perfect_foresight_solver(maxit=100, tolf=1e-25, stack_solve_algo=0, linear_approximation);

//write_latex_dynamic_model;
results.oo_ = oo_ ;
results.M_ = M_;
save simulated_results_base;

%@# include "saveResults.mod"   // added by epi-mmb team

