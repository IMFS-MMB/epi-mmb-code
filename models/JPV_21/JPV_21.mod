// infections.mod
// Jones, Philippon, Venkateswaran (2021)
//

var 
    c1 l1 m1 chi1t mbar1
    Vmbar lambda lambdae lambdas lambdad lambdai Vi Vd Vs e L1 C1 M1 
	in d s delta f Vf a Va abar 
    r;                         // added by epi-mmb team
  
@#include "commonVar.mod"    // added by epi-mmb team

varexo	
	ee_i0 ee_vac ee_e0 ee_d ;

parameters 
    beta ec1 el1 chibar1 
	N gamma deltabar phi kappa rho e0 eta  
    rhom Deltachi us ud psim deltascale
    rhof omega1 omega2 rhoa varphi1 varphi2 amean;

load input ;
load inf_ini
ee_T_i(1)=helper;

beta		= params(1) ;
ec1			= params(2) ;
el1			= params(3) ;
chibar1		= params(4) ;
N			= params(5) ;
gamma		= params(6) ;
deltabar	= params(7) ;
phi         = params(8) ;
kappa		= params(9) ;
rho			= params(10) ;
e0			= params(11) ;
eta			= params(12) ;
rhom		= params(13) ;
Deltachi	= params(14) ;
us			= params(15) ;
ud			= params(16) ;
psim		= params(17) ;
deltascale	= params(18) ;
rhof        = params(19) ;
omega1      = params(20) ;
omega2      = params(21) ;
rhoa        = params(22) ;
varphi1     = params(23) ;
varphi2     = params(24) ;
amean       = params(25) ;

model ;

    // c1 FOC
    1/c1 = lambda + lambdae*ec1*C1 ;

    // l1 FOC
    l1^eta = lambda - lambdae*(1-m1)*el1*(1-M1)*L1 ;

    // m1 FOC
    lambda*chi1t*m1 + omega1*f^omega2*M1 = 1/(1-d(-1)-kappa*in(-1)) * (beta * Vmbar(+1) + beta*Vf(+1)) + lambdae*el1*l1*(1-M1)*L1 ;
    
    // e FOC
    lambdae = (1-ee_vac)*(lambdai - lambdas)*gamma*in(-1)*s(-1) ;

    // i FOC
    lambdai = -beta*Vi(+1) ;
    
    // d FOC
    lambdad = -beta*Vd(+1) ;

    // s FOC
    lambdas = -beta*Vs(+1) ;

    // Vi
    Vi = kappa * l1^(1+eta) / (1+eta) - us*kappa - kappa*lambda*(l1 - chi1t/2 * m1^2) + 
        kappa*lambdae*(1-m1)*el1*l1*L1*(1-M1) - 
        (1-rho-delta*kappa)*lambdai - 
        delta*kappa*lambdad ;

    // Vd
    Vd = l1^(1+eta) / (1+eta) - ud - log(c1) - lambda*(l1 - chi1t/2 * m1^2 - c1) + lambdae*(ec1*c1*C1 + (1-m1)*el1*l1*L1*(1-M1)) - lambdad ;

    // Vs
    Vs = -(1-ee_vac)*lambdai*gamma*e*in(-1) - lambdas*(1-(1-ee_vac)*gamma*e*in(-1)) ;

    // Vmbar
    Vmbar = beta*psim*Vmbar(+1) + (1/2)*lambda*(1-d(-1)-kappa*in(-1))*chibar1*Deltachi*rhom*exp(-rhom*mbar1(-1))*m1^2 ;

    // chi1t and mbar
    chi1t = chibar1*(1-Deltachi*(1-exp(-rhom*mbar1(-1)))) ;
    mbar1 =  psim*mbar1(-1) + m1 ;

    // goods market
    c1 = (1-d(-1)-kappa*in(-1))*(l1-chi1t/2*m1^2)/(1-d(-1)) ;

    // e 
    e = e0*(1-ee_e0)*(1-a) + (1-d(-1))*ec1*c1*C1 + (1-d(-1)-kappa*in(-1))*(1-m1)*el1*l1*L1*(1-M1) ;

    C1 = (N-d(-1)) * c1 ;
    L1 = (N-kappa*in(-1)-d(-1)) * l1 ;

    M1 = m1 ;

    // f
    f = 0;//rhof*f(-1) + M1(-1) ;

    // Vf
    Vf = 0;//-(1-d(-1)-kappa*in(-1))*omega1/2*omega2*f^(omega2-1)*M1^2 + beta*rhof*Vf(+1) ;

    // a
    0 = a;//-(1-d(-1))*varphi1*abar^varphi2*a + lambdae*e0*(1-ee_e0) + beta*Va(+1) ;

    // Va
    Va = 0;//-(1-d(-1))*varphi1*varphi2/2*abar^(varphi2-1)*a^2 + beta*rhoa*Va(+1) ;

    // abar
    abar = (1-rhoa)*amean + rhoa*abar(-1) + rhoa*a(-1) ;

    //*****************//
    // Contagion Block //
    //*****************//

    in = in(-1) + (1-ee_vac)*gamma * e * in(-1) * s(-1) - rho * in(-1) - delta * kappa * in(-1) + ee_i0  ;
    d = d(-1) + delta * kappa * in(-1)   ;
    s = s(-1) - (1-ee_vac)*gamma * e * in(-1) * s(-1) - ee_i0 ;    
    delta = (1-ee_d)*deltascale*(deltabar + (exp(phi * kappa * in(-1))-1)) ;
    r = 1 - in - s - d;         // added by epi-mmb team

@# include "commonVarEq.mod"  // added by epi-mmb team
end;


initval;

c1			= xinitial(1) ;
l1			= xinitial(2) ;
m1			= xinitial(3) ;
chi1t		= xinitial(4) ;
mbar1		= xinitial(5) ;
Vmbar		= xinitial(6) ;
lambda		= xinitial(7) ;
lambdae		= xinitial(8) ;
lambdas		= xinitial(9) ;
lambdad		= xinitial(10) ;
lambdai		= xinitial(11) ;
Vi			= xinitial(12) ;
Vd			= xinitial(13) ;
Vs			= xinitial(14) ;
e			= xinitial(15) ;
L1			= xinitial(16) ;
C1			= xinitial(17) ;
M1			= xinitial(18) ;
in			= xinitial(19) ;
d			= xinitial(20) ;
s 			= xinitial(21) ;
delta 		= xinitial(22) ;
f           = xinitial(23) ;
Vf          = xinitial(24) ;
a           = xinitial(25) ;
Va          = xinitial(26) ;
abar        = xinitial(27) ;
r           = 1 - xinitial(19) - xinitial(20) - xinitial(21);   // added by epi-mmb team

@# include "commonVarSS.mod"
end;

resid 

endval ;

c1			= xfinal(1) ;
l1			= xfinal(2) ;
m1			= xfinal(3) ;
chi1t		= xfinal(4) ;
mbar1		= xfinal(5) ;
Vmbar		= xfinal(6) ;
lambda		= xfinal(7) ;
lambdae		= xfinal(8) ;
lambdas		= xfinal(9) ;
lambdad		= xfinal(10) ;
lambdai		= xfinal(11) ;
Vi			= xfinal(12) ;
Vd			= xfinal(13) ;
Vs			= xfinal(14) ;
e			= xfinal(15) ;
L1			= xfinal(16) ;
C1			= xfinal(17) ;
M1			= xfinal(18) ;
in			= xfinal(19) ;
d			= xfinal(20) ;
s 			= xfinal(21) ;
delta 		= xfinal(22) ;
f           = xfinal(23) ;
Vf          = xfinal(24) ;
a           = xinitial(25) ;
Va          = xinitial(26) ;
abar        = xinitial(27) ;
r           = 1 - xfinal(19) - xfinal(20) - xfinal(21);  //added by epi-mmb team

@# include "commonVarSS.mod" //added by epi-mmb team
end ;

shocks;

    var ee_i0 ;
        periods 1:1000 ;
        values (ee_T_i) ;

    var ee_vac ;
        periods 1:1000 ;
        values (ee_T_vac) ;

    var ee_e0 ;
        periods 1:1000 ;
        values (ee_T_e0) ;

    var ee_d ;
        periods 1:1000 ;
        values (ee_T_d) ;

end;


resid 

model_diagnostics;

simul(periods=1000, maxit=10,noprint) ;
//forecast(periods=20);

//save allresults_infections;
results.oo_ = oo_ ;
results.M_ = M_;
save simulated_results_base;
@# include "saveResults.mod"   //added by epi-mmb team
