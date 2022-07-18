// Endogenous variables
var Pi Q Qd mu Y w_n Fb nun nua nuu Qb Zb lbd Cb Cs D E mkt;
var mb ms T rot C sprqb V Gk ck Psik p_a N_a w_a Cs_a GDP unemp; 
var Cs_n Cb_n C_n C_a Cb_a;
// State variables
var Bb N_n Bg a mp G nw tau_l varsigma util_a R govt_wage transfer;
var bailout spend_G spend_tau_l spend_varsigma spend_transfer spend_govt_wage income;
var m f GDP_fix;
// Shocks
varexo ee_a ee_mp ee_G ee_nw ee_taul ee_varsigma ee_util_a;
varexo ee_govt_wage ee_transfer ee_bailout;

// Parameters
parameters pbeta peis pveps peta ppitgt pmu_eps pphitau;
parameters ptau_l pvarsigma pRtgt prho_r pphipi pphiy ptheta plambda prho_a;
parameters psigma_a prho_mp psigma_mp pbetak pbetab pgamma psig_epsa psig_epsn psig_epsu ptau_k;
parameters pLTV psdep pBgss pGss pcr prho_nw psigma_nw pYss;
parameters prho_pol psigma_pol pchi pNss;
parameters prho_covid psigma_covid pzeta psigmak pwss psigmaw_n pshare_a pphipi_a;
parameters ppref_a ppass pGDPss peis_a pkappa pmss pe_kappa_m pNass plss pYass pYnss;
parameters sho_a sho_mp sho_G sho_util_a sho_taul sho_varsigma sho_govt_wage sho_transfer sho_bailout;
parameters punempss;

load params

set_param_value('pbeta',par.beta);
set_param_value('peis',par.eis);
set_param_value('pveps',par.veps);
set_param_value('peta',par.eta);
set_param_value('ppitgt',detss.Pi);
set_param_value('pmu_eps',par.mu_eps);
set_param_value('pphitau',par.phitau);

set_param_value('ptau_l',par.tau_l);
set_param_value('pvarsigma',par.varsigma);
set_param_value('pRtgt',par.Rtgt);
set_param_value('prho_r',par.rho_r);
set_param_value('pphipi',par.phipi);
set_param_value('pphiy',par.phiy);
set_param_value('ptheta',par.theta);
set_param_value('plambda',par.lambda);
set_param_value('prho_a',par.rho_a);

set_param_value('psigma_a',par.sigma_a);
set_param_value('prho_mp',par.rho_mp);
set_param_value('psigma_mp',par.sigma_mp);
set_param_value('pbetak',par.betak);
set_param_value('pbetab',par.betab);
set_param_value('pgamma',par.gamma);
set_param_value('psig_epsa',par.sig_epsa);
set_param_value('psig_epsn',par.sig_epsn);
set_param_value('psig_epsu',par.sig_epsu);
set_param_value('ptau_k',par.tau_k);

set_param_value('pLTV',par.LTV);
set_param_value('psdep',par.sdep);
set_param_value('pcr',par.cr);
set_param_value('prho_nw',par.rho_nw);
set_param_value('psigma_nw',par.sigma_nw);
set_param_value('pYss',detss.Y);

set_param_value('prho_pol',par.rho_pol);
set_param_value('psigma_pol',par.sigma_pol);
set_param_value('pchi',par.chi);
set_param_value('pNss',detss.N_a);
set_param_value('punempss',detss.unemp);

set_param_value('prho_covid',par.rho_covid);
set_param_value('psigma_covid',par.sigma_covid);
set_param_value('pzeta',par.zeta);
set_param_value('psigmak',par.sigmak);
set_param_value('psigmaw_n',par.sigmaw_n);
set_param_value('pshare_a',par.share_a);
set_param_value('pphipi_a',par.phipi_a);
set_param_value('pe_kappa_m',par.e_kappa_m);

set_param_value('ppref_a',par.pref_a);
set_param_value('ppass',detss.p_a);
set_param_value('peis_a',par.eis_a);
set_param_value('pkappa',par.kappa);
set_param_value('pGDPss',detss.GDP);
set_param_value('pBgss',detss.Bg);
set_param_value('pGss',detss.G);
set_param_value('pwss',detss.w_a);
set_param_value('pNass',detss.N_a);
set_param_value('pmss',detss.m^2);
set_param_value('plss',par.l);
set_param_value('pYass',detss.Y_a);
set_param_value('pYnss',detss.Y_n);

set_param_value('sho_mp',shock.mp);
set_param_value('sho_a',shock.a);
set_param_value('sho_G',shock.G);
set_param_value('sho_taul',shock.tau_l);
set_param_value('sho_varsigma',shock.varsigma);
set_param_value('sho_govt_wage',shock.govt_wage);
set_param_value('sho_transfer',shock.transfer);
set_param_value('sho_bailout',shock.bailout);
set_param_value('sho_util_a',shock.util_a);

// Declare model

model; 
    // (1) Default threshold employed, N sector
    nun = exp(Bb(-1))/pchi - exp(w_n)*(1-exp(tau_l));

    // (2) Default threshold employed, A sector
    nua = exp(Bb(-1))/pchi - exp(w_a)*(1-exp(tau_l));
    
    // (3) default threshold unemployed
    nuu = exp(Bb(-1))/pchi - exp(varsigma);

    // (4) Default Rate
    exp(Fb) = exp(N_n)*normcdf(nun/psig_epsn) + exp(N_a)*normcdf(nun/psig_epsa) + (1-exp(N_n)-exp(N_a))*normcdf(nuu/psig_epsu);

    // (5) Borrower SDF
    exp(mb) = pbetab*(exp(Cb_n(-1))/exp(Cb_n))^peis;

    // (6) Borrower Euler
    exp(Qb)*(1-max(0,lbd)^2) = (exp(mb(+1)))*(1-exp(Fb(+1)));

    // (7) Borrower budget constraint
    exp(Cb_n) + exp(p_a)*exp(Cb_a) + (exp(Bb(-1))/pchi)*(1-exp(Fb)) = exp(N_n)*exp(w_n)*(1-exp(tau_l)) + exp(N_a)*exp(w_a)*(1-exp(tau_l)) + (1 - exp(N_a) - exp(N_n))*exp(varsigma) + exp(Qb)*(exp(Bb)/pchi) + transfer - T;

    // (8) Borrower constraint
    exp(Qb)*(exp(Bb)/pchi) + max(0,-lbd)^2 = pLTV;

    // (9) Bank lending
    exp(ms(+1))*(1-ptheta+ptheta*exp(mkt(+1)))*(exp(Zb(+1))/exp(Qb) - 1/exp(Qd)) = pcr*max(0,mu)^2;

    // (10) Bank borrowing
    exp(ms(+1))*(1-ptheta+ptheta*exp(mkt(+1))) = exp(mkt)*(1-max(0,mu)^2)*exp(Qd);

    // (11) Bank balance sheet
    exp(Qb)*exp(Bb) = exp(E) + bailout + exp(Qd)*exp(D);

    // (12) Bank leverage constraint
    pcr*exp(Qb)*exp(Bb) + max(0,-mu)^2  = exp(mkt)*(exp(E) + bailout);

    // (13) LoM for bank equity
    exp(E) = exp(nw)*ptheta*(exp(Zb)*exp(Bb(-1)) - exp(D(-1))) + plambda ;

    // (14) Payoff on loans
    exp(Zb) = (1-exp(Fb));

    // (15) Saver SDF
    exp(ms) = pbeta*(exp(Cs_n(-1))/exp(Cs_n))^peis;

    // (16) Deposit Euler Equation
    exp(Qd) = (1+psdep)*(exp(ms(+1)));

    // (17) Govt Bond Euler Equation
    exp(Q) = (exp(ms(+1))/exp(Pi(+1)));

    // (18) NKPC
    peta*(exp(Pi)/ppitgt)*(exp(Pi)/ppitgt-1) + pveps*((pveps-1)/pveps - exp(w_n)) = peta*exp(ms(+1))*(exp(Y(+1))/exp(Y))*(exp(Pi(+1))/ppitgt)*(exp(Pi(+1))/ppitgt-1);

    // (19) Mkt clearing, non-affected sector
    // exp(C) + exp(G) + pshare_a*exp(Psik)*exp(f(-1)) + pshare_a*(max(0,m)^2)*pkappa*(max(0,m)^2 / pmss)^pe_kappa_m = exp(Y)*(1-rot);
    exp(C_n) + exp(G) + pshare_a*exp(Psik)*exp(f(-1)) = exp(Y)*(1-rot);
    // exp(C) + exp(G) + pshare_a*exp(Psik)*exp(f(-1)) = exp(Y);

    // (20) Rotemberg resource costs
    rot = 0.5*peta*(exp(Pi)/ppitgt-1)^2;

    // (21) Output, non-affected sector
    exp(Y) = exp(a)*exp(N_n);

    // (22) Consumption, non-affected sector
    exp(C_n) = (1-pchi)*exp(Cs_n) + pchi*exp(Cb_n);

    // (23) Government BC
    exp(G) + exp(Bg(-1))/exp(Pi) + (1-exp(N_a)-exp(N_n))*exp(varsigma) + govt_wage*exp(N_a)*exp(w_a) + transfer + bailout = exp(N_n)*exp(w_n)*exp(tau_l) +  exp(N_a)*exp(w_a)*exp(tau_l) + ptau_k*(exp(Y)*(1-rot) - exp(w_n)*exp(N_n) + exp(N_a)*(exp(p_a) - exp(w_a)) - pshare_a*exp(Psik)*exp(f(-1))) + exp(Q)*exp(Bg) + T;
    // exp(G) + exp(Bg(-1))/exp(Pi) + (1-exp(N_a)-exp(N_n))*exp(varsigma) + govt_wage*exp(N_a)*exp(w_a) + transfer + bailout = exp(N_n)*exp(w_n)*exp(tau_l) +  exp(N_a)*exp(w_a)*exp(tau_l) + ptau_k*(exp(Y) - exp(w_n)*exp(N_n) + exp(N_a)*(exp(p_a) - exp(w_a)) - pshare_a*exp(Psik)*exp(f(-1))) + exp(Q)*exp(Bg) + T;

    // (24) Tax rule
    T = (exp(Bg(-1))/pBgss)^pphitau-1;

    // (25) Taylor Rule
    // 1/exp(Q) = max(1,pRtgt*((exp(Pi)/ppitgt)^pphipi)*((exp(p_a)/exp(p_a(-1)))^pphipi_a)*(((1-exp(unemp))/(1-punempss))^pphiy)*exp(mp));
    1/exp(Q) = max(1,((1/exp(Q(-1)))^prho_r)*((pRtgt*((exp(Pi)/ppitgt)^pphipi)*((exp(p_a)/exp(p_a(-1)))^pphipi_a)*(((1-exp(unemp))/(1-punempss))^pphiy))^(1-prho_r))*exp(mp));

    // (26) Wage, non-affected sector
    exp(w_n) = pwss*psigmaw_n*exp(a)*(exp(N_n)+exp(N_a))^pzeta;

    // (27) Labor mkt clearing, affected sector
    exp(N_a) = pshare_a*exp(f);

    // (28) LoM for number of firms
    exp(f)  = exp(f(-1))*exp(Gk) + max(0,m)^2;

    // (29) Firm survival rate
    exp(Gk) = normcdf(log(exp(ck))/psigmak + 0.5*psigmak);

    // (30) Costs paid by active firms, affected sector
    exp(Psik) = normcdf(log(exp(ck))/psigmak - 0.5*psigmak);

    // (31) Value of firm
    exp(ck) = exp(p_a)*exp(a) - exp(w_a) + govt_wage*exp(w_a) + exp(ms(+1))*(exp(Gk(+1))*exp(ck(+1)) - exp(Psik(+1)));

    // (32) Free-entry condition
    exp(ck) + max(0,-m)^2 = pkappa*(max(0,m)^2 / pmss)^pe_kappa_m;

    // (33) Wage, a-sector
    exp(w_a) = exp(w_n);

    // (34) Mkt clearing, a-sector
    (1-pchi)*exp(Cs_a) + (pchi)*exp(Cb_a) = exp(a)*pshare_a*exp(f);

    // (35) Saver consumption, a-sector
    exp(Cs_a) = (ppref_a*exp(util_a)*exp(Cs_n)^peis / exp(p_a))^(1/peis_a);

    exp(Cb_a) = (ppref_a*exp(util_a)*exp(Cb_n)^peis / exp(p_a))^(1/peis_a);

    // (36) TFP shock
    a = prho_a*a(-1) - psigma_a*ee_a;

    // (37) MP shock
    mp = prho_mp*mp(-1) + psigma_mp*ee_mp;

    // (38) Govt spending shock
    G = (1-prho_pol)*log(pGss) + prho_pol*G(-1) + psigma_pol*ee_G;

    // (39) Net worth shock
    nw = prho_nw * nw(-1) + psigma_nw*ee_nw;

    // (40) Payroll Tax
    tau_l = (1-prho_pol)*log(ptau_l) + prho_pol*tau_l(-1) - psigma_pol*ee_taul;

    // (41) UI
    varsigma = (1-prho_pol)*log(pvarsigma) + prho_pol*varsigma(-1) + psigma_pol*ee_varsigma;

    // (42) Spread
    exp(sprqb) = 1/exp(Qb) - 1/exp(Q);

    // (43) Bank Value
    exp(V) = exp(E)*exp(mkt);

    // (44) Utility shock
    util_a = prho_covid*util_a(-1) - psigma_covid*ee_util_a;

    // (45) Definition of GDP
    exp(GDP) = exp(a)*(exp(N_n) + exp(p_a)*exp(N_a));

    // (46) Interest rate
    exp(R) = 1/exp(Q);

    // (47) Liquidity assistance
    govt_wage = prho_pol*govt_wage(-1) + psigma_pol*ee_govt_wage;

    // (48) Transfer
    transfer  = prho_pol*transfer(-1)  + psigma_pol*ee_transfer;

    // (49) Bank bailout
    bailout = prho_pol*bailout(-1)  + psigma_pol*ee_bailout;

    // Deficit variables
    spend_G         = exp(G);
    spend_tau_l     = -(exp(N_n)*exp(w_n) +  exp(N_a)*exp(w_a))*exp(tau_l);
    spend_varsigma  = (1-exp(N_a)-exp(N_n))*exp(varsigma);
    spend_transfer  = transfer;
    spend_govt_wage = govt_wage*exp(w_a)*exp(N_a);

    // Labor income net of transfers
    income       = exp(N_n)*exp(w_n)*(1-exp(tau_l)) + exp(N_a)*exp(w_a)*(1-exp(tau_l)) + (1-exp(N_a)-exp(N_n))*exp(varsigma) + transfer;
    exp(GDP_fix) = exp(GDP)/(pYnss + exp(p_a)*pYass);
    exp(unemp)   = 1 - exp(N_a) - exp(N_n);
    exp(Cb) = exp(Cb_n) + exp(p_a)*exp(Cb_a);
    exp(Cs) = exp(Cs_n) + exp(p_a)*exp(Cs_a);
    exp(C) = exp(C_n) + exp(p_a)*exp(C_a);
    exp(C_a) = pchi*exp(Cb_a) + (1-pchi)*exp(Cs_a);
end;


%% Declare the exogenous shocks: they are iid, standard normal
shocks;
    var ee_a;
    periods 1;
    values (sho_a);
    var ee_mp;
    periods 1;
    values (sho_mp);
    var ee_bailout;
    periods 1;
    values (sho_bailout);
    var ee_taul;
    periods 1;
    values (sho_taul);
    var ee_G;
    periods 1:3;
    values (sho_G);
    var ee_varsigma;
    periods 1:3;
    values (sho_varsigma);
    var ee_govt_wage;
    periods 1:3;
    values (sho_govt_wage);
    var ee_transfer;
    periods 1:3;
    values (sho_transfer);
    var ee_util_a;
    periods 1:5;
    values (sho_util_a);
end;


%% Initial values correspond to the steady state
resid(1);
check;
model_diagnostics;
steady;


%% Simulate the model
%simul(periods=100,maxit=100,stack_solve_algo=0);
%simul(periods=100,maxit=100,stack_solve_algo=1);
perfect_foresight_setup(periods=100);
perfect_foresight_solver(stack_solve_algo=7, solve_algo=9, maxit=100);