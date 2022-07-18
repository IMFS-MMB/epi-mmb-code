%%% This file calibrates the model and computes shocks that match the path
%%% of US unemployment
clear
clc
close all

par       = struct; % structure for parameters
detss     = struct; % structure to save deterministic steady state

%% (1) Calibration & Steady State
% Standard parameters
par.beta     = (1/1.02)^(1/4);    % discount factor for saver
par.betak    = (1/1.02)^(1/4);    % discount factor for bank  
par.betab    = 0.98*par.betak;    % borrower impatience relative to bank
par.eis      = 1.00;              % EIS
par.veps     = 6;                 % CES
par.calvo    = 0.75 ;
par.eta      = (par.veps-1)*par.calvo/((1-par.calvo)*(1-par.beta*par.calvo)); % This ensures that the slope of the mktllips Curve is the same up to a 1st order
par.pitgt    = 1.00^(1/4);         % inflation target
par.chi      = 0.45;               % fraction of borrowers         
par.share_a  = 0*0.30 + 0.40;      % employment share of affected sector
par.eis_a    = 0.50;

% Labor market
par.l         = 1 - 0.03833;          % employment at steady state
par.defk      = 1 - par.l;            % firm default/exit rate at SS
par.mu_eps    = 0.00;                 % Mean of liquidity shocks
par.zeta      = 0.05;                 % elasticity of wage to tightness
par.markupk   = 1.10;                 % markup in the services sector
par.kappa     = 0.531;                % cost of entry in the services sector
par.e_kappa_m = 1.0;                  % congestion elasticity 

% Household finance
par.pti = 0.30;    % = gamma*B/w*l, implies gamma
par.Fb  = 0.08/4;  % implies Fb parameter
par.Fbu = 0.40/4;  % Default rate for unemployed

% Government
par.goy        = 0.20;               % spending equals 20% of GDP
par.bgoy       = 4*0.90;             % Govt debt/GDP
par.phitau     = 0.01;               % Speed of adjustment for taxes
par.tau_l      = 0.15;
par.varsigma_w = 0.325;              % percentage of unemployed who receive insurance
par.sdep       = 1.00^(1/4)-1;       % deposit wedge vs govt debt
par.Rtgt       = par.pitgt/par.beta; % interest rate target
par.rho_r      = 0.80;               % Taylor Rule persistence
par.phipi      = 2.00;               % Taylor rule parameter on inflation
par.phiy       = 0.25;               % Taylor rule parameter on unemployment
par.phipi_a    = 0.00;

% Banks
par.gzspr    = 1.01^(1/4); % target for the lending spread (GZ spread) Zb/Qb / (1/Qd) --> implies betak
par.cr       = 0.10  ;     % capital requirement
par.theta    = 0.90;

% Shocks
par.abar    = 1;       % mean of TFP shock
par.rho_a   = 0.75;    % persistence of TFP shock
par.sigma_a = 0.01;    % standard deviation of TFP shock

par.mpbar    = 1;      % Monetary policy shock
par.rho_mp   = 0.75;
par.sigma_mp = 0.01;

par.rho_nw = 0.75;     % net worth shock on banks
par.sigma_nw = 0.01;

par.rho_pol = 0.00;    % CARES Act policies
par.sigma_pol = 1;

par.rho_covid   = 0.00;   % COVID shock
par.sigma_covid = 1;

% Inflation and risk-free rate
detss.Pi  = par.pitgt;
detss.Q   = par.beta/detss.Pi;     % Risk-free rate
detss.Qd  = par.beta*(1+par.sdep); % Price of deposits
detss.unemp = 1 - par.l;           % unemployment

% Banking: avg return, mkt value, and multiplier on constraint
par.ret = (1/detss.Qd)*(par.gzspr-1);
fmu = @(x) (par.betak).*(1-par.theta+par.theta*x).*par.ret/par.cr;
f1  = @(x) x.*(1-fmu(x)).*detss.Qd - (par.betak).*(1-par.theta+par.theta*x);
out = fsolve(f1, 1.1);
detss.mkt = out; 
detss.mu  = sqrt(fmu(out));

% Labor markets
detss.N_n = (1-par.share_a)*par.l;
detss.Y   = detss.N_n;
detss.N_a = par.share_a*par.l;
detss.w_n = (par.veps-1)/par.veps;
detss.w_a = detss.w_n;
detss.f   = detss.N_a/par.share_a; % number of firms
detss.p_a = par.markupk*detss.w_a;
par.sigmaw_n = (1/(detss.N_n+detss.N_a))^par.zeta; % constant for wage rule
par.varsigma = par.varsigma_w*detss.w_n;           % unemployment subsidy


detss.ck = par.kappa; % free-entry condition imposes this

% Cost shock parameters, target entry rate in the services sector
fGk   = @(x) normcdf(log(detss.ck)/x(1) + x(1)/2);
fPsik = @(x) normcdf(log(detss.ck)/x(1) - x(1)/2);
f1 = @(x) detss.p_a  - detss.w_a + par.beta*(fGk(x)*detss.ck - fPsik(x)) - detss.ck;
% f2 = @(x) fGk(x) - detss.Gk;
% F  = @(x) [f1(x); f2(x)];
% x0  = [0.05; 0.05];
x0 = [4.9675];
fsolveopts = optimoptions('fsolve','OptimalityTolerance',1e-8);
out        = fsolve(f1, x0, fsolveopts);
par.sigmak = out;
detss.Psik = fPsik(out);
detss.Gk   = fGk(out);
detss.m    = sqrt(detss.f*(1-detss.Gk)); % entry rate
fprintf('Number of entrants is %4.4f\n', detss.m^2);

% Expenditure components
detss.G = par.goy*detss.Y;
%detss.C = detss.Y - detss.G - par.share_a*detss.Psik*detss.f - par.kappa*(detss.m^2)*par.share_a; 
detss.C = detss.Y - detss.G - par.share_a*detss.Psik*detss.f ; 

% Borrower, household finance
detss.Fb  = par.Fb;
% Household debt, based on PTI ratio
detss.Bb     = par.pti*par.chi*(detss.w_n*detss.N_n*(1-par.tau_l) + detss.w_a*detss.N_a*(1-par.tau_l) + par.varsigma*(1-par.l));
% Set debt maturity to 1 for the baseline
% par.gamma = par.pti*detss.w*detss.N*(1-par.tau_l)/detss.Bb;
par.gamma = 1.0;

detss.nuu = par.gamma*detss.Bb/detss.Pi/par.chi - par.varsigma;
detss.nun = par.gamma*detss.Bb/detss.Pi/par.chi - detss.w_n*(1-par.tau_l);
detss.nua = par.gamma*detss.Bb/detss.Pi/par.chi - detss.w_a*(1-par.tau_l);

par.sig_epsu = (detss.nuu - par.mu_eps)/norminv(par.Fbu);
par.sig_epsa = (detss.nun - par.mu_eps)/norminv((par.Fb - (1-detss.N_a-detss.N_n)*par.Fbu)/(detss.N_a+detss.N_n));
par.sig_epsn = par.sig_epsa;

detss.Qb = par.gamma*(1-par.Fb)/(par.gzspr/detss.Qd - (1-par.gamma)*(1-par.Fb));
detss.Zb = (1-par.gamma)*detss.Qb*(1-par.Fb) + par.gamma*(1-par.Fb);

par.LTV   = detss.Qb*detss.Bb*(1-(1-par.gamma)*(1-par.Fb)/detss.Pi)/par.chi;
detss.lbd = sqrt(1 - (1/detss.Qb)*((par.betab/detss.Pi)*par.gamma*(1-detss.Fb))/(1-(par.betab/detss.Pi)*(1-par.gamma)*(1-par.Fb)));

% Consumption and preference weights
detss.G    = par.goy*detss.Y;
detss.Cb   = detss.N_n*detss.w_n*(1-par.tau_l) + detss.N_a*detss.w_a*(1-par.tau_l) + (1-detss.N_n-detss.N_a)*par.varsigma + detss.Qb*detss.Bb/par.chi - (detss.Bb/detss.Pi/par.chi)*(par.gamma*(1-detss.Fb) + (1-par.gamma)*(1-par.Fb)*detss.Qb);
detss.Cs   = (detss.C - par.chi*detss.Cb)/(1-par.chi);
par.pref_a = (detss.N_a/(1-par.chi))^par.eis_a * (detss.p_a/detss.Cs^par.eis);

detss.C    = par.chi*detss.Cb + (1-par.chi)*detss.Cs;
detss.Cs_a = (par.pref_a*(detss.Cs^par.eis)/detss.p_a)^(1/par.eis_a);

% Banking quantities
detss.E    = par.cr*detss.Qb*detss.Bb/detss.mkt;
detss.D    = (detss.Qb*detss.Bb-detss.E)/detss.Qd;
par.lambda = detss.E*(detss.Pi-par.theta/detss.Qd)/detss.Qb/detss.Bb - par.theta*par.ret;

% Government
detss.rot = 0;
detss.Bg  = par.bgoy*detss.Y;
par.tau_k = (detss.G + detss.Bg*(1/detss.Pi - detss.Q) + (1-detss.N_n-detss.N_a)*par.varsigma - (detss.N_n*detss.w_n+detss.N_a*detss.w_a)*par.tau_l)/(detss.Y*(1-detss.rot) - detss.w_n*detss.N_n + detss.N_a*(detss.p_a-detss.w_a) - par.share_a*detss.Psik*detss.f);

% Shocks at the steady state
detss.a  = par.abar;
detss.mp = par.mpbar;
detss.nw = 1;

% Other variables that are useful to define
detss.mb = par.betab; %SDF
detss.ms = par.beta;  % SDF
detss.T  = 0;         % Lump-sum tax
detss.tau_l     = par.tau_l;
detss.varsigma  = par.varsigma;
detss.LTV       = par.LTV;
detss.sprqb     = 1/detss.Qb - 1/detss.Q;
par.lambda      = par.lambda*detss.Qb*detss.Bb;
detss.V         = detss.E*detss.mkt;
detss.util_a    = 1;
detss.GDP       = detss.Y + detss.p_a*detss.N_a;
detss.R         = 1/detss.Q;
detss.govt_wage = 0;
detss.transfer  = 0;
detss.bailout   = 0;
detss.spend_G   = detss.G ;
detss.spend_tau_l       = - (detss.N_n*detss.w_n+detss.N_a*detss.w_a)*par.tau_l;
detss.spend_varsigma    = par.varsigma*(1-detss.N_a-detss.N_n);
detss.spend_transfer    = 0;
detss.spend_govt_wage   = 0;
detss.income            = (1-par.tau_l)*(detss.w_a*detss.N_a + detss.w_n*detss.N_n) + (1-detss.N_a-detss.N_n)*par.varsigma;
detss.GDP_fix = 1;
detss.Y_a = detss.N_a;
detss.Y_n = detss.Y;

save detss.mat detss
save params.mat par detss

%% (2) Find util_a shocks to match path of unemployment given fiscal policy response

% Data, 1 = 2020Q1, 2 = 2020Q2, etc
% Observations through 202Q3 are from FRED (Data/Unemployment/UNRATE.csv).
% Median SPF forecasts used after that, through 2021Q2
unemp_data = [0.03833; 0.1303; 0.0883; 0.068; 0.070; 0.067];

% Flow percentages of ANNUAL GDP - Equivalent to old percentages/4 
% From Data/CARES Act/cares_spending_v3.xlsx
annual_GDP = 4*detss.GDP;

% data_G         = [0.0; 0.009748778; 0.001039286; 0.0];
data_varsigma  = [0.0; 0.009380556; 0.007786111; 0.0];
data_govt_wage = [0.0; 0.012556429; 0.014181524; 0.000462/100];
data_transfer  = [0.0; 0.012834921; 0.000184921; 0.0];
data_G         = [0.0; 0.023/3; 0.023/3; 0.023/3];     % Data on G is low quality, use CARES Law forecast

% Initial guesses for the shocks, only restriction is shock.util_a > 0
shock = struct;
shock.a         = 0;
shock.mp        = 0;
shock.G         = 0*(0.023 + 0*0.0107)*annual_GDP/detss.G/4 + 0*(0.023)*annual_GDP/detss.G/3;
shock.tau_l     = 0;
shock.varsigma  = 0*(0.012 + 0*0.017)*annual_GDP/(par.varsigma*(1-detss.N_a-detss.N_n))/4 + 0*(0.017)*annual_GDP/(par.varsigma*(1-detss.N_a-detss.N_n))/3;
shock.govt_wage = 0*(0.020 + 0*0.0267)*annual_GDP/(detss.N_a*detss.w_a)/4 + 0*(0.0267)*annual_GDP/(detss.N_a*detss.w_a)/3;
shock.transfer  = 0*(0.012)*annual_GDP/3;
shock.bailout   = 0;
shock.util_a    = (1.0^(1/par.eis_a)); % initial guess for pandemic shock
save params.mat par detss shock

% Run Dynare once
dynare model.mod noclearall;

% load data_shocks
% oo_.exo_simul = data_shocks;
% perfect_foresight_solver;
% 
% return

% Indices for variables of interest
ind_unemp  = strmatch('unemp', M_.endo_names, 'exact');
ind_ua_var = strmatch('util_a', M_.endo_names, 'exact');
ind_ua     = strmatch('ee_util_a', M_.exo_names, 'exact');

ind_G         = strmatch('ee_G', M_.exo_names, 'exact');
ind_varsigma  = strmatch('ee_varsigma', M_.exo_names, 'exact');
ind_govt_wage = strmatch('ee_govt_wage', M_.exo_names, 'exact');
ind_transfer  = strmatch('ee_transfer', M_.exo_names, 'exact');

% Load the fiscal shocks as in the data
oo_.exo_simul(1:4,ind_G)         = data_G*annual_GDP/detss.G;
oo_.exo_simul(1:4,ind_varsigma)  = data_varsigma*annual_GDP/(par.varsigma*(1-detss.N_a-detss.N_n));
oo_.exo_simul(1:4,ind_govt_wage) = data_govt_wage*annual_GDP/(detss.N_a*detss.w_a);
oo_.exo_simul(1:4,ind_transfer)  = data_transfer*annual_GDP;

% These are the "old" preference shocks to be updated
old_shocks = oo_.exo_simul(1:6,ind_ua);

% Control for the while loop
dist    = 1000;
tol     = 1e-4;
maxiter = 50;
iter    = 0;

% Loop to update shocks
while dist > tol && iter < maxiter
    
    % Run solver given a guess for the exogenous shocks in oo_.exo_simul
    perfect_foresight_solver;
    
    % Compute path of unemployment in the model
    unemp_model = exp(oo_.endo_simul(ind_unemp,1:end)');
    
    % Difference with respect to data (20Q1-21Q2)
    diff_2020 = unemp_data(1:6) - unemp_model(1:6);
    
    % Distance measure and iteration count
    dist      = sum(diff_2020.^2);
    iter      = iter + 1;
    
    % Adjust shocks to match path of unemployment given data
    if dist > tol
        
        adjustment = max(-0.50, min(0.50, diff_2020./unemp_data(1:6)));
        
        % Update shock series
        new_shocks = old_shocks.*(1+adjustment);
        oo_.exo_simul(1:6,ind_ua) = new_shocks;
        
        old_shocks = new_shocks;
    end
    
    % print loop info
    fprintf('-------------------\n')
    fprintf('Iteration %4.0f\n', iter)
    fprintf('Distance: %4.6f\n', dist)  
end

% Save the shock
data_shocks = oo_.exo_simul;
save data_shocks.mat data_shocks

% Plot data vs. model
Tirf = 9;
time = (1:Tirf)';
% dt_start = datetime(2019,12,31);
dt_start = datetime(2020,03,31);
dt_end   = dt_start+calmonths(3*Tirf-1);
dateplot = dt_start:calmonths(3):dt_end;

unemp_data = [unemp_data; NaN(Tirf-length(unemp_data),1)];

figure
subplot(2,1,1)
plot((1:Tirf), 100*[unemp_model(1:Tirf), unemp_data], 'Linewidth', 2)
axis tight
grid minor
title('Unemployment Rate, %')
legend('Model', 'Data', 'Location', 'Northeast')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21', 'Q2-21', 'Q3-21', 'Q4-21', 'Q1-22'})

subplot(2,1,2)
plot((1:Tirf), 100*(exp(oo_.endo_simul(ind_ua_var,1:Tirf)')/detss.util_a-1), 'Linewidth', 2)
axis tight
grid minor
title('Shock, % deviation from SS')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21', 'Q2-21', 'Q3-21', 'Q4-21', 'Q1-22'})

print -depsc figures/shock_model_data
print -dpng figures/shock_model_data


% Plot fiscal policy impulses
figure
subplot(2,2,1)
plot((1:5), 100*[data_G; 0], '-s', 'Linewidth', 2)
axis tight
grid minor
title('G spending, % of GDP')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21'})

subplot(2,2,2)
plot((1:5), 100*[data_transfer; 0], '-s', 'Linewidth', 2)
axis tight
grid minor
title('Tb spending, % of GDP')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21'})

subplot(2,2,3)
plot((1:5), 100*[data_varsigma; 0], '-s', 'Linewidth', 2)
axis tight
grid minor
title('ui spending, % of GDP')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21'})

subplot(2,2,4)
plot((1:5), 100*[data_govt_wage; 0], '-s', 'Linewidth', 2)
axis tight
grid minor
title('Ta spending, % of GDP')
% datetick('x','QQ-YY','keeplimits','keepticks')
xticklabels({'Q1-20', 'Q2-20', 'Q3-20', 'Q4-20', 'Q1-21'})

print -depsc figures/fiscal_policy_data
print -dpng figures/fiscal_policy_data

%% (end) Delete temporary files and other Dynare-generate clutter
delete *.log *_dynamic.m *_results.mat *_variables.m *_static.m
rmdir('model', 's')
rmdir('+model', 's')