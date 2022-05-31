%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters, calibration targets and other settings %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta=0.96^(1/52);  %Weekly household discount factor
pid=7*0.005/18;    %Weekly probability of dying
% pid=0;    %Weekly probability of dying
pir=7*1/18-pid;    %Weekly probability of recovering
phi=phiin;          %Productivity of infected people

alpha = alphain; % Home bias
sigma = 6; % Elasticity of Substitution (from Hassan JF)

%Calibration targets for hours and income
n_target=28;         %Weekly hours
inc_target=58000/52; %weekly income

RplusD_target=0.60;

pop0 = 1; % Initial Population
popstar0 = 1;
I0=I0in;%Initial infected
Istar0=0.000; %Initial infected

% Choose which country to use for calibrating pi
pop_ini = pop0; %Initial population
I_ini = I0;     %Initial infected

HH=52*3;         %Number of periods to solve and simulate the model

use_parallel=0;

%nonlinear solver and minimizer settings
opts_fsolve=optimoptions('fsolve','Display','off','TolFun',1e-9,...
    'MaxFunctionEvaluations',1e10,...
    'MaxIterations',2000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steady State Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa=1/n_target^2;     %calculate disutility of labor parameter kappa
%so that pre-infection steady state labor is equal to
%desired target (using n=(1/kappa)^(1/2) pre-infection
%steady state equation)
z=inc_target/n_target;  %calculate parameter A such that income is equal to
%desired target, c=inc_target=A*n_target
z=z;
zstar=z;

phss = 1;
pfss = 1;

deltaHss = (alpha.^sigma*(1).^(-sigma))./...
    (alpha.^sigma*(1).^(-(sigma-1)) +...
    (1-alpha).^sigma*(1).^(-(sigma-1))).^(sigma/(sigma-1));
deltaFss = ((1-alpha).^sigma*(1).^(-sigma))./...
    (alpha.^sigma*(1).^(-(sigma-1)) + (1-alpha).^sigma*(1).^(-(sigma-1))).^(sigma/(sigma-1));





wss=z;
wstarss=zstar;

% steady states
lrss=(1/kappa)^(1/2);          %labor recovered (same as post-infection steady state)
crss=wss*lrss;                %consumption recovered = cH+cF in SS since price is 1
xrss = (alpha*(crss*deltaHss/(deltaHss+deltaFss)).^((sigma-1)/sigma) +...
    (1-alpha)*(crss*deltaFss/(deltaHss+deltaFss)).^((sigma-1)/sigma)).^(sigma/(sigma-1));

urss=log(xrss)-kappa/2*lrss^2; %utility recovered
Urss=1/(1-beta)*urss;          %PV utility recovered

liss=(1/kappa)^(1/2);          %labor infected
ciss=phi*wss*liss;             %consumption infected
xiss = (alpha*(ciss*deltaHss/(deltaHss+deltaFss)).^((sigma-1)/sigma) +...
    (1-alpha)*(ciss*deltaFss/(deltaHss+deltaFss)).^((sigma-1)/sigma)).^(sigma/(sigma-1));

uiss=log(xiss)-kappa/2*liss^2; %utility infected

LRV = 10*1e6;
LRV_per_period = LRV*(1-beta);
Udeath = 0;

Uiss=1/(1-beta*(1-pir-pid))*(uiss+beta*pir*Urss+beta*pid*Udeath);  %PV utility infected

% Check level of present value utility
if Uiss-Urss>0, error(['Error: parameterization implies Uiss>Urss: ',num2str(Uiss-Urss)]);end

% calibrate the pi's in T-function
go_calibrate_pis;

pi4 = pi4_level/(crss*deltaFss);

% policy variables
muh = zeros(HH,1);
muf = zeros(HH,1);
muhstar = zeros(HH,1);
mufstar = zeros(HH,1);

policylength = HH;
granularity = 1;
n_per = policylength/granularity;
