%% COMMON FIXED PARAMETERS
% counterfactual 1: no occ / 2: no home/fear
p.c = 0;

% sk v uk
% sk gdp drop (data): 0.0374
% sk drop: -0.047
% uk drop: -0.0410

% fear factor
p.chi = 5000;	%1500;
%p.chi = 500;

%quarantineenforcement
qe0 = 177;
qe1 = 2;

%lockdown
ld0 = 270;
ld1 = 4;

% antibody: useless
p.AB = 0; 
%p.ABt = 181;
p.ABt = Inf;

%% computation
maxit = 10000;
tol = eps;

%% # of states
p.Ns = 2;
p.No = 3;
p.NE = 3;
p.Ne = 6;
% indices
p.sick = [2;4;6];

%% # of days
p.Nt = 	365; 	% 730;	%

%% Prices
% sectoral price (for now)
e.P0 = 1.0;
e.price0 = zeros(2,3);

p.z = [1.0;1.3]*ones(1,p.No);

%% sick discounts
p.phi = ones(p.Ns,p.No);	% s, and m-w
p.kappa = 0.03;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALIBRATION PARAMETERSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should probably change this to 1?
Ia0 = 1./p.pop;

p.ygamm = ones(p.Ns,1)*gpara;
p.ydelt = ones(p.Ns,1)*0.00;
p.ym = ones(p.Ns,1)* mpara .*p.ygamm /30;

% should probably keep recovery rate equal for old?
p.ogamm = gpara/2;
p.odelt = 1.0/365 * 0.02;
p.om = mpara *p.ogamm;

%% managers care of workers' fear of infection?
p.wchi = ones(p.Ns,p.No);
p.wchi(:,3) = ones(2,1).*0;

%% testing params
p.omeg = 0.8;

% flu incidence from paper
p.yf = 0.03;%*ones(p.Ns,p.No);p.wchi
p.of = 0.03;

% symptoms: nature medicine paper and cruise
p.yeta = 0.3;   %*ones(p.Ns,p.No);
p.oeta = 0.6;
