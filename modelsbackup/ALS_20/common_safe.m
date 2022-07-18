% set:
% 1. p.z(:,1), S productivity
% 2. p.kappa: work sick disutility
% 3. e: prices

p.psihome = sum(p.L0.*p.psi0,'all');
[e,p] = getparams(e,p,maxit,eps);

%% POINTERS TO SAVE + 
% beginning of period: make home choice
safe.y0 = zeros(p.Ns,p.No,p.NE,p.Ne);
safe.choice = safe.y0; %1 means stay home
safe.earn = safe.y0;
safe.util = safe.y0;
safe.prod = zeros(p.Ns,p.No-1);

safe.o0 = zeros(p.NE,p.Ne);

% after home choice, SIR evolves
safe.y1 = zeros(p.Ns,p.No,p.NE,2);
safe.o1 = zeros(p.NE);

% then test
safe.tested = 0;
safe.y2 = safe.y0;
safe.o2 = safe.o0;

%% government policies; p.omeg conditional success rate
safe.tau = p.omeg * zeros(2,1);
safe.Q = 0;     %quarantine
safe.q = 0.0;   %intensity
safe.j = 0.0;   %job intensity

% time series
safe(1:p.Nt,1) = safe;

%% init steady state

%% SOLVING THE MODEL
%% safe: steady state (get initial distribution)
% initialize
Ia = 0.00000000;

safe(1).y0(:,:,1,1) = p.L0.*(1.0-p.yf);
safe(1).y0(:,:,1,2) = p.L0.*p.yf;
safe(1).o0(1,1) = p.o0;

safe(1).P = e.P0;
safe(1).price = e.price0;%

%% switch possibilities: get logit locations
% only need to run once!
p.trans = ones(p.Ns,p.No,p.No);
if runsafe==1
	[safe,p] = solvesafe(safe,Ia,p,maxit);
else
	safe = solvemodel(safe,Ia,p);
end

% maybe calibrate this?
p.trans = p.trans./365;	%*0.01;
if p.c==1
    p.trans = zeros(p.Ns,p.No,p.No);
end

safe(1:p.Nt) = safe(p.Nt);
e.P0 = safe(1).P;
e.price0 = safe(1).price;
