global changes_vec changes_names

alphap = 0.3;

kappap = 0.6;
rhop_l = 0.95;
rhop_a = 0.95;

%rho = 1/(elast-1)
rhop = 1/(2-1);   % elasticity is 1/(1+rhop)
%rho = 3;

nup  = 0.01;

%nup = 1000;
zetap = 10;


betap = 1-(0.01/3);
deltap = 0.025/3;
share1_ss = 0.27;


omegap = 0.9;
l1_ss = 0.38*0.65;
l2_ss = (1-0.38)*0.65;


% share1_ss = 0.10
% l1_ss = 0.10
% l2_ss = 0.55


%share1_ss = 0.0001;
%l1_ss = 0.0001;
%l2_ss = 0.65;

%gammap = 1/.85/(1/.85-1);  % marginal product of labor in sector 1
%chip = l1_ss/(1/.85);   % minimum scale parameter

threshold  = 0.6;  % 0.8 for sensitivity analysis
%threshold = 0.01;
gammap = 1/threshold/(1/threshold-1);  % marginal product of labor in sector 1
chip = l1_ss*threshold;   % minimum scale parameter

% threshold = 0.5;
% chip = threshold*l1_ss;
% gammap = l1_ss/(l1_ss-chip);


l1 = l1_ss;
l2 = l2_ss;

if ~isempty(changes_vec)
    nchanges = length(changes_vec);
    for this_parm = 1:nchanges
        eval([changes_names(this_parm,:),'= changes_vec(this_parm);'])
    end
end


algo = char('active-set', 'trust-region-reflective', 'interior-point', 'interior-point-convex', 'levenberg-marquardt', 'trust-region-dogleg', 'lm-line-search','sqp');
options = optimset('Display','none','Jacobian','off','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-9,'Algorithm',algo(5,:));

kguess = 5;

omegap = fsolve(@(x) function_share(x, share1_ss, l1, l2, betap, deltap, alphap, rhop,gammap,chip,options), omegap,options);


k = fsolve(@(x) function_k(x,l1, l2, betap, deltap, alphap, rhop, omegap,gammap,chip),kguess,options);
  
[resid, rk,v1,v2,y] = function_k(k, l1, l2, betap, deltap, alphap, rhop, omegap,gammap,chip);
in = deltap*k;
in_ss = in;
nu0p = rk*k;
u = 1;
c = y-in-nu0p*u^(1+nup)/(1+nup);
lambdac = 1/((1-kappap)*c)-betap*kappap/((1-kappap)*c);

lambdak = lambdac;
p1 = ((1-omegap)*y/v1)^(rhop/(1+rhop));
w1 = gammap*p1;
w2 = (1-alphap)*(omegap*y/v2)^(rhop/(1+rhop))*v2/l2;
lambdai = 0;

a = 1;

cpc = c;
inpc = in;
ypc= y;

pop = 1;

m = 0;
n = 0;
pop_shock = 0;

m1 = 0;
m2 = 0;
m3 = 0;
m4 = 0;
m5 = 0;
m6 = 0;
m7 = 0;
m8 = 0;
m9 = 0;
m10 = 0;
m11 = 0;
m12 = 0;
m13 = 0;
m14 = 0;
m15 = 0;
m16 = 0;
m17 = 0;
m18 = 0;
m19 = 0;
m20 = 0;
m21 = 0;
m22 = 0;
m23 = 0;
m24 = 0;
m25 = 0;
m26 = 0;
m27 = 0;
m28 = 0;
m29 = 0;
m30 = 0;
m31 = 0;
m32 = 0;
m33 = 0;
m34 = 0;
m35 = 0;
m36 = 0;


n1 = 0;
n2 = 0;
n3 = 0;
n4 = 0;
n5 = 0;
n6 = 0;
n7 = 0;
n8 = 0;
n9 = 0;
n10 = 0;
n11 = 0;
n12 = 0;
n13 = 0;
n14 = 0;
n15 = 0;
n16 = 0;
n17 = 0;
n18 = 0;
n19 = 0;
n20 = 0;
n21 = 0;
n22 = 0;
n23 = 0;
n24 = 0;
n25 = 0;
n26 = 0;
n27 = 0;
n28 = 0;
n29 = 0;
n30 = 0;
n31 = 0;
n32 = 0;
n33 = 0;
n34 = 0;
n35 = 0;
n36 = 0;

pop_shock1 = 0;
pop_shock2 = 0;
pop_shock3 = 0;
pop_shock4 = 0;
pop_shock5 = 0;
pop_shock6 = 0;
pop_shock7 = 0;
pop_shock8 = 0;
pop_shock9 = 0;
pop_shock10 = 0;
pop_shock11 = 0;
pop_shock12 = 0;
pop_shock13 = 0;
pop_shock14 = 0;
pop_shock15 = 0;
pop_shock16 = 0;
pop_shock17 = 0;
pop_shock18 = 0;
pop_shock19 = 0;
pop_shock20 = 0;
pop_shock21 = 0;
pop_shock22 = 0;
pop_shock23 = 0;
pop_shock24 = 0;
pop_shock25 = 0;
pop_shock26 = 0;
pop_shock27 = 0;
pop_shock28 = 0;
pop_shock29 = 0;
pop_shock30 = 0;
pop_shock31 = 0;
pop_shock32 = 0;
pop_shock33 = 0;
pop_shock34 = 0;
pop_shock35 = 0;
pop_shock36 = 0;