function [resid,rk,v1,v2,y] = function_k(k, l1, l2, betap, deltap, alphap, rhop, omegap,gammap,chip)

rk = 1-betap*(1-deltap);
v1 = gammap*(l1-chip);
v2 = k^alphap*l2^(1-alphap);
y = ((1-omegap)^(rhop/(1+rhop))*v1^(1/(1+rhop))+omegap^(rhop/(1+rhop))*v2^(1/(1+rhop)))^(1+rhop);

resid = -rk+ alphap*(omegap*y/v2)^(rhop/(1+rhop))*v2/k;