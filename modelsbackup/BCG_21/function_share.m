function resid=function_share(omegap, share1, l1, l2, betap, deltap, alphap, rhop,gammap,chip,options)

kguess = 17;
k = fsolve(@(x) function_k(x,l1, l2, betap, deltap, alphap, rhop, omegap,gammap,chip),kguess,options);
  
[residk, rk,v1,v2,y] = function_k(k, l1, l2, betap, deltap, alphap, rhop, omegap,gammap,chip);
p1 = ((1-omegap)*y/l1)^(rhop/(1+rhop));
resid = share1-p1*v1/y;