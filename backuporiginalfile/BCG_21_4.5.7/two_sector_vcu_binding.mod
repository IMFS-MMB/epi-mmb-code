var y, ypc, v1, v2, c, cpc, in, inpc, k, lambdac, lambdak, lambdai, p1, w1, w2, rk, l1, l2, a, u, pop; 


var m,   m1,  m2,  m3,  m4,  m5,  m6,  m7,  m8,  m9,  m10, m11, m12, m13, m14, m15, m16, 
    m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29, m30, m31, m32, m33, m34, m35, m36; 

var n, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, 
    n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36; 

var pop_shock, pop_shock1, pop_shock2, pop_shock3, pop_shock4, pop_shock5, pop_shock6, pop_shock7, pop_shock8, pop_shock9, pop_shock10, pop_shock11, pop_shock12, pop_shock13, pop_shock14, pop_shock15, pop_shock16, 
    pop_shock17, pop_shock18, pop_shock19, pop_shock20, pop_shock21, pop_shock22, pop_shock23, pop_shock24,  pop_shock25, pop_shock26, pop_shock27, pop_shock28, pop_shock29, pop_shock30, pop_shock31, pop_shock32, pop_shock33, pop_shock34, pop_shock35, pop_shock36;  
     
// innovations to shock processes
varexo eps_l1, eps_l2, eps_a;
varexo eps_m36, eps_n36, eps_pop_shock36;
   

parameters alphap, kappap, rhop, rhop_l, rhop_a, betap, zetap,
           deltap, omegap, nup, nu0p, gammap, chip, in_ss, l1_ss, l2_ss;


// see rbca_steadystate for parameter definitions

model;
/////////////////////////////////////////////////////////////////


// 1. lambdac
lambdac =  1/(cpc-kappap*cpc(-1))-betap*kappap/(cpc(1)-kappap*cpc)*pop(1)/pop;

// 2. lambdak 
lambdac*(1+zetap*(in-in(-1))/in(-1))-betap*lambdac(1)*(zetap*(in(1)-in)/in+zetap/2*(in(1)-in)^2/in^2) = lambdak+lambdai;

// 3. k
lambdac(1)*rk(1)*u(1)-lambdak+betap*(1-deltap)*lambdak(1) = 0;

// 4. in
k=(1-deltap)*k(-1) + in;

// 5. u
rk*k(-1) = nu0p*u^nup;

// 6. p1
w1 = gammap*p1;

// 7. v1
v1 = gammap*(l1-chip);

// 8. rk
rk = alphap*(omegap*y/v2)^(rhop/(1+rhop))*v2/(u*k(-1));

// 9. w2
w2 = (1-alphap)*(omegap*y/v2)^(rhop/(1+rhop))*v2/l2;


// 10. p1
p1 = ((1-omegap)*y/v1)^(rhop/(1+rhop));

// 11. y
y  = ( (1-omegap)^(rhop/(1+rhop))*v1^(1/(1+rhop)) + omegap^(rhop/(1+rhop))*v2^(1/(1+rhop)) )^(1+rhop);

// 12. v2
v2 = (u*k(-1))^alphap*(a*l2)^(1-alphap);

// 13. c
y = c+in+nu0p*u^(1+nup)/(1+nup);

// 15. l1
l1-l1_ss = rhop_l*(l1(-1)-l1_ss)+m+eps_l1;

// 16. l2
l2-l2_ss = rhop_l*(l2(-1)-l2_ss)+n+eps_l2;

// 17. lambdai
in = 0;

// 18. cpc
c = cpc*pop;

// 19. pop
pop - 1 = rhop_l*(pop(-1)-1)+pop_shock;

// 20. a
log(a) = rhop_a *log(a(-1))+eps_a;

in = inpc*pop;
y = ypc*pop;

// structure for setting exogenous paths for l1 and l2
m = m1(-1);
m1 = m2(-1);
m2 = m3(-1);
m3 = m4(-1);
m4 = m5(-1);
m5 = m6(-1);
m6 = m7(-1);
m7 = m8(-1);
m8 = m9(-1);
m9 = m10(-1);
m10 = m11(-1);
m11 = m12(-1);
m12 = m13(-1);
m13 = m14(-1);
m14 = m15(-1);
m15 = m16(-1);
m16 = m17(-1);
m17 = m18(-1);
m18 = m19(-1);
m19 = m20(-1);
m20 = m21(-1);
m21 = m22(-1);
m22 = m23(-1);
m23 = m24(-1);
m24 = m25(-1);
m25 = m26(-1);
m26 = m27(-1);
m27 = m28(-1);
m28 = m29(-1);
m29 = m30(-1);
m30 = m31(-1);
m31 = m32(-1);
m32 = m33(-1);
m33 = m34(-1);
m34 = m35(-1);
m35 = m36(-1);
m36 = eps_m36;


n = n1(-1);
n1 = n2(-1);
n2 = n3(-1);
n3 = n4(-1);
n4 = n5(-1);
n5 = n6(-1);
n6 = n7(-1);
n7 = n8(-1);
n8 = n9(-1);
n9 = n10(-1);
n10 = n11(-1);
n11 = n12(-1);
n12 = n13(-1);
n13 = n14(-1);
n14 = n15(-1);
n15 = n16(-1);
n16 = n17(-1);
n17 = n18(-1);
n18 = n19(-1);
n19 = n20(-1);
n20 = n21(-1);
n21 = n22(-1);
n22 = n23(-1);
n23 = n24(-1);
n24 = n25(-1);
n25 = n26(-1);
n26 = n27(-1);
n27 = n28(-1);
n28 = n29(-1);
n29 = n30(-1);
n30 = n31(-1);
n31 = n32(-1);
n32 = n33(-1);
n33 = n34(-1);
n34 = n35(-1);
n35 = n36(-1);
n36 = eps_n36;


pop_shock = pop_shock1(-1);
pop_shock1 = pop_shock2(-1);
pop_shock2 = pop_shock3(-1);
pop_shock3 = pop_shock4(-1);
pop_shock4 = pop_shock5(-1);
pop_shock5 = pop_shock6(-1);
pop_shock6 = pop_shock7(-1);
pop_shock7 = pop_shock8(-1);
pop_shock8 = pop_shock9(-1);
pop_shock9 = pop_shock10(-1);
pop_shock10 = pop_shock11(-1);
pop_shock11 = pop_shock12(-1);
pop_shock12 = pop_shock13(-1);
pop_shock13 = pop_shock14(-1);
pop_shock14 = pop_shock15(-1);
pop_shock15 = pop_shock16(-1);
pop_shock16 = pop_shock17(-1);
pop_shock17 = pop_shock18(-1);
pop_shock18 = pop_shock19(-1);
pop_shock19 = pop_shock20(-1);
pop_shock20 = pop_shock21(-1);
pop_shock21 = pop_shock22(-1);
pop_shock22 = pop_shock23(-1);
pop_shock23 = pop_shock24(-1);
pop_shock24 = pop_shock25(-1);
pop_shock25 = pop_shock26(-1);
pop_shock26 = pop_shock27(-1);
pop_shock27 = pop_shock28(-1);
pop_shock28 = pop_shock29(-1);
pop_shock29 = pop_shock30(-1);
pop_shock30 = pop_shock31(-1);
pop_shock31 = pop_shock32(-1);
pop_shock32 = pop_shock33(-1);
pop_shock33 = pop_shock34(-1);
pop_shock34 = pop_shock35(-1);
pop_shock35 = pop_shock36(-1);
pop_shock36 = eps_pop_shock36;

end;

