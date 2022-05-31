%%%%%%%%% Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all,

%%%%%%%%%%%% model 2 sectors, nominal frictions, roundabout

%%%%%%%%%%%(change also inside the fuction if you modify something)

%%%%%%%%% iota fixed to 1

chi=0.5;

eta=1.5;

theta=3.8;   %elasticity of sub

kappa=6;     %tail Pareto 

zmin=1;   %min Pareto

Z=1;    %aggregate productivity

fx=0.47;     %fixed costs production (check test)

thetaw=4;    %elasticity labor inputs

delta=0.00211;   %exit rate: (1-delta)^52 approx 90%

betta=0.98^(1/52);    %Weekly household discount factor

phi=4;   %Frisch elasticity labor supply

tau=(kappa/(kappa-(theta-1)))^(1/(theta-1));   %ztilde=tau*zc

alphaw=0;   %fraction wage bill in advance

nu=1;   %disutility labor multiplier

alpha=1/3;   %cobb-douglas

alphatilde=0.98;  %calvo union wage resetting, change price once per year

tr_pi=1.5;  %Taylor rule

tr_smut=0.8;  %Taylor rule

tr_y=0.5/52;  %Taylor rule

pid=7*0.002/14;    %Weekly probability of dying

pir=7*1/14-pid;    %Weekly probability of recovering  

pi1=0.05; %%% just to check steady state

pi2=0.1;  %%%% same

pi3=0.5; %%%% same

ud=10;             %Psycological costs death

psi0=1;         %%%%%%%%% normalization, fe=1 if psi1=0

psi1=1000;      %%%%%%%%%%%%%%   ratio ne*fe/GDP approx 15%

gamma=1.5;          %%%%% elasticity entry as in data?



%%%%%%%%%%%%%%%%%%%%% non-linear solver with lower and upper bounds


fun=@equilibrium;

x0=[1.5,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15];
    
lb=[0,0,0,0,0,0,0,0,0,0,0];
   
ub=[1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000];

options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',100000,'StepTolerance',1e-18,'FunctionTolerance',1e-18,'MaxIterations',10000);

y=lsqnonlin(fun,x0,lb,ub,options);

save SS_for_dynare y;


C=y(1);

C1=y(2);

lambda=y(3);

zc1=y(4);

N1=y(5);

zc2=y(6);

N2=y(7);

w=y(8);

X=y(9);

Y=y(10);

C2=y(11);


rho1=(theta/(theta-1))*(1/(Z*tau*zc1))*(((alphaw/betta+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

rho2=(theta/(theta-1))*(1/(Z*tau*zc2))*(((alphaw/betta+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

e1=(1/theta)*chi*rho1^(1-eta)*((zmin/zc1)^kappa*N1)^((theta-eta)/(1-theta))*(lambda*C)^(-eta)*Y-fx;

e2=(1/theta)*(1-chi)*rho2^(1-eta)*((zmin/zc2)^kappa*N2)^((theta-eta)/(1-theta))*(lambda*C)^(-eta)*Y-fx;

mutilde=(thetaw)/(thetaw-1);

ls=((w*lambda)/(nu*mutilde))^(1/phi);

fe1=psi0+psi1*(delta/(1-delta)*N1)^gamma;

fe2=psi0+psi1*(delta/(1-delta)*N2)^gamma;

d=0;

i=0;

s=1;

t=0;

lambdas=0;

lambdad=(betta*log(C)-betta*nu*(ls^(1+phi))/(1+phi)+betta*ud-betta*lambda*C+betta*lambda*w/mutilde*ls)/(1-betta);

lambdai=betta*pid*lambdad/(1-betta*(1-pir)+betta*pid);

lambdat=lambdai;

R=1/betta;

Ld=ls;

Ne1=N1*delta/(1-delta);

Ne2=N2*delta/(1-delta);

stilde=1;

pi=1;

wtilde=w;

f1=Ld*lambda*w^thetaw/(1-alphatilde*betta);

f2=Ld*w^thetaw*nu*ls^phi/(1-alphatilde*betta);


%%%%%%%%%%%%%%%%%%%%%%% check equilibrium is satisfied


%%%%%%%%%%%%%%%%%%% Households


a(1)=t-pi1*s*i*C1^2/(1-d)-s*i*pi2*ls*Ld-pi3*s*i;

a(2)=i-i-t+(pir+pid)*i;

a(3)=s-s+t;

a(4)=d-d-pid*i;

a(5)=nu*ls^phi-lambda*w/mutilde+lambdat*s*i/(1-d)*pi2*Ld;

a(6)=lambdat-lambdai+lambdas;

a(7)=lambdai+betta*(-lambdai*(1-pir)+pid*(lambdai-lambdad));

a(8)=lambdas+betta*(lambdat*(-i/(1-d)*pi1*C1^2-i*pi2*ls*Ld-pi3*i)-lambdas);

a(9)=lambdad+betta*(-log(C/(1-d))+nu*ls^(1+phi)/(1+phi)-ud+lambda*C/(1-d)-lambda*w/mutilde*ls-lambdad);

a(10)=C-C1^(1/(1-eta))*(1-d)^(eta/(eta-1))*chi^(1/(eta-1))*(lambda*rho1*((zmin/zc1)^kappa*N1)^(1/(1-theta))+lambdat*s*i/(1-d)*pi1*C1)^(eta/(1-eta));

a(11)=C2-(1-chi)*(1-d)^(eta)*(lambda*rho2*((zmin/zc2)^kappa*N2)^(1/(1-theta)))^(-eta)*C^(1-eta);

a(12)=fe1-betta*(1-delta)*(lambda/lambda)*(fe1+((zmin/zc1)^kappa)*e1);

a(13)=fe2-betta*(1-delta)*(lambda/lambda)*(fe2+((zmin/zc2)^kappa)*e2);

a(14)=1-betta*(lambda/lambda)*R/pi;

%%%%%%%%%%%%%%%%%%%%%%%%% Union

a(15)=wtilde-thetaw/(thetaw-1)*f2/f1;

a(16)=f1-Ld*w^(thetaw)*lambda-alphatilde*betta*(pi)^(thetaw-1)*f1;

a(17)=f2-Ld*w^(thetaw)*(nu*ls^(phi)+lambdat*s*i/(1-d)*pi2*Ld)-alphatilde*betta*(pi)^(thetaw)*f2;

%%%%%%%%%%%%%%%%%%%%%%%%% Firms 

a(18)=rho1-theta/(theta-1)*(1/(Z*tau*zc1))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

a(19)=rho2-theta/(theta-1)*(1/(Z*tau*zc2))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

a(20)=e1-(1/theta)*chi*rho1*((zmin/zc1)^kappa*N1)^(theta/(1-theta))*(lambda*rho1*((zmin/zc1)^kappa*N1)^(1/(1-theta))+lambdat*s*i/(1-d)*pi1*C1)^(-eta)*(C/(1-d))^(-eta)*Y+fx;

a(21)=e2-(1/theta)*(1-chi)*rho2*((zmin/zc2)^kappa*N2)^(theta/(1-theta))*(lambda*rho2*((zmin/zc2)^kappa*N2)^(1/(1-theta)))^(-eta)*(C/(1-d))^(-eta)*Y+fx;

%%%%%%%%%%%%%%%%%%% Entry and Exit

a(22)=N1-(1-delta)*(N1+Ne1);

a(23)=N2-(1-delta)*(N2+Ne2);

a(24)=zc1-(theta^(theta/(theta-1)))/(theta-1)*(1/Z)*(fx/(rho1^theta*((zmin/zc1)^kappa*N1)^(theta/(1-theta))*(lambda*rho1*((zmin/zc1)^kappa*N1)^(1/(1-theta))+lambdat*s*i/(1-d)*pi1*C1)^(-eta)*(C/(1-d))^(-eta)*Y))^(1/(theta-1))*chi^(1/(1-theta))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

a(25)=zc2-(theta^(theta/(theta-1)))/(theta-1)*(1/Z)*(fx/(rho2^theta*((zmin/zc2)^kappa*N2)^(theta/(1-theta))*(lambda*rho2*((zmin/zc2)^kappa*N2)^(1/(1-theta)))^(-eta)*(C/(1-d))^(-eta)*Y))^(1/(theta-1))*(1-chi)^(1/(1-theta))*(((alphaw*R+1-alphaw)*w)/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

%%%%%%%%%%%%%%%%%%%%% Taylor rule

a(26)=R/R-((pi/pi)^tr_pi*(Y/Y)^tr_y)^(1-tr_smut)*(R/R)^tr_smut;

%%%%%%%%%%%%%%%%% Aggregation

a(27)=C+Ne1*fe1+Ne2*fe2-(alphaw*R+1-alphaw)*w*Ld-e1*N1*(zmin/zc1)^kappa-e2*N2*(zmin/zc2)^kappa;

a(28)=C^eta-(1-d)^eta*(chi*rho1*((zmin/zc1)^kappa*N1)^(1/(1-theta))*(lambda*((zmin/zc1)^kappa*N1)^(1/(1-theta))*rho1+lambdat*s*i/(1-d)*pi1*C1)^(-eta)+(1-chi)*lambda^(-eta)*(rho2*((zmin/zc2)^kappa*N2)^(1/(1-theta)))^(1-eta));

a(29)=Y-C-X-Ne1*fe1-Ne2*fe2-fx*N1*(zmin/zc1)^kappa-fx*N2*(zmin/zc2)^kappa;

a(30)=X-chi*((zmin/zc1)^kappa*N1)^(1/(1-theta))/(Z*tau*zc1)*(alpha*w*(alphaw*R+1-alphaw)/(1-alpha))^(1-alpha)*(lambda*((zmin/zc1)^kappa*N1)^(1/(1-theta))*rho1+lambdat*s*i/(1-d)*pi1*C1)^(-eta)*(C/(1-d))^(-eta)*Y-(1-chi)*((zmin/zc2)^kappa*N2)^(1/(1-theta))/(Z*tau*zc2)*(alpha*w*(alphaw*R+1-alphaw)/(1-alpha))^(1-alpha)*(lambda*((zmin/zc2)^kappa*N2)^(1/(1-theta))*rho2)^(-eta)*(C/(1-d))^(-eta)*Y;

a(31)=w^(1-thetaw)-(1-alphatilde)*wtilde^(1-thetaw)-alphatilde*(w/pi)^(1-thetaw);

a(32)=(1-d)*ls-stilde*Ld;

a(33)=stilde-(1-alphatilde)*(wtilde/w)^(-thetaw)-alphatilde*(w/w)^(-thetaw)*pi^(thetaw)*stilde;

%%%%%% C2

a(34)=C-(chi^(1/eta)*C1^((eta-1)/eta)+(1-chi)^(1/eta)*C2^((eta-1)/eta))^(eta/(eta-1));

%%%%%%% fe

a(35)=fe1-psi0-psi1*(Ne1)^gamma;

a(36)=fe2-psi0-psi1*(Ne2)^gamma;



%%%%% checks

test_rhosec1=((zmin/zc1)^kappa*N1)^(1/(1-theta))*rho1;
test_rhosec2=((zmin/zc2)^kappa*N2)^(1/(1-theta))*rho2;
test_lambda=C*lambda;

test_fe=(fe1*Ne1+fe2*Ne2)/(Y-X);
test_fx=(fe1+fe2)/(2*fx);

No1=N1*(zmin/zc1)^kappa;
No2=N2*(zmin/zc2)^kappa;

income=w*ls+e1*No1+e2*No2+fe1*N1+fe2*N2;
lifetime_costs=(1/(1-betta)*ud)/(4*income)*5000;

income_2=w*ls;
lifetime_costs_2=(1/(1-betta)*ud)/(4*income_2)*2000;


function F=equilibrium(x)

%%%%%%%%% Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%(change also above if you modify something)

chi=0.5;

eta=1.5;

theta=3.8;   %elasticity of sub

kappa=6;     %tail Pareto 

zmin=1;   %min Pareto

Z=1;    %aggregate productivity

fx=0.47;     %fixed costs production (check test)

thetaw=4;    %elasticity labor inputs

delta=0.00211;   %exit rate: (1-delta)^52 approx 90%

betta=0.98^(1/52);    %Weekly household discount factor

phi=4;   %Frisch elasticity labor supply

tau=(kappa/(kappa-(theta-1)))^(1/(theta-1));   %ztilde=tau*zc

alphaw=0;   %fraction wage bill in advance

nu=1;   %disutility labor multiplier

alpha=1/3;   %cobb-douglas

psi0=1;         %%%%%%%%% normalization, fe=1 if psi1=0

psi1=1000;      %%%%%%%%%%%%%%   ratio ne*fe/GDP approx 15%

gamma=1.5;          %%%%% elasticity entry as in data?


%%%%%%%%%%%%%%%%%%%%%%%%% supporting equations

rho1=(theta/(theta-1))*(1/(Z*tau*x(4)))*(((alphaw/betta+1-alphaw)*x(8))/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

rho2=(theta/(theta-1))*(1/(Z*tau*x(6)))*(((alphaw/betta+1-alphaw)*x(8))/(1-alpha))^(1-alpha)*(1/alpha)^(alpha);

e1=(1/theta)*chi*rho1^(1-eta)*((zmin/x(4))^kappa*x(5))^((theta-eta)/(1-theta))*(x(3)*x(1))^(-eta)*x(10)-fx;

e2=(1/theta)*(1-chi)*rho2^(1-eta)*((zmin/x(6))^kappa*x(7))^((theta-eta)/(1-theta))*(x(3)*x(1))^(-eta)*x(10)-fx;

mutilde=(thetaw)/(thetaw-1);

ls=((x(3)*x(8))/(nu*mutilde))^(1/phi);

fe1=psi0+psi1*(delta/(1-delta)*x(5))^gamma;

fe2=psi0+psi1*(delta/(1-delta)*x(7))^gamma;

%%%%%%%%%%%%% non linear solver


F(1)=x(1)-x(2)^(1/(1-eta))*chi^(1/(eta-1))*(x(3)*rho1*((zmin/x(4))^kappa*x(5))^(1/(1-theta)))^(eta/(1-eta));

F(2)=x(11)-(1-chi)*(x(3)*rho2*((zmin/x(6))^kappa*x(7))^(1/(1-theta)))^(-eta)*x(1)^(1-eta);

F(3)=fe1-betta*(1-delta)*(fe1+(zmin/x(4))^kappa*e1);

F(4)=fe2-betta*(1-delta)*(fe2+(zmin/x(6))^kappa*e2);

F(5)=x(4)-(1/Z)*(fx/(((zmin/x(4))^kappa*x(5))^((theta-eta)/(1-theta))*rho1^(theta-eta)*(x(3)*x(1))^(-eta)*x(10)))^(1/(theta-1))*(chi)^(1/(1-theta))*(((alphaw/betta+1-alphaw)*x(8))/(1-alpha))^(1-alpha)*(1/alpha)^(alpha)*theta^(theta/(theta-1))/(theta-1);

F(6)=x(6)-(1/Z)*(fx/(((zmin/x(6))^kappa*x(7))^((theta-eta)/(1-theta))*rho2^(theta-eta)*(x(3)*x(1))^(-eta)*x(10)))^(1/(theta-1))*(1-chi)^(1/(1-theta))*(((alphaw/betta+1-alphaw)*x(8))/(1-alpha))^(1-alpha)*(1/alpha)^(alpha)*theta^(theta/(theta-1))/(theta-1);

F(7)=x(1)+(delta/(1-delta))*fe1*x(5)+(delta/(1-delta))*fe2*x(7)-(alphaw/betta+1-alphaw)*x(8)*ls-e1*(zmin/x(4))^kappa*x(5)-e2*(zmin/x(6))^kappa*x(7);

F(8)=x(1)^eta*x(3)^eta-chi*(x(5)*(zmin/x(4))^kappa)^((1-eta)/(1-theta))*rho1^(1-eta)-(1-chi)*(x(7)*(zmin/x(6))^kappa)^((1-eta)/(1-theta))*rho2^(1-eta);

F(9)=x(10)-x(1)-x(9)-fx*((zmin/x(4))^kappa*x(5)+(zmin/x(6))^kappa*x(7))-(delta/(1-delta))*fe1*x(5)-(delta/(1-delta))*fe2*x(7);

F(10)=x(9)-chi*((zmin/x(4))^kappa*x(5))^(1/(1-theta))/(Z*tau*x(4))*(alpha*x(8)*(alphaw/betta+1-alphaw)/(1-alpha))^(1-alpha)*(x(3)*x(1))^(-eta)*x(10)*(((zmin/x(4))^kappa*x(5))^(1/(1-theta))*rho1)^(-eta)-(1-chi)*((zmin/x(6))^kappa*x(7))^(1/(1-theta))/(Z*tau*x(6))*(alpha*x(8)*(alphaw/betta+1-alphaw)/(1-alpha))^(1-alpha)*(x(3)*x(1))^(-eta)*x(10)*(((zmin/x(6))^kappa*x(7))^(1/(1-theta))*rho2)^(-eta);

F(11)=x(1)^((eta-1)/eta)-chi^(1/eta)*x(2)^((eta-1)/eta)-(1-chi)^(1/eta)*x(11)^((eta-1)/eta);

end