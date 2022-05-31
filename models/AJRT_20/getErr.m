function [err,ls,lsstar,li,listar,lr,lrstar,xs,xsstar,xi,xistar,xr,xrstar,chs,cfs,chsstar,cfsstar,chi,cfi,chistar,cfistar,chr,cfr,chrstar,cfrstar,...
    ph,pfstar,S,Sstar,I,Istar,R,Rstar,T,Tstar,aggC,aggCstar,aggL,aggLstar,g,gstar,aggU,aggUstar,pop,popstar] = ...
    getErr(guess,HH,policylength,beta,z,zstar,kappa,alpha,sigma,I0,Istar0,pop0,popstar0,pi1,pi2,pi3,pi4,pir,pid,pid_accel,phi,Uiss,Urss,Udeath,muh,muf,muhstar,mufstar,delta_mu,delta_nu)

rho = 1;

% guesses for labor
ls=guess(1:HH);
lsstar=guess(HH+1:2*HH);
chs=guess(2*HH+1:3*HH);
cfs=guess(3*HH+1:4*HH);
cfsstar=guess(4*HH+1:5*HH);
chsstar=guess(5*HH+1:6*HH);

if (sum(ls<=0)+sum(lsstar<=0)+sum(chs<=0)+sum(cfsstar<=0)+...
        sum(cfs<=0)+sum(chsstar<=0)>0)
    err = 1e8*ones(7*HH,1);
    return;
end

% guesses for pfstar
pfstar=guess(6*HH+1:7*HH);
% initialize ph normalized to one
ph=ones(HH,1);


%%%%%%%%%%%%%%%%%%%%%%%%
% equilibrium equations 
%%%%%%%%%%%%%%%%%%%%%%%%

% wages
w = ph.*z;
wstar = pfstar.*zstar;


% consumption of susceptible
% we conjectured chs, cfsstar. 
% Solve for g and gstar by budget constraint
g = (1+muf).*pfstar.*cfs + (1+muh).*ph.*chs - w.*ls;
gstar = (1+muhstar).*ph.*chsstar + (1+mufstar).*pfstar.*cfsstar - wstar.*lsstar;


xs = (alpha*(chs).^((sigma-1)/sigma) + (1-alpha)*(cfs).^((sigma-1)/sigma)).^(sigma/(sigma-1));
xsstar = (alpha*(cfsstar).^((sigma-1)/sigma) + (1-alpha)*(chsstar).^((sigma-1)/sigma)).^(sigma/(sigma-1));


% spending shares for I and R households
deltaH = (alpha.^sigma*((1+muh).*ph).^(-sigma))./...
    (alpha.^sigma*((1+muh).*ph).^(-(sigma-1)) +...
    (1-alpha).^sigma*((1+muf).*pfstar).^(-(sigma-1))).^(sigma/(sigma-1));
deltaF = ((1-alpha).^sigma*((1+muf).*pfstar).^(-sigma))./...
    (alpha.^sigma*((1+muh).*ph).^(-(sigma-1)) +...
    (1-alpha).^sigma*((1+muf).*pfstar).^(-(sigma-1))).^(sigma/(sigma-1));

deltaFstar = (alpha.^sigma*((1+mufstar).*pfstar).^(-sigma))./...
    (alpha.^sigma*((1+mufstar).*pfstar).^(-(sigma-1)) +...
    (1-alpha)^sigma*((1+muhstar).*ph).^(-(sigma-1))).^(sigma/(sigma-1));
deltaHstar = ((1-alpha).^sigma*((1+muhstar).*ph).^(-sigma))./...
    (alpha.^sigma*((1+mufstar).*pfstar).^(-(sigma-1)) +...
    (1-alpha)^sigma*((1+muhstar).*ph).^(-(sigma-1))).^(sigma/(sigma-1));


% labor of infected and recovered
chi1 = alpha*(deltaH).^(-1/sigma)./(kappa*ph.*(1+muh));
chi2 = (1+muh).*ph.*deltaH + (1+muf).*pfstar.*deltaF;
lr = (-g + ((g).^2 + 4*w.^2.*chi1.*chi2).^(1/2))./(2*w);
li = (-g + ((g).^2 + 4*(phi*w).^2.*chi1.*chi2).^(1/2))./(2*phi*w);

chi1star = alpha*(deltaFstar).^(-1/sigma)./(kappa*pfstar.*(1+mufstar));
chi2star = (1+mufstar).*pfstar.*deltaFstar + (1+muhstar).*ph.*deltaHstar;
lrstar = (-gstar + ((gstar).^2 + 4*wstar.^2.*chi1star.*chi2star).^(1/2))./(2*wstar);
listar = (-gstar + ((gstar).^2 + 4*(phi*wstar).^2.*chi1star.*chi2star).^(1/2))./(2*phi*wstar);

% consumption of infected and recovered
xi = (phi*w.*li + g)./chi2;
xistar = (phi*wstar.*listar + gstar)./chi2star;
xr = (w.*lr + g)./chi2;
xrstar = (wstar.*lrstar + gstar)./chi2star;

chi = deltaH.*xi;
chr = deltaH.*xr;
cfi = deltaF.*xi;
cfr = deltaF.*xr;

chistar = deltaHstar.*xistar;
chrstar = deltaHstar.*xrstar;
cfistar = deltaFstar.*xistar;
cfrstar = deltaFstar.*xrstar;


% pre-allocate
I=NaN*ones(HH,1);
S=NaN*ones(HH,1);
R=NaN*ones(HH,1);
pop=NaN*ones(HH,1);
T=NaN*ones(HH,1);

Istar=NaN*ones(HH,1);
Sstar=NaN*ones(HH,1);
Rstar=NaN*ones(HH,1);
popstar=NaN*ones(HH,1);
Tstar=NaN*ones(HH,1);

%initial conditions
I(1)=I0;
S(1)=pop0-I(1);
R(1)=0;
pop(1)=pop0;

Istar(1)=Istar0;
Sstar(1)=popstar0-Istar(1);
Rstar(1)=0;
popstar(1)=popstar0;

%iterate on SIR equations
for t=1:(HH-1)
    T(t)=(pi1*chs(t).*chi(t)+pi1*cfs(t).*cfi(t)+pi2*ls(t)*li(t)+pi3).*S(t).*I(t) ...
         + (pi4*chs(t).*chistar(t)+pi4*cfs(t).*cfistar(t)).*S(t).*Istar(t);
    
    Tstar(t)=(pi1*cfsstar(t).*cfistar(t)+pi1*chsstar(t).*chistar(t)+pi2*lsstar(t)*listar(t)+pi3)*Sstar(t)*Istar(t)...
             + (pi4*cfsstar(t).*cfi(t)+pi4*chsstar(t).*chi(t)).*Sstar(t).*I(t);
    
    S(t+1)=S(t)-T(t);
    Sstar(t+1)=Sstar(t)-Tstar(t);
    
    % I(t+1)=I(t) + T(t) - (pir + pid)*I(t);
    % Istar(t+1)=Istar(t) + Tstar(t) - (pir + pid)*Istar(t);    
    I(t+1)=I(t) + T(t) - (pir + pid + pid_accel*I(t))*I(t);   
    Istar(t+1)=Istar(t) + Tstar(t) - (pir + pid + pid_accel*Istar(t))*Istar(t);        
    R(t+1)=R(t)+pir*I(t);
    Rstar(t+1)=Rstar(t)+pir*Istar(t);

    pop(t+1) = pop(t) - (pid + pid_accel*I(t))*I(t);
    popstar(t+1) = popstar(t) - (pid + pid_accel*Istar(t))*Istar(t);
    % pop(t+1) = pop(t) - pid*I(t);
    % popstar(t+1) = popstar(t) - pid*Istar(t);
end

% vaccine arrives at t=HH, perfect foresight
R(HH) = R(HH) + I(HH);
I(HH) = 0;
Rstar(HH) = Rstar(HH) + Istar(HH);
Istar(HH) = 0;

% T(HH-1) = 0; % perfect foresight
% pop(HH) = pop(HH-1);
% S(HH) = S(HH-1);
% R(HH) = pop(HH)-S(HH);
% Tstar(HH-1) = 0;
% popstar(HH) = popstar(HH-1);
% Sstar(HH) = Sstar(HH-1);
% Rstar(HH) = popstar(HH)-Sstar(HH);

T(HH) = 0;
Tstar(HH) = 0;

%%%%%%%%%%%
% Utility %
%%%%%%%%%%%

% I should still be R due to labor productivity difference

%Recovered 
ur=log(xr)-kappa/2*lr.^2;
Ur=NaN*ones(HH+1,1);
Ur(HH+1)=Urss;
for tt=HH:-1:1
    Ur(tt)=ur(tt)+beta*Ur(tt+1);    
end


urstar=log(xrstar)-kappa/2*lrstar.^2;
Urstar=NaN*ones(HH+1,1);
Urstar(HH+1)=Urss;
for tt=HH:-1:1
    Urstar(tt)=urstar(tt)+beta*Urstar(tt+1);    
end


%Infected 
ui=log(xi)-kappa/2*li.^2;
Ui=NaN*ones(HH+1,1);
Ui(HH+1)=Urss; % this is Urss since at HH+1 they are cured
Ui(HH)=ur(HH) + beta*Urss; % this is Urss since at HH they are cured
for tt=HH-1:-1:1
    Ui(tt)=ui(tt)+beta*(pir*Ur(tt+1)+(pid + pid_accel*I(tt))*Udeath+...
                        (1-pir-(pid + pid_accel*I(tt)))*Ui(tt+1)); 
    % Ui(tt)=ui(tt)+beta*(pir*Ur(tt+1)+pid*Udeath+(1-pir-pid)*Ui(tt+1));    
end

uistar=log(xistar)-kappa/2*listar.^2;
Uistar=NaN*ones(HH+1,1);
Uistar(HH+1)=Urss; % this is Urss since at HH+1 they are cured
Uistar(HH)=urstar(HH) + beta*Urss; % this is Urss since at HH they are cured
for tt=HH-1:-1:1
    Uistar(tt)=uistar(tt)+beta*(pir*Urstar(tt+1)+(pid + pid_accel*Istar(tt))*Udeath+...
        (1-pir-(pid + pid_accel*Istar(tt)))*Uistar(tt+1));        
    % Uistar(tt)=uistar(tt)+beta*(pir*Urstar(tt+1)+pid*Udeath+(1-pir-pid)*Uistar(tt+1));    
end


%Susceptible
Usss=Urss; % PV utility of susceptibles same as recovered in steady state

us=log(xs)-kappa/2*ls.^2;
Us=NaN*ones(HH+1,1);
Us(HH+1)=Usss;
Us(HH)=us(HH) + beta*Urss;
for tt=HH-1:-1:1
    taus=T(tt)/S(tt);
    Us(tt)=us(tt)+beta*((1-taus)*Us(tt+1) + taus*Ui(tt+1));    
end

usstar=log(xsstar)-kappa/2*lsstar.^2;
Usstar=NaN*ones(HH+1,1);
Usstar(HH+1)=Usss;
Usstar(HH)=usstar(HH) + beta*Urss;
for tt=HH-1:-1:1    
    tausstar=Tstar(tt)/Sstar(tt);
    Usstar(tt)=usstar(tt)+beta*((1-tausstar)*Usstar(tt+1) + tausstar*Uistar(tt+1));
end


aggU = S(1)*Us(1) + I(1)*Ui(1);
aggUstar = Sstar(1)*Usstar(1) + Istar(1)*Uistar(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium condition residuals %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = S(1:HH);
Sstar = Sstar(1:HH);
I = I(1:HH);
Istar = Istar(1:HH);
R = R(1:HH);
Rstar = Rstar(1:HH);
pop = pop(1:HH);
popstar = popstar(1:HH);

% aggregate consumption
aggC = S.*xs + I.*xi + R.*xr;
aggCstar = Sstar.*xsstar + Istar.*xistar + Rstar.*xrstar;

% of each good
aggch = S.*chs + I.*chi + R.*chr;
aggchstar = Sstar.*chsstar + Istar.*chistar + Rstar.*chrstar;

aggcf = S.*cfs + I.*cfi + R.*cfr;
aggcfstar = Sstar.*cfsstar + Istar.*cfistar + Rstar.*cfrstar;

aggL=S.*ls + I.*phi.*li + R.*lr;
aggLstar=Sstar.*lsstar + Istar.*phi.*listar + Rstar.*lrstar;

err=NaN*ones(7*HH,1);

% Susceptible FOC
%h Eq. 30 and 31
err(1:HH) = ...
    w.*alpha.*xs.^(1/sigma-rho).*chs.^(-1/sigma) - kappa*(1+muh).*ph.*ls ...
    - beta.*(pi1.*w.*chi.*I + pi4.*w.*chistar.*Istar +...
    pi2.*(1+muh).*ph.*li.*I).*(Us(2:end)-Ui(2:end));

err(HH+1:2*HH) = ...
    w.*(1-alpha).*xs.^(1/sigma-rho).*cfs.^(-1/sigma) - kappa*(1+muf).*pfstar.*ls ...
    - beta.*(pi1.*w.*cfi.*I + pi4.*w.*cfistar.*Istar +...
    pi2.*(1+muf).*pfstar.*li.*I).*(Us(2:end)-Ui(2:end));

%f
err(2*HH+1:3*HH) = ...
    wstar.*alpha.*xsstar.^(1/sigma-rho).*cfsstar.^(-1/sigma) - kappa*(1+mufstar).*pfstar.*lsstar ...
    - beta.*(pi1.*wstar.*cfistar.*Istar + pi4.*wstar.*cfi.*I +...
    pi2.*(1+mufstar).*pfstar.*listar.*Istar).*(Usstar(2:end)-Uistar(2:end));

err(3*HH+1:4*HH) = ...
    wstar.*(1-alpha).*xsstar.^(1/sigma-rho).*chsstar.^(-1/sigma) - kappa*(1+muhstar).*ph.*lsstar ...
    - beta.*(pi1.*wstar.*chistar.*Istar + pi4.*wstar.*chi.*I +...
    pi2.*(1+muhstar).*ph.*listar.*Istar).*(Usstar(2:end)-Uistar(2:end));

% Market Clearing
err(4*HH+1:5*HH) = aggch + aggchstar + ...
    (1-delta_mu).*(muh.*aggch + mufstar.*aggchstar) - z*aggL;
tmp_check = aggcf + aggcfstar + ...
    (1-delta_mu).*(muh.*aggcf + mufstar.*aggcfstar) - zstar*aggLstar;

% Government budget
err(5*HH+1:6*HH) = (delta_mu.*muh.*(ph.*aggch + pfstar.*aggcf) + ...
                    delta_nu.*(muf-muh).*pfstar.*aggcf) - g.*pop;
                
err(6*HH+1:7*HH) = (delta_nu.*(muhstar-mufstar).*ph.*aggchstar + ...
                    delta_mu.*mufstar.*(ph.*aggchstar + pfstar.*aggcfstar)) - gstar.*popstar;

% err(5*HH+1:6*HH) = delta.*(muh.*ph.*aggch + ...
%          muf.*pfstar.*aggcf) - g.*pop;
% 
% err(6*HH+1:7*HH) = delta.*(muhstar.*ph.*aggchstar + ...
%          mufstar.*pfstar.*aggcfstar) - gstar.*popstar;


% scaling of errors
err(1:4*HH) = err(1:4*HH)*1e3;
% err(4*HH+1:5*HH) = err(4*HH+1:5*HH)*1e-4;
% err(5*HH+1:6*HH) = err(5*HH+1:6*HH)*1e-4;
% err(6*HH+1:7*HH) = err(6*HH+1:7*HH)*1e-4;

% mean(abs(err(1:4*HH)))
% mean(abs(err(4*HH+1:7*HH)))

end