function [err,I,S,R,D,T,Pop,cs,ns,Us,RnotSIRmacro,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U,pid_endo] = get_err(guess,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa,cbar,Loan,wah,beds,tau,eta)

%back out guesses for ns,ni,nr
ns=guess(1:HH);
ni=guess(HH+1:2*HH);
nr=guess(2*HH+1:3*HH); 

%equilibrium equations

%Recovered people
lambr=(theta*nr)./A;
cr=((1+muc).*lambr).^(-1)+cbar;
ur=log(cr-cbar)-theta/2*nr.^2;
%%%%%%%
% if min(cr)>cbar
%     ur=log(cr-cbar)-theta/2*nr.^2;
% else
%     ur = -1e20*ones(size(cr));
% end
%%%%%%
Ur=NaN*ones(HH+1,1);Ur(HH+1)=Urss;
for tt=HH:-1:1    
    Ur(tt,1)=ur(tt)+betta*Ur(tt+1,1);    
end
Gamma=(1+muc).*cr-A.*nr-Loan;

%Infected People
lambi=(theta*ni)./(phii*A);
ci=((1+muc).*lambi).^(-1)+cbar;
ui=log(ci-cbar)-theta/2*ni.^2;

%%%%%%%
% if min(ci)>cbar
%     ui=log(ci-cbar)-theta/2*ni.^2;
% else
%     ui = -1e20*ones(size(ci));
% end
%%%%%%%




%Susceptible People
cs=1./(1+muc).*(A.*ns+Gamma+Loan);
us=log(cs-cbar)-theta/2*ns.^2;

%%%%%%%
% if min(cs)>cbar
%     us=log(cs-cbar)-theta/2*ns.^2;
% else
%     us = -1e20*ones(size(cs));
% end
%%%%%%%

%pre-allocate
I=NaN*ones(HH+1,1);
S=NaN*ones(HH+1,1);
D=NaN*ones(HH+1,1);
R=NaN*ones(HH+1,1);
Pop=NaN*ones(HH+1,1);
T=NaN*ones(HH,1);

%initial conditions
Pop(1)=pop_ini;
I(1)=i_ini;
S(1)=Pop(1)-I(1);
D(1)=0;
R(1)=0;

%%Endogenous death probability
pid_endo=NaN*ones(HH,1);

%iterate on SIR equations
for j=1:1:HH
    T(j,1)=pis1*S(j)*cs(j)*I(j)*ci(j)/crss^2+pis2*S(j)*ns(j)*I(j)*ni(j)*(1-wah)+pis3*S(j)*I(j);
    T(j,1) = (1-eta)*T(j,1);
    pid_endo(j,1)=pid+kappa*I(j)^2;
    S(j+1,1)=S(j)-T(j);
    I(j+1,1)=I(j)+T(j)-(pir+pid_endo(j,1))*I(j);
    R(j+1,1)=R(j)+pir*I(j);
    if eta*I(j)<beds
        D(j+1,1)=D(j)+pid_endo(j,1)*I(j);
    else
        D(j+1,1)=D(j)+pid_endo(j,1)*beds+pid_endo(j,1)/eta*tau*(eta*I(j)-beds);
    end
    Pop(j+1,1)=1-D(j+1,1);
end

%Infected People (continued)
Ui=NaN*ones(HH+1,1); 
Ui(HH+1)=(1-deltac)^HH*Uiss+(1-(1-deltac)^HH)*Urss;%terminal condition 
%Ui(HH+1)=Uiss;%
%Ui(HH+1)=Urss;%
for tt=HH:-1:1    
    %Ui(tt,1)=ui(tt)+(1-deltac)*betta*((1-pir-pid_endo(tt))*Ui(tt+1,1)+pir*Ur(tt+1,1))+deltac*betta*Ur(tt+1,1);    
    Ui(tt,1)=(1-eta)*ui(tt)+(1-deltac)*betta*((1-pir-pid_endo(tt))*Ui(tt+1,1)+pir*Ur(tt+1,1))+deltac*betta*Ur(tt+1,1);    
end

%Susceptible People (continued)
Us=NaN*ones(HH+1,1);
Usss=Urss; %PV utility of susceptibles same as recovered in steady state
Us(HH+1)=(1-deltav)^HH*Usss+(1-(1-deltav)^HH)*Urss;%terminal condition
for tt=HH:-1:1    
    Us(tt,1)=us(tt)+(1-deltav)*betta*(1-T(tt)/S(tt)).*Us(tt+1,1)+(1-deltav)*betta*T(tt)/S(tt).*Ui(tt+1,1)+deltav*betta*Ur(tt+1,1);
end

%Lagrange multipliers susceptibles
lamtau=(1-deltav)*betta*(Ui(2:HH+1)-Us(2:HH+1));
lambs=((cs-cbar).^(-1)+lamtau*pis1.*I(1:HH).*ci)./(1+muc);

%equation residuals
err=NaN*ones(3*HH,1);
err(1:HH)=(1+muc).*ci-phii*A.*ni-Gamma-Loan; % FOC consumption infected
err(HH+1:2*HH)=muc.*(S(1:HH).*cs+I(1:HH).*ci+R(1:HH).*cr)-Gamma.*(S(1:HH)+R(1:HH)+I(1:HH)); % Gov bc
err(2*HH+1:3*HH)=-theta*ns+A.*lambs+lamtau*pis2.*I(1:HH).*ni; % S - FOC wrt n



%Aggregate consumption and hours
aggC=S(1:HH).*cs+I(1:HH).*ci+R(1:HH).*cr;
aggH=S(1:HH).*ns+I(1:HH).*ni*phii+R(1:HH).*nr;

%Present value of society utility 
%U=S(1:HH).*Us(1:HH)+I(1:HH).*Ui(1:HH)*(1-eta)+R(1:HH).*Ur(1:HH);
U=S(1:HH).*Us(1:HH)+I(1:HH).*Ui(1:HH)+R(1:HH).*Ur(1:HH);

RnotSIRmacro=T(1)/I(1)/(pir+pid);
end