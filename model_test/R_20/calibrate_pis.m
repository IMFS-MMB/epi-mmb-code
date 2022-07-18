function [err,pis1,pis2,pis3,RnotSIR,I,S,D,R,T] =calibrate_pis(pis_guess,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,C,N,scale1,scale2,wah,beds,tau,eta)

%back out initial guesses
pis1=pis_guess(1)/scale1;
pis2=pis_guess(2)/scale2;
pis3=pis_guess(3);

%pre-allocate
I=NaN*ones(HH+1,1);
S=NaN*ones(HH+1,1);
D=NaN*ones(HH+1,1);
R=NaN*ones(HH+1,1);
T=NaN*ones(HH,1);

%initial conditions
I(1)=i_ini;
S(1)=pop_ini-I(1);
D(1)=0;
R(1)=0;

%iterate on SIR model equations
for j=1:1:HH
    T(j,1)=pis1*S(j)*I(j)+pis2*S(j)*N^2*I(j)*(1-wah)+pis3*S(j)*I(j);
    T(j,1)=(1-eta)*T(j,1);
    S(j+1,1)=S(j)-T(j);
    I(j+1,1)=I(j)+T(j)-(pir+pid)*I(j);
    R(j+1,1)=R(j)+pir*I(j);
    if I(j)<beds
        D(j+1,1)=D(j)+pid*I(j);
    else
        D(j+1,1)=D(j)+pid*beds+pid*tau*(I(j)-beds);
    end
        
end

%calculate equation residuals for target equations
err(1)=pis1_shr_target-(pis1*C^2)/(pis1*C^2+pis2*N^2+pis3);
err(2)=pis2_shr_target-(pis2*N^2)/(pis1*C^2+pis2*N^2+pis3);
err(3)=RplusD_target-(R(end)+D(end));


RnotSIR=T(1)/I(1)/(pir+pid);

end