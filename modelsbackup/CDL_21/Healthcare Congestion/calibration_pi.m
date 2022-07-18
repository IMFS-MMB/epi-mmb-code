% Infection Dynamics needed to calibrate the infection probabilities in
% model with congestion in health care
% The new parameter kappa is introduced leading to a new calibration
function [err,pic,varepsilon,pi,RnotSIR,I,S,D,R,T] = calibration_pi(pis_guess,periods,theta_r,pol_mkt_s,pol_mkt_r,pol_mkt_i ...
    ,C_ss,n_occ,tau_m,wage_s,wage_r,wage_i,risk_prop,epsilon,pir,pid,kappa,ini_pop,scale1,scale2,pis1_shr_target,pis2_shr_target,RplusD_target)

%back out initial guesses
pi=pis_guess(3);
varepsilon=pis_guess(2)/scale2;
pic=pis_guess(1)/scale1;

mkt_risk = varepsilon.*(1+risk_prop);

%pre-allocate
Cs = NaN*ones(periods,n_occ);
Ci = NaN*ones(periods,n_occ);
Cr = NaN*ones(periods,n_occ);
C = NaN*ones(periods,1);
G = NaN*ones(periods,1);
I=NaN*ones(periods+1,1);
S=NaN*ones(periods+1,1);
R=NaN*ones(periods+1,1);
P = NaN*ones(periods+1,1);
D = NaN*ones(periods+1,1);
T = NaN*ones(periods,1);

Si = NaN*ones(periods+1,n_occ);
Ii = NaN*ones(periods+1,n_occ);
Ri = NaN*ones(periods+1,n_occ);
Ti = NaN*ones(periods,n_occ);


Sh=NaN*ones(periods+1,1);
Sm=NaN*ones(periods+1,1);
Ih=NaN*ones(periods+1,1);
Im=NaN*ones(periods+1,1);
Rh=NaN*ones(periods+1,1);
Rm=NaN*ones(periods+1,1);

 %number of infected each period after pandemic starts

% Time period 0 at start of pandemic
P(1,1) = 1; %initial population
R(1,1)=0;        %initial recovered
D(1,1) = 0; % deaths
Pi(1,:) = ini_pop;
Di(1,:) = zeros(1,n_occ);
Ri(1,:) = zeros(1,n_occ);
Rm(1,1) = 0; % recovered working at market
Rh(1,1) = 0; % recovered working at home

I(1,1)=epsilon;    %initial infected
Ii(1,:) = epsilon.*ini_pop;
Ih(1,1) = sum(Ii(1,:).*(pol_mkt_i(1,:)==1));% Infected working at home
Im(1,1) = sum(Ii(1,:).*(pol_mkt_i(1,:)==2));  % Infected working at market
S(1,1)=P(1,1)-epsilon;        %initial susceptible population
Si(1,:) = (P(1,1)-epsilon).*ini_pop;
Sh(1,1) = sum(Si(1,:).*(pol_mkt_s(1,:)==1)); % susceptible working at home
Sm(1,1) = sum(Si(1,:).*(pol_mkt_s(1,:)==2)); % susceptible working at market

for t = 1:1:periods
    
    G(t,1) = sum( tau_m(t,:).*( wage_s(t,:).*Si(t,:) +  wage_r(t,:).*Ri(t,:) + wage_i(t,:).*Ii(t,:)) );
    
    Cs(t,:) = (1-tau_m(t,:)).*wage_s(t,:).*Si(t,:) + G(t,1).*Si(t,:);
                        
    Ci(t,:) = (1-tau_m(t,:)).*wage_i(t,:).*Ii(t,:) + G(t,1).*Ii(t,:);
    
    Cr(t,:) = (1-tau_m(t,:)).*wage_r(t,:).*Ri(t,:) + G(t,1).*Ri(t,:);
     
    C(t,1) = sum(Cs(t,:) + Ci(t,:) + Cr(t,:));
    
    Ti(t,:) = pic.*Cs(t,:).*sum(Ci(t,:)) + pi.*Si(t,:).*I(t,1) + mkt_risk.*Si(t,:).*Im(t,1).*(pol_mkt_s(t,:)==2);
       
    Si(t+1,:) = Si(t,:) - Ti(t,:);
    pid_endo = pid+ kappa.*(I(t,1).^2);
    Ii(t+1,:) = (1-pir-pid_endo)*Ii(t,:) + Ti(t,:);
    Ri(t+1,:) = Ri(t,:) + pir*Ii(t,:);
    Pi(t+1,:) = Pi(t,:) - pid_endo*Ii(t,:);
    Di(t+1,:) = Di(t,:) + pid_endo*Ii(t,:);
       
    % Aggregate
    T(t,1) = sum(Ti(t,:));
    R(t+1,1) = sum(Ri(t+1,:));
    P(t+1,1) = sum(Pi(t+1,:));
    D(t+1,1) = sum(Di(t+1,:));
    I(t+1,1) = sum(Ii(t+1,:));
    S(t+1,1) = sum(Si(t+1,:));
    
    Ih(t+1,1) = sum(Ii(t+1,:).*(pol_mkt_i(t+1,:)==1));
    Im(t+1,1) = sum(Ii(t+1,:).*(pol_mkt_i(t+1,:)==2));
    Sh(t+1,1) = sum(Si(t+1,:).*(pol_mkt_s(t+1,:)==1));
    Sm(t+1,1) = sum(Si(t+1,:).*(pol_mkt_s(t+1,:)==2));
    Rh(t+1,1) = sum(Ri(t+1,:).*(pol_mkt_r(t+1,:)==1));
    Rm(t+1,1) = sum(Ri(t+1,:).*(pol_mkt_r(t+1,:)==2));
        
end

%calculate equation residuals for target equations
err(1)=pis1_shr_target-(pic*(C_ss(1,1).^2))./( (pic*(C_ss(1,1).^2)) + pi + varepsilon.*(theta_r(1,1).^2));

err(2)=pis2_shr_target- (varepsilon.*(theta_r(1,1).^2))./( (pic*(C_ss(1,1).^2)) + pi + varepsilon.*(theta_r(1,1).^2));

err(3)=RplusD_target-(R(end)+D(end));

%err(3) = ( T(1)/I(1)/(pir+pid)) - 2;

RnotSIR=T(1)/I(1)/(pir+pid);

end