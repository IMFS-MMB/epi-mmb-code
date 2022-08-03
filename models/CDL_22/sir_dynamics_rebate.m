% SIR dynamics and value functions of Susceptible, infected and recovered
% in model with tax and rebate by government
function [pol_mkt_s_opt,wage_s_opt,RnotSIR,I,Ih,Im,S,Sh,Sm,R,D,V_r,V_i,V_sh,V_sm,Cs,Ci,Cr,C,Gi,Y,theta_total] = sir_dynamics_rebate(pol_mkt_s,...
    wage_s,sec_weight,sec_elas,occ_weight,occ_elas,tau,below_med,above_med,wage_sm,wage_sh,wage_r,wage_ih,wage_im,wage_i,pol_mkt_r,pol_mkt_i,...
    risk_prop,ini_pop,betta,pi_par,periods)

n_occ = size(occ_weight,2); % number of occupations
n_sec = size(sec_weight,1); % number of sectors

% infection parameters
pir = pi_par(1); %Weekly probability of recovering
pi = pi_par(2); %random infection rate between susceptible and infected
varepsilon = pi_par(3);%infection rate between susceptible and infected in market work
epsilon = pi_par(4);  %initial infected
pic = pi_par(5); % infection rate between consumption of susceptible and infected
pid = pi_par(6); % Weekly probability of dying

mkt_risk = varepsilon.*(1+risk_prop);

%pre-allocate
Cs = NaN*ones(periods+1,n_occ);
Ci = NaN*ones(periods+1,n_occ);
Cr = NaN*ones(periods+1,n_occ);
C = NaN*ones(periods+1,1);
Gi = NaN*ones(periods+1,n_occ);
Y = NaN*ones(periods+1,1);
y_sec = NaN*ones(n_sec,1);
y_occ = NaN*ones(n_sec,n_occ); % output at time t, in sector j in occupation i

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
    
   Gi(t,below_med) = sum(tau.*( wage_s(t,:).*Si(t,:) +  wage_r(t,:).*Ri(t,:) ...
        + wage_i(t,:).*Ii(t,:)),2)./sum(Si(t,below_med) + Ri(t,below_med) + Ii(t,below_med),2);
    
    Gi(t,[1:4,above_med]) = 0;
        
    Cs(t,:) = (1-tau).*wage_s(t,:).*Si(t,:) + Gi(t,:).*Si(t,:);
                        
    Ci(t,:) = (1-tau).*wage_i(t,:).*Ii(t,:) + Gi(t,:).*Ii(t,:);
    
    Cr(t,:) = (1-tau).*wage_r(t,:).*Ri(t,:) + Gi(t,:).*Ri(t,:);
     
    C(t,1) = sum(Cs(t,:) + Ci(t,:) + Cr(t,:));
    
    for j = 1:n_sec
        y_occ(j,:) = wage_r(t,:).*Ri(t,:) + wage_s(t,:).*Si(t,:) + wage_i(t,:).*Ii(t,:);
        y_sec(j) = sum(occ_weight(j,:).*(y_occ(j,:).^(occ_elas(j))) ).^(1/occ_elas(j));
    end

    Y(t,1) = sum(sec_weight.*(y_sec.^sec_elas)).^(1/sec_elas);
    
    Ti(t,:) = pic.*Cs(t,:).*sum(Ci(t,:)) + pi.*Si(t,:).*I(t,1) + mkt_risk.*Si(t,:).*Im(t,1).*(pol_mkt_s(t,:)==2);
       
    Si(t+1,:) = Si(t,:) - Ti(t,:);
    Ii(t+1,:) = (1-pir-pid)*Ii(t,:) + Ti(t,:);
    Ri(t+1,:) = Ri(t,:) + pir*Ii(t,:);
    Pi(t+1,:) = Pi(t,:) - pid*Ii(t,:);
    Di(t+1,:) = Di(t,:) + pid*Ii(t,:);
       
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

% Last period
%Gi(end,below_med) = sum(tau.*( wage_s(end,:).*Si(end,:) +  wage_r(end,:).*Ri(end,:) ...
%        + wage_i(end,:).*Ii(end,:)),2)./sum(Si(end,below_med) + Ri(end,below_med) + Ii(end,below_med),2);
    
%Gi(end,[1:4,above_med]) = 0;
Gi(end,:) = 0; % last period return to baseline no benefits scenario
        
Cs(end,:) = (1-tau).*wage_s(end,:).*Si(end,:) + Gi(end,:).*Si(end,:);
                        
Ci(end,:) = (1-tau).*wage_i(end,:).*Ii(end,:) + Gi(end,:).*Ii(end,:);
    
Cr(end,:) = (1-tau).*wage_r(end,:).*Ri(end,:) + Gi(end,:).*Ri(end,:);

     
C(end,1) = sum(Cs(end,:) + Ci(end,:) + Cr(end,:));
    
for j = 1:n_sec
    y_occ(j,:) = wage_r(end,:).*Ri(end,:) + wage_s(end,:).*Si(end,:) + wage_i(end,1).*Ii(end,:);
    y_sec(j) = sum(occ_weight(j,:).*(y_occ(j,:).^(occ_elas(j))) ).^(1/occ_elas(j));
end

Y(end,1) = sum(sec_weight.*(y_sec.^sec_elas)).^(1/sec_elas);

theta_total =  Sm + Rm+ Im;

pim = repmat(mkt_risk,periods,1).*repmat(Im(1:end-1,1),1,n_occ);
pi_t = pi*I(1:end-1,1) + pic.*sum(Ci(1:end-1,:),2);

%pre-allocate
V_r = NaN*ones(periods+1,n_occ); % Value function of recovered
V_i = NaN*ones(periods+1,n_occ);% Value function of infected

V_r(end,:) = max([log( (1-tau).*wage_sh(end,:) + Gi(end,:) ) ;log( (1-tau).*wage_sm(end,:) + Gi(end,:) )],[],1)/(1-betta); % last period value of recovered
V_i(end,:) = max([log( (1-tau).*wage_ih(end,:) + Gi(end,:) ) ;log( (1-tau).*wage_im(end,:) + Gi(end,:) )],[],1)/(1-betta*(1-pir-pid)) + V_r(end,:)/(1-betta*pir);% last period value of infected

for t = periods:-1:1
   V_r(t,:) = max([log( (1-tau).*wage_sh(t,:) + Gi(t,:) );log( (1-tau).*wage_sm(t,:)  + Gi(t,:) )],[],1) + betta*V_r(t+1,:); % max between wage of market vs. home + next period value
   V_i(t,:) = max([log( (1-tau).*wage_ih(t,:) + Gi(t,:) );log( (1-tau).*wage_im(t,:)  + Gi(t,:) )],[],1) + betta*pir*V_r(t+1,:) + betta*(1-pir-pid)*V_i(t+1,:); % max between wage of market vs. home + next period value
end

%pre-allocate
V_sh = NaN*ones(periods+1,n_occ); % Value function of susceptible at home
V_sm = NaN*ones(periods+1,n_occ); % Value function of susceptible at market
pol_mkt_s_opt = NaN*ones(periods+1,n_occ); % optimal choice of susceptible individuals work at market or home
wage_s_opt = NaN*ones(periods+1,n_occ); % optimal wage of susceptibel individuals

%last period
V_sh(end,:) = log( (1-tau).*wage_sh(end,:) + Gi(end,:) )/(1-betta);
V_sm(end,:) = log( (1-tau).*wage_sm(end,:) + Gi(end,:) )/(1-betta);

for i = 1:n_occ
    [~,pol_mkt_s_opt(end,i)] = max([V_sh(end,i); V_sm(end,i)]);
    wage_s_opt(end,i) = wage_sh(end,i).*(pol_mkt_s_opt(end,i)==1) + wage_sm(end,i).*(pol_mkt_s_opt(end,i)==2);
end

for t = periods:-1:1
    for i = 1:n_occ
        V_sh(t,i) = log( (1-tau).*wage_sh(t,i) + Gi(t,i) ) + betta*(1-pi_t(t,1))*max([V_sh(t+1,i),V_sm(t+1,i)]) ...
            + betta*pi_t(t,1)*V_i(t+1,i);
        V_sm(t,i) = log( (1-tau).*wage_sm(t,i) + Gi(t,i) ) + betta*(pim(t,i)+pi_t(t,1))*V_i(t+1,i) ...
            + betta*(1-pim(t,i)-pi_t(t,1))*max([V_sh(t+1,i),V_sm(t+1,i)]);
        
        [~,pol_mkt_s_opt(t,i)] = max([V_sh(t,i); V_sm(t,i)]);
         wage_s_opt(t,i) = wage_sh(t,i).*(pol_mkt_s_opt(t,i)==1) + wage_sm(t,i).*(pol_mkt_s_opt(t,i)==2);
    end       
end

RnotSIR=T(1)/I(1)/(pir+pid);
end