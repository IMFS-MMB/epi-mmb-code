%calibrate steeady state

%guess theta, Ah, cbar, csl1 and P1
guess=[0.001343458194544   0.008733050233044   0.029289912030052   0.141323015469625   0.288799541242928];

%new guess
guess=[  0.001343458194544   0.008733050233044   0.029289912030052   0.111615051180050   87];%sh=0.1

scalee=10000; %scaling factor

fsolve_opts=optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1000000,'MaxIterations',100000);

%find parameters to match desired targets
sol=fsolve(@get_calib_err,guess,fsolve_opts,n_target,inc_target,winc_share_target,sh_ss,betta,P2_ss_numeraire_normalization,eta,bhstar_ss,pidh,pidl,mh,ml,pirh,pirl,scalee);

% theta_calib=sol(1);
% Ah_calib=sol(2)*scalee;
% cbar_calib=sol(3)*scalee;
% Al_calib=Ah_calib;
% csl1_calib=sol(4)*scalee;
% P1_calib=sol(5);

%Solve for Al to maintain same P1 as in baseline 
theta_calib=sol(1);
Ah_calib=sol(2)*scalee;
cbar_calib=sol(3)*scalee;
csl1_calib=sol(4)*scalee;
Al_calib=sol(5);

P1_calib=0.288799541242928; %baseline P1


[err]=get_calib_err(sol,n_target,inc_target,winc_share_target,sh_ss,betta,P2_ss_numeraire_normalization,eta,bhstar_ss,pidh,pidl,mh,ml,pirh,pirl,scalee);


function [err]=get_calib_err(guess,n_target,inc_target,winc_share_target,sh_ss,betta,P2_ss_numeraire_normalization,eta,bhstar_ss,pidh,pidl,mh,ml,pirh,pirl,scalee)
%back out guess
%     theta=abs(guess(1));
%     Ah=abs(guess(2))*scalee;
%     cbar=abs(guess(3))*scalee;
%     csl1_ss=abs(guess(4))*scalee;
%     P1_ss=abs(guess(5));
%
%     Al=Ah;

theta=abs(guess(1));
Ah=abs(guess(2))*scalee;
cbar=abs(guess(3))*scalee;
csl1_ss=abs(guess(4))*scalee;
Al=abs(guess(5));

P1_ss=0.288799541242928;%baseline P1


%steady state equations
sl_ss=1-sh_ss;
rstar=1/betta-1;
P2_ss=P2_ss_numeraire_normalization;

nsl_ss=Al/theta*(1-eta)/(csl1_ss-cbar);
csh1_ss=1/sh_ss*(Al*sl_ss*nsl_ss-sl_ss*csl1_ss);
nsh_ss=Ah/theta*(1-eta)/(csh1_ss-cbar)*P2_ss/P1_ss;

wl_ss=P1_ss*Al;
wh_ss=P2_ss*Ah;

csl2_ss=eta/(P2_ss/P1_ss*(1-eta)/(csl1_ss-cbar));
csh2_ss=eta/(P2_ss/P1_ss*(1-eta)/(csh1_ss-cbar));

C_ss=P1_ss*(sh_ss*csh1_ss+sl_ss*csl1_ss)+P2_ss*(sh_ss*csh2_ss+sl_ss*csl2_ss);
N_ss=sl_ss*nsl_ss+sh_ss*nsh_ss;

inch_ss=wh_ss*(sh_ss*nsh_ss)+rstar*P2_ss*bhstar_ss;
incl_ss=wl_ss*(sl_ss*nsl_ss);

%targets
err(1,1)=theta-0.001343458194544;%n_target-N_ss;
err(2,1)=Ah-87.330502330435110;%inc_target-(incl_ss+inch_ss);
err(3,1)=cbar-292.8991203005227;%winc_share_target-(sh_ss*wh_ss*nsh_ss)/(sh_ss*wh_ss*nsh_ss+sl_ss*wl_ss*nsl_ss);
err(4,1)=csl1_ss+eta/(1-eta)*(csl1_ss-cbar)-Al*nsl_ss;
err(5,1)=eta/(1-eta)*(csh1_ss-cbar)-P2_ss*Ah/sh_ss/P1_ss*(sh_ss*nsh_ss)+csh1_ss-rstar*P2_ss*bhstar_ss/sh_ss/P1_ss;
end
