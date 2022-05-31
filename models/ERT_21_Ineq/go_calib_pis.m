%calibrate pi's using standard SIR model

%nonlinear solver settings
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve
HH=100;

scale1=1000000;scale2=1000;%scale pi's for numerical solver
[sol,fval,exitflag]=fsolve(@calib_pis,[0.17;0.17;0.27;0.27;0.85],opts_fsolve,HH,i_ini,pirl,pirh,pidl,pidh,pi12_shr_target,pi3_shr_target,RplusD_target,scale1,scale2,sh_ss,sl_ss,csl1_ss,csl2_ss,csh1_ss,csh2_ss,nsl_ss,nsh_ss,cil1_ss,cil2_ss,cih1_ss,cih2_ss,nil_ss,nih_ss,pi1vs2factor,pi3lvshfactor);

[errSIR,pi1,pi2,pih3,pil3,pi4,ihSIR,ilSIR,shSIR,slSIR,dhSIR,dlSIR,rhSIR,rlSIR,taulSIR,tauhSIR] =calib_pis(sol,HH,i_ini,pirl,pirh,pidl,pidh,pi12_shr_target,pi3_shr_target,RplusD_target,scale1,scale2,sh_ss,sl_ss,csl1_ss,csl2_ss,csh1_ss,csh2_ss,nsl_ss,nsh_ss,cil1_ss,cil2_ss,cih1_ss,cih2_ss,nil_ss,nih_ss,pi1vs2factor,pi3lvshfactor);

disp(['Max. abs. error in calibration targets:',num2str(max(abs(errSIR)))]);
disp([' ']);
pi1=sol(1)/scale1
pi2=sol(2)/scale1
pih3=sol(3)/scale2
pil3=sol(4)/scale2
pi4=sol(5)
    