%calibrate pi1, pi2 and pi3 using SIR model
%nonlinear solver settings
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve
HH=200;

scale1=1000000;scale2=1000;%scale pi's for numerical solver
[sol,fval,exitflag]=fsolve(@calibrate,[0.17;0.27;0.85],opts_fsolve,HH,i_ini,pir,pid,pi1_shr_target,pi2_shr_target,RplusD_target,inc_target,n_target,scale1,scale2);

if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else
    [errSIR,pi1,pi2,pi3,RnotSIR] =calibrate(sol,HH,i_ini,pir,pid,pi1_shr_target,pi2_shr_target,RplusD_target,inc_target,n_target,scale1,scale2);
    
    disp(['Max. abs. error in calibration targets:',num2str(max(abs(errSIR)))]);
    disp([' ']);
    pi1=sol(1)/scale1
    pi2=sol(2)/scale2
    pi3=sol(3)
    RnotSIR
end