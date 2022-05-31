%calibrate pi1, pi2 and pi3 using SIR model
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-8,'MaxFunctionEvaluations',1000000,'MaxIterations',40000); %options for fsolve
HH=520;
scale1=1000000;scale2=1000;%scale pi for numerical solver

[sol,fval,exitflag]=fsolve(@calibrate_pi,[0.02;0.02;0.02],opts_fsolve,HH,i_ini,pir,pid,pi1_shr_target,pi2_shr_target,RplusD_target,C1bar,Ldbar,lsbar,scale1,scale2);

if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else
    [err,pi1,pi2,pi3,RnotSIR] =calibrate_pi(sol,HH,i_ini,pir,pid,pi1_shr_target,pi2_shr_target,RplusD_target,C1bar,Ldbar,lsbar,scale1,scale2);
    
    disp(['Max. abs. error in calibration targets:',num2str(max(abs(err)))]);
    disp([' ']);
    pi1=sol(1)/scale1
    pi2=sol(2)/scale2
    pi3=sol(3)
    RnotSIR
end 