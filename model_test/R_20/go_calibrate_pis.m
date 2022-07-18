%calibrate pis1, pis2 and pis3 using SIR model
scale1=1000000;scale2=1000;%scale pis for numerical solver
[sol,fval,exitflag]=fsolve(@calibrate_pis,[0.2;0.2;0.2],opts_fsolve,HH,i_ini,pop_ini,pir,pid*eta,pis1_shr_target,pis2_shr_target,RplusD_target,phii,crss,nrss,scale1,scale2,wah,beds,tau,eta);

if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else
    [err,pis1,pis2,pis3,RnotSIR,I,S,D,R,T] =calibrate_pis(sol,HH,i_ini,pop_ini,pir,pid*eta,pis1_shr_target,pis2_shr_target,RplusD_target,phii,crss,nrss,scale1,scale2,wah,beds,tau,eta);
    
    disp(['Max. abs. error in calibration targets:',num2str(max(abs(err)))]);
    disp([' ']);
    pis1=sol(1)/scale1;
    pis2=sol(2)/scale2;
    pis3=sol(3);
    RnotSIR;
end

clear err RnotSIR I S D R T
