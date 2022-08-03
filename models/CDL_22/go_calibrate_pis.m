%calibrate pis1, pis2 and pis3 using SIR model
scale1=1000000;scale2=1;%scale pis for numerical solver
% Scale 1 corresponds to consumption and will change if level of income
% changes


 [sol,fval,exitflag]=fsolve(@calibration_pi,[0.12;0.15;0.41],opts_fsolve,periods,theta_r,pol_mkt_s,pol_mkt_r,pol_mkt_i ...
    ,C_ss,n_occ,tau_m,wage_s,wage_r,wage_i,risk_prop,epsilon,pir,pid,ini_pop,scale1,scale2,pis1_shr_target,pis2_shr_target,RplusD_target);

if exitflag~=1
error('Fsolve could not calibrate the SIR model');

else
    [err,pi,varepsilon,pic,RnotSIR,~,~,~,R,~] =calibration_pi(sol,periods,theta_r,pol_mkt_s,pol_mkt_r,pol_mkt_i ...
    ,C_ss,n_occ,tau_m,wage_s,wage_r,wage_i,risk_prop,epsilon,pir,pid,ini_pop,scale1,scale2,pis1_shr_target,pis2_shr_target,RplusD_target);
    
    %disp(['Max. abs. error in calibration targets:',num2str(max(abs(err)))]);
    %disp([' ']);
    pi=sol(3);
    varepsilon=sol(2)/scale2;
    pic=sol(1)/scale1;
    RnotSIR;
end
