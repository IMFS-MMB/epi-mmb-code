pis_vec_guess=.2*ones(7,1);
FF=250; 
scale1=1000000;scale2=1000;%scale pis for numerical solver
opts_fsolve_pis=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve

Io1=eps*fo;
Iy1=eps-Io1;
So1=(1-eps)*fo;
Sy1=1-eps-So1;
N1=1;

[sol,fval,exitflag]=fsolve(@calibrate_pis,pis_vec_guess,opts_fsolve_pis,...
    FF, cry_ss,nry_ss,cro_ss,Zyy,Zyo,Zoy,Zoo,scale1,scale2,...
    Iy1,Io1,Sy1,So1,fy,fo,N1, ...
    piy_d, pio_d,piy_r,pio_r,eta, gamma, kappa, other_share);
if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else

[err,pisy1,pisy2,pisy3,pisy4,piso1,piso2,piso3,Iy,Io,Sy,So,Ry,Ro,Dy,Do,Ty,To,Ny,No,N]=calibrate_pis(...
    sol, FF,cry_ss,nry_ss,cro_ss,Zyy,Zyo,Zoy,Zoo,scale1,scale2,...
    Iy1,Io1,Sy1,So1,fy,fo,N1, ...
    piy_d, pio_d,piy_r,pio_r,eta, gamma, kappa, other_share);
pisy1=sol(1)/scale1;
pisy2=sol(2)/scale1;
pisy3=sol(3)/scale1;
pisy4=sol(4)/scale2;
piso1=sol(5)/scale1;
piso2=sol(6)/scale1;
piso3=sol(7)/scale2;
end
pisy1_base=pisy1;
pisy2_base=pisy2;
pisy3_base=pisy3;
pisy4_base=pisy4;
piso1_base=piso1;
piso2_base=piso2;
piso3_base=piso3;

%{
ia=2;
ib=3;
figure;
subplot(ia,ib,1);
plot(Iy), hold on, plot(Io), hold on, plot(Iy+Io)
legend('young', 'old', 'total')
title('I');
subplot(ia,ib,2);
plot(Sy), hold on, plot(So), hold on, plot(Sy+So)
title('S');
subplot(ia,ib,3);
plot(Ry), hold on, plot(Ro), hold on, plot(Ry+Ro)
title('R');
subplot(ia,ib,4);
plot(Dy), hold on, plot(Do), hold on, plot(Dy+Do)
title('D');
subplot(ia,ib,5);
plot(Ty), hold on, plot(To), hold on, plot(Ty+To)
title('T');
suptitle('SIR Model, Calibration')
%}
