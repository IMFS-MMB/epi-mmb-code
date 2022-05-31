function [err,pisy1,pisy2,pisy3,pisy4,piso1,piso2,piso3,Iy,Io,Sy,So,Ry,Ro,Dy,Do,Ty,To,Ny,No,N] =calibrate_pis(...
    pis_guess,FF, cry_ss,nry_ss,cro_ss,Zyy,Zyo,Zoy,Zoo,scale1,scale2,...
    Iy1,Io1,Sy1,So1,fy,fo,N1, ...
    piy_d,pio_d,piy_r,pio_r,eta, gamma, kappa, other_share)

pisy1=pis_guess(1)/scale1;
pisy2=pis_guess(2)/scale1;
pisy3=pis_guess(3)/scale1;
pisy4=pis_guess(4)/scale2;
piso1=pis_guess(5)/scale1;
piso2=pis_guess(6)/scale1;
piso3=pis_guess(7)/scale2;


%pre-allocate and initial conditions
Iy= [Iy1; zeros(FF,1)];
Io= [Io1; zeros(FF,1)];
Sy= [Sy1; zeros(FF,1)];
So= [So1; zeros(FF,1)];
Ty= zeros(FF,1);
To= zeros(FF,1);
Ry= zeros(FF+1,1);
Ro= zeros(FF+1,1);
Dy= zeros(FF+1,1);
Do= zeros(FF+1,1);
No= [fo; zeros(FF,1)];
Ny= [fy; zeros(FF,1)];
N= [N1; zeros(FF,1)];

%%Endogenous death probability
piy_d_endo=NaN*ones(FF,1);
pio_d_endo=NaN*ones(FF,1);

%iterate on SIR equations
for i=1:FF
    Ty(i)=eta*Sy(i)*(...
        pisy1*Zyy*(Iy(i)/fy)*cry_ss^2+...
        pisy2*Zyo*(Io(i)/fo)*cry_ss*cro_ss+...
        pisy3*Zyy*(Iy(i)/fy)*nry_ss^2+...
        pisy4*(Zyy*(Iy(i)/fy)+Zyo*(Io(i)/fo)));
    To(i)=eta*So(i)*(...
        piso1*Zoo*(Io(i)/fo)*cro_ss^2+...
        piso2*Zoy*(Iy(i)/fy)*cro_ss*cry_ss+...
        piso3*(Zoy*(Iy(i)/fy)+Zoo*(Io(i)/fo))); 
    Sy(i+1)=Sy(i)-Ty(i);
    So(i+1)=So(i)-To(i);
    piy_d_endo(i,1)=piy_d+kappa*Iy(i)^2;
    pio_d_endo(i,1)=pio_d+kappa*Io(i)^2;
    Iy(i+1)=Iy(i)+Ty(i)-(piy_d_endo(i,1)+piy_r)*Iy(i);
    Io(i+1)=Io(i)+To(i)-(pio_d_endo(i,1)+pio_r)*Io(i);
    Ry(i+1)=Ry(i)+piy_r*Iy(i);
    Ro(i+1)=Ro(i)+pio_r*Io(i);
    Dy(i+1)=Dy(i)+piy_d_endo(i,1)*Iy(i);
    Do(i+1)=Do(i)+pio_d_endo(i,1)*Io(i);
    Ny(i+1)=Sy(i+1)+Iy(i+1)+Ry(i+1);
    No(i+1)=So(i+1)+Io(i+1)+Ro(i+1);
    N(i+1)=Ny(i+1)+No(i+1);
    % check Pop = Pop1
%     Pop1(i+1)=Pop(i)-(piy_d*Iy(i)+pio_d*Io(i));
end

Piy=pisy1*Zyy*cry_ss^2+...
    pisy2*Zyo*cry_ss*cro_ss+...
    pisy3*Zyy*nry_ss^2+...
    pisy4*(Zyy+Zyo);

Pio=piso1*Zoo*cro_ss^2+...
    piso2*Zoy*cry_ss*cro_ss+...
    piso3*(Zoy+Zoo);

% other_share=2/3;
% other_share=1;
% load('C:\Users\papet\Dropbox\SIRage_macro\model\canonicalSIRage_Towers\RplusD_Towers')
load RplusD_Towers.mat
RplusD_y_target=RplusD_y_Towers;
RplusD_o_target=RplusD_o_Towers;

err(1)=pisy1*Zyy*cry_ss^2/Piy-(1-other_share)/3;
err(2)=pisy2*Zyo*cry_ss*cro_ss/Piy-(1-other_share)/3;
err(3)=pisy3*Zyy*nry_ss^2/Piy-(1-other_share)/3;
err(4)=RplusD_y_target-(Ry(end)+Dy(end));
% err(4)=RplusD_target-(Ry(end)+Ro(end)+Do(end)+Dy(end));
% err(4)=pisy4*(Zyy+Zyo)/Piy-other_share;
err(5)=piso1*Zoo*cro_ss^2/Pio-(1-other_share)/2;
err(6)=piso2*Zoy*cry_ss*cro_ss/Pio-(1-other_share)/2;
% err(6)=piso3*(Zoy+Zoo)/Pio-other_share;
% err(7)= RplusD_target-(Dy(end)*gamma/piy_d+Do(end)*gamma/pio_d);
err(7)=RplusD_o_target-(Ro(end)+Do(end));
end
