function [err,pi1,pi2,pih3,pil3,pi4,ih,il,sh,sl,dh,dl,rh,rl,taul,tauh] =calib_pis(pi_guess,HH,i_ini,pirl,pirh,pidl,pidh,pi12_shr_target,pi3_shr_target,RplusD_target,scale1,scale2,sh_ss,sl_ss,csl1,csl2,csh1,csh2,nsl,nsh,cil1,cil2,cih1,cih2,nil,nih,pi1vs2factor,pi3lvshfactor)


%back out initial guesses
pi1=pi_guess(1)/scale1;
pi2=pi_guess(2)/scale1;
pih3=pi_guess(3)/scale2;
pil3=pi_guess(4)/scale2;
pi4=pi_guess(5);

%pre-allocate
ih=NaN*ones(HH+1,1);
il=NaN*ones(HH+1,1);
sh=NaN*ones(HH+1,1);
sl=NaN*ones(HH+1,1);
dh=NaN*ones(HH+1,1);
dl=NaN*ones(HH+1,1);
rh=NaN*ones(HH+1,1);
rl=NaN*ones(HH+1,1);
taul=NaN*ones(HH,1);
tauh=NaN*ones(HH,1);

if i_ini==0
    i_ini=0.001;
    disp(' ');
    disp('Setting i_ini=0.001');
    disp(' ');
end

%initial conditions
ih(1)=sh_ss*i_ini;
il(1)=(1-sh_ss)*i_ini;
sh(1)=sh_ss-ih(1);
sl(1)=(1-sh_ss)-il(1);
dh(1)=0;
dl(1)=0;
rh(1)=0;
rl(1)=0;

%iterate on SIR model equations
for j=1:1:HH
 
    taul(j,1)=pi1*sl(j)*csl1*( ih(j)*cih1 + il(j)*cil1 ) + pi2*sl(j)*csl2*( ih(j)*cih2 + il(j)*cil2 )+pil3*sl(j)*nsl*il(j)*nil + pi4*sl(j)*( ih(j)+il(j) );
    tauh(j,1)=pi1*sh(j)*csh1*( ih(j)*cih1 + il(j)*cil1 ) + pi2*sh(j)*csh2*( ih(j)*cih2 + il(j)*cil2 )+pih3*sh(j)*nsh*ih(j)*nih + pi4*sh(j)*( ih(j)+il(j) );
         
    sh(j+1,1)=sh(j)-tauh(j);
    sl(j+1,1)=sl(j)-taul(j);
    
    ih(j+1,1)=(1-pirh-pidh)*ih(j)+tauh(j);
    il(j+1,1)=(1-pirl-pidl)*il(j)+taul(j);
    
    rh(j+1,1)=rh(j)+pirh*ih(j);
    rl(j+1,1)=rl(j)+pirl*il(j);
    
    dh(j+1,1)=dh(j)+pidh*ih(j);
    dl(j+1,1)=dl(j)+pidl*il(j);
    
end

err=NaN*zeros(5,1);

%calculate equation residuals for target equations
den=pi1*(sl_ss*csl1+sh_ss*csh1)*(sh_ss*cih1+sl_ss*cil1)+pi2*(sl_ss*csl2+sh_ss*csh2)*(sh_ss*cih2+sl_ss*cil2)+pil3*sl_ss*nsl*nil*sl_ss+pih3*sh_ss*nsh*nih*sh_ss+pi4;
err(1)=pi12_shr_target-(pi1*(sl_ss*csl1+sh_ss*csh1)*(sh_ss*cih1+sl_ss*cil1)+pi2*(sl_ss*csl2+sh_ss*csh2)*(sh_ss*cih2+sl_ss*cil2))/den;
err(2)=pi3_shr_target-(pil3*sl_ss*nsl*nil*sl_ss+pih3*sh_ss*nsh*nih*sh_ss)/den;
err(3)=pil3*sl_ss*nsl*nil*sl_ss-pi3lvshfactor*pih3*sh_ss*nsh*nih*sh_ss;
err(4)=pi1*(sl_ss*csl1+sh_ss*csh1)*(sh_ss*cih1+sl_ss*cil1)-pi1vs2factor*pi2*(sl_ss*csl2+sh_ss*csh2)*(sh_ss*cih2+sl_ss*cil2);
err(5)=RplusD_target-(rh(end)+rl(end)+dh(end)+dl(end));

