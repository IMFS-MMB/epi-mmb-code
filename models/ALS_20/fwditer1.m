function [time,Ia] = fwditer1(time,p)

%% initialize and save for testing
% guys who are still confirmed or recovered
    ysir = time.y0;
    yc = sum(ysir(:,:,:,3:4),[3,4]);
    yc1 = (1.0-p.ydelt).*(1.0-p.ygamm).*(1.0-p.ym).*yc;
    
    yr = sum(ysir(:,:,:,5:6),[3,4]);
    yr1 = p.ygamm.*(1.0-p.ym).*yc + yr;
    yr1 = (1.0-p.ydelt).*yr1;
    
    osir = time.o0;
    oc = sum(osir(:,3:4),'all');
    oc1 = (1.0-p.odelt) *(1.0-p.ogamm) *(1.0-p.om) *oc;
    
    or = sum(osir(:,5:6),'all');
    or1 = p.ogamm *(1.0-p.om) *oc + or;
    or1 = (1.0-p.odelt) *or1;
    
% guys who stay home or not
    qsir = ysir .*time.choice;
    ysir = ysir - qsir;

%% true law of motion
%% no need to keep track of confirmed status
    ysir = sum(ysir,4);
    qsir = sum(qsir,4);

%% control for quarantines
    Ia = sum(ysir(:,:,2),'all') + (1.0-time.q)*(sum(qsir(:,:,2),'all') + osir(2));
    Na = sum(ysir,'all') + sum(qsir,'all') + sum(osir,'all');
    Ia = Ia/Na;
    Ia(isnan(Na)) = 0.0;

% aggregate forces of infection
    oforce = p.oq*Ia;
    qforce = p.yq*Ia;
	if (time.Q>0)
%        oforce = oforce *(1.0-time.q);
%    	qforce = qforce *(1.0-time.q);
	end

% the groups: initialize first    
%    Nbar = sum(time.y0(:,:,:,:),[3,4]);
%    Ifrac = sum(time.y0(:,2:3,2,:),4)./Nbar(:,2:3);
%    Ifrac(isnan(Ifrac)) = 0.0;

    Ifrac = ysir(:,2:3,2,:) + (1.0-time.q)*qsir(:,2:3,2,:);
    Nbar = sum(ysir(:,2:3,:)+qsir(:,2:3,:),3);
    Ifrac = Ifrac./Nbar;
    Ifrac(isnan(Ifrac)) = 0.0;

% actual infection rates
% control for quaratines
%    yforce = p.yb(:,:,1).*Ifrac(:,1);
%    yforce = p.yb(:,:,2).*Ifrac(:,2) + yforce;
%    yforce = (1.0-sum(p.yb,3))*Ia + yforce;
%    yforce = [ones(p.Ns,1)*Ia,yforce];

%    yforce = p.yv.*yforce;
    yforce = p.yv.*Ia;
    
% cap infectionable to susceptible
    yforce = min(yforce,1.0);
    qforce = min(qforce,1.0);
    oforce = min(oforce,1.0);
    
%% transition matrices
% young
    trans = zeros(p.Ns,p.No,p.NE,p.NE);

    trans(:,:,2,2) = (1.0-p.ygamm).*(1.0-p.ym).*ones(p.Ns,p.No);
    trans(:,:,3,2) = p.ygamm.*(1.0-p.ym).*ones(p.Ns,p.No);
    trans(:,:,3,3) = 1.0;

    for io = 1:p.No
        for is = 1:p.Ns
            trans(is,io,1,1) = 1.0-yforce(is,io);
            trans(is,io,2,1) = yforce(is,io);

            ysir(is,io,:) = squeeze(trans(is,io,:,:)) ...
                            *squeeze(ysir(is,io,:));

            trans(is,io,1,1) = 1.0-qforce;
            trans(is,io,2,1) = qforce;

            qsir(is,io,:) = squeeze(trans(is,io,:,:)) ...
                            *squeeze(qsir(is,io,:));
        end
    end
    ysir = ysir.*(1.0-p.ydelt);
    qsir = qsir.*(1.0-p.ydelt);

% old 
    osir = sum(osir,2);
    trans = zeros(p.NE,p.NE);

    trans(2,2) = (1.0-p.ogamm).*(1.0-p.om);
    trans(3,2) = p.ogamm.*(1.0-p.om);
    trans(3,3) = 1.0;

    trans(1,1) = 1.0-oforce;
    trans(2,1) = oforce;

    osir = trans * osir;
    osir = osir * (1.0-p.odelt);

%% save and return Ia for next eriod
    time.y1 = [ysir;qsir];
    time.o1 = osir;

    ysir = ysir+qsir;

%    Ia = sum(ysir(:,:,2),'all') + osir(2);
%    Na = sum(ysir,'all') + sum(osir);
%    Ia = Ia/Na;
%    Ia(isnan(Na)) = 0.0;
    

%% testing
    ysir1 = zeros(p.Ns,p.No,p.NE,p.Ne);
    osir1 = zeros(p.NE,p.Ne);
    taua = time.tau(1);
    taus = time.tau(2);
    
    atest = 0;
    stest = 0;

% S: can only be 0
    ysir1(:,:,1,1) = ysir(:,:,1).*(1.0-p.yf);
    ysir1(:,:,1,2) = ysir(:,:,1).*p.yf;

    osir1(1,1) = osir(1).*(1.0-p.of);
    osir1(1,2) = osir(1).*p.of;

    atest = sum(ysir1(:,:,1,1),'all') + osir1(1,1);
    stest = sum(ysir1(:,:,1,2),'all') + osir1(1,2);
    
% I: can be 0/c
% subtract already confirmed and not recovered (Ihat)
    ysir(:,:,2) = ysir(:,:,2) - yc1;
    osir(2) = osir(2) - oc1;

    
% infected and not confirmed / confirmed
    ysir1(:,:,2,1) = ysir(:,:,2).*(1.0-p.yeta);
    ysir1(:,:,2,2) = ysir(:,:,2).*p.yeta;
    ysir1(:,:,2,3) = ysir1(:,:,2,1)*taua + yc1.*(1.0-p.yf);
    ysir1(:,:,2,4) = ysir1(:,:,2,2)*taus + yc1.*p.yf;

    atest = sum(ysir1(:,:,2,1),'all') + atest;
    stest = sum(ysir1(:,:,2,2),'all') + stest;

    ysir1(:,:,2,1) = ysir1(:,:,2,1).*(1.0-taua);
    ysir1(:,:,2,2) = ysir1(:,:,2,2).*(1.0-taus);
    
    osir1(2,1) = osir(2).*(1.0-p.oeta);
    osir1(2,2) = osir(2).*p.oeta;
    osir1(2,3) = osir1(2,1)*taua + oc1 *(1.0-p.of);
    osir1(2,4) = osir1(2,2)*taus + oc1 *p.of;
    
    atest = osir1(2,1) + atest;
    stest = osir1(2,2) + stest;
    
    osir1(2,1) = osir1(2,1).*(1.0-taua);
    osir1(2,2) = osir1(2,2).*(1.0-taus);
    
% R: can be 0/r
% subtract already confirmed and recovered (Rhat)
    ysir(:,:,3) = ysir(:,:,3) - yr1;
    osir(3) = osir(3) - or1;

% recovered and not confirmed / confirmed
    ysir1(:,:,3,1) = ysir(:,:,3).*(1.0-p.yf);
    ysir1(:,:,3,2) = ysir(:,:,3).*p.yf;
    ysir1(:,:,3,5) = yr1.*(1.0-p.yf);
    ysir1(:,:,3,6) = yr1.*p.yf;

	atest = sum(ysir1(:,:,3,1),'all') + atest;
    stest = sum(ysir1(:,:,3,2),'all') + stest;
    if (p.AB==1)
        ysir1(:,:,3,5) = ysir1(:,:,3,1).*taua + ysir1(:,:,3,5);
        ysir1(:,:,3,6) = ysir1(:,:,3,2).*taus + ysir1(:,:,3,6);

        ysir1(:,:,3,1) = ysir1(:,:,3,1).*(1.0-taua);
        ysir1(:,:,3,2) = ysir1(:,:,3,2).*(1.0-taus);
    end

	osir1(3,1) = osir(3).*(1.0-p.of);
    osir1(3,2) = osir(3).*p.of;
    osir1(3,5) = or1.*(1.0-p.of);
    osir1(3,6) = or1.*p.of;

	atest = osir1(3,1) + atest;
    stest = osir1(3,2) + stest;
    if (p.AB==1)
        osir1(3,5) = osir1(3,5) + osir1(3,1).*taua;
        osir1(3,6) = osir1(3,6) + osir1(3,2).*taus;

        osir1(3,1) = osir1(3,1).*(1.0-taua);
        osir1(3,2) = osir1(3,2).*(1.0-taus);
    end
 
%% save post test distributions
    time.y2 = ysir1;
    time.o2 = osir1;

% save to match total tested

    Ia = sum(ysir1(:,:,2,:),'all') + sum(osir1(2,:),'all');
    Na = sum(ysir1,'all') + sum(osir1,'all');
    Ia = Ia/Na;
    Ia(isnan(Na)) = 0.0;

    time.tested = taua*atest + taus*stest;

    
end
