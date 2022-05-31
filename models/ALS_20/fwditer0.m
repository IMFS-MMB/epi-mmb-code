function [time,util] = fwditer0(time,PP,price,Ia,p)

% perceived aggregate forces of infection
    qforce = p.yq*Ia;

% the groups: initialize first    
    Nbar = sum(time.y0(:,:,:,:),[3,4]);
    Ifrac = sum(time.y0(:,2:3,2,:),[3 4])./Nbar(:,2:3);
    Ifrac(isnan(Ifrac)) = 0.0;

% subjective infection rates
%    yforce = p.yb(:,:,1).*Ifrac(:,1);
%    yforce = p.yb(:,:,2).*Ifrac(:,2) + yforce;
%    yforce = (1.0-sum(p.yb,3))*Ia + yforce;
%    yforce = [ones(p.Ns,1)*Ia,yforce];

%    yforce = p.yv.*yforce;
    yforce = p.yv.*Ia;
    
% cap infectionable to susceptible
    qforce = min(qforce,1.0);
    yforce = min(yforce,1.0);
    
%% work from home choices
% effective productivities, including sickness
    eff = p.z.*ones(p.Ns,p.No,p.NE,p.Ne);
    eff(:,:,:,p.sick) = p.phi.*eff(:,:,:,p.sick);
    
% earnings
    offc(:,1,:,:) = price(:,1).*eff(:,1,:,:);
    offc(:,2,:,:) = price(:,2).*eff(:,2,:,:);
    offc(:,3,:,:) = price(:,3).*eff(:,3,:,:);

% home discount and sick utilities
    home = p.psi.*offc;
    home = log(1.0+home);	%
    offc = log(1.0+offc);	%1.0+
    offc(:,1:2,:,p.sick) = offc(:,1:2,:,p.sick) - p.kappa;
    
% infection disutilities    
    if (p.c<2)
%        home = home - p.wchi.*f.chi(qforce,p.chi).*ones(p.Ns,p.No,p.NE,p.Ne);
%        offc = offc - p.wchi.*f.chi(yforce,p.chi).*ones(p.Ns,p.No,p.NE,p.Ne);
        home(:,:,:,1:4) = home(:,:,:,1:4) - p.wchi.* ...
            f.chi(qforce,p.chi).*ones(p.Ns,p.No,p.NE,4);
        offc(:,:,:,1:4) = offc(:,:,:,1:4) - p.wchi.* ...
            f.chi(yforce,p.chi).*ones(p.Ns,p.No,p.NE,4);
    end
        
% SE and manager choices    
    ichoice = (home>=offc);

% adjust for continuity
    dchoice = 1.0/(1.0+exp((offc-home)/p.qloc));

    if time.Q==1
        dchoice(:,:,:,2:4) = max(time.q,dchoice(:,:,:,2:4));
    end
% lock-down
    if time.Q==2
		qj = time.j*p.lockj.*ones(p.Ns,p.No,p.NE,p.Ne);
        dchoice = max(qj,dchoice);
        dchoice(:,:,:,2:4) = max(time.q,dchoice(:,:,:,2:4));
    end
% virus visa
    if time.Q==3
		qj = time.j*p.lockj.*ones(p.Ns,p.No,p.NE,p.Ne);
        dchoice(:,:,:,3:4) = max(qj(:,:,:,3:4),dchoice(:,:,:,3:4));
        dchoice(:,:,:,2:4) = max(time.q,dchoice(:,:,:,2:4));
    end
    
%% equilibrium price and wage
% effective units
    temp = p.psi .*dchoice;
    temp = temp + (1.0-dchoice);
    eff = eff.*temp;

    teff = eff.*time.y0;
    teff = sum(teff,[3,4]);

% total production    
    prod(:,1) = teff(:,1);
    prod(:,2) = teff(:,2).^p.alph.*teff(:,3).^(1.0-p.alph);

% relative price update
    p21 = (1.0-p.thet)/p.thet * sum(prod(1,:))/sum(prod(2,:));

    px(1) = PP * p.thet^p.thet;
    px(1) = px(1)/(p21/(1.0-p.thet))^(1.0-p.thet);
    px(2,1) = px(1)*p21;

% manager and worker wage update
    mx = p.alph.*px.*(teff(:,3)./teff(:,2)).^(1.0-p.alph);
    wx = (1.0-p.alph).*px.*(teff(:,2)./teff(:,3)).^p.alph;

    price = [px,mx,wx];
    
%% new earnings for occupation choice
% earnings
    offc(:,1,:,:) = price(:,1).*eff(:,1,:,:);
    offc(:,2,:,:) = price(:,2).*eff(:,2,:,:);
    offc(:,3,:,:) = price(:,3).*eff(:,3,:,:);

% home discount and sick utilities
    home = p.psi.*offc;
    home = log(1.0+home);	%
    offc = log(1.0+offc);	%
    offc(:,:,:,p.sick) = offc(:,:,:,p.sick) - p.kappa;
    
%% save results
    time.choice = dchoice;
    time.earn = eff.*price;
    time.prod = prod;

    time.P = PP;
    time.price = price;

    util.home = home;
    util.offc = offc;
    
end
