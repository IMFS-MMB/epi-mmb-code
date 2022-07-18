function time = solvemodel(time,Ia,p);

% first period is safe
% prices are fixed from last period
    PP = time(1).P;
    price = time(1).price;

%% time iteration from here
    for t = 1:p.Nt

%% 1. compute earnings/utilities/home choice
        [time(t),util] = fwditer0(time(t),PP,price,Ia,p);

% set prices to last period
        PP = time(t).P;
        price = time(t).price;

%% 2. SIR update and testing

        if t<p.ABt
            p.AB = 0;
        else
            p.AB = 1;
        end
        [time(t),Ia] = fwditer1(time(t),p);

%% 3. occupation choice: takes utilites as given
% no infection utilities
        dchoice = time(t).choice;
        home = util.home;
        offc = util.offc;

% aggregate forces of infection
        oforce = p.oq*Ia;
        qforce = p.yq*Ia;
    	if (time(t).Q>0)
            oforce = oforce *(1.0-time(t).q);
            qforce = qforce *(1.0-time(t).q);
        end

% initialize groups
        Nbar = sum(time(t).y2(:,:,:,:),[3,4]);
        Ifrac = sum(time(t).y2(:,2:3,2,:),[3 4])./Nbar(:,2:3);
        Ifrac(isnan(Ifrac)) = 0.0;

% subjective infection rates
%        yforce = p.yb(:,:,1).*Ifrac(:,1);
%        yforce = p.yb(:,:,2).*Ifrac(:,2) + yforce;
%        yforce = (1.0-sum(p.yb,3))*Ia + yforce;
%        yforce = [ones(p.Ns,1)*Ia,yforce];

%        yforce = p.yv.*yforce;        
        yforce = p.yv.*Ia;        
        
% infection disutilities    
        if (p.c<2) 
%            home = home - f.chi(qforce,).*ones(p.Ns,p.No,p.NE,p.Ne);
%            offc = offc - f.chi(yforce,p.chi).*ones(p.Ns,p.No,p.NE,p.Ne);
            home(:,:,:,1:4) = home(:,:,:,1:4) - ...
                f.chi(qforce,p.chi).*ones(p.Ns,p.No,p.NE,4);
            offc(:,:,:,1:4) = offc(:,:,:,1:4) - ...
                f.chi(yforce,p.chi).*ones(p.Ns,p.No,p.NE,4);
        end
            
        val = dchoice.*home + (1.0-dchoice).*offc;
        time(t).util = val;
        
        if (t==p.Nt)
            break
        end
        
% transition matrices
% mass of people who move
        move = zeros(p.Ns,p.No,p.NE,p.Ne);
        move(:,1,:,:) = 1.0/(1.0+exp(val(:,2,:,:)-val(:,1,:,:) + p.loc(:,2)) ...
                                +exp(val(:,3,:,:)-val(:,1,:,:) + p.loc(:,3)));
        move(:,2,:,:) = 1.0/(1.0+exp(val(:,3,:,:)-val(:,2,:,:) + p.loc(:,3)-p.loc(:,2)) ...
                                +exp(val(:,1,:,:)-val(:,2,:,:) - p.loc(:,2)));
        move(:,3,:,:) = 1.0/(1.0+exp(val(:,1,:,:)-val(:,3,:,:) - p.loc(:,3)) ...
                                +exp(val(:,2,:,:)-val(:,3,:,:) + p.loc(:,2)-p.loc(:,3)));

        trans = zeros(p.No,p.No);
        new = zeros(size(move));
        for ie = 1:p.Ne
            for iE = 1:p.NE
                for is = 1:p.Ns
                    trans(:,1) = squeeze(p.trans(is,1,:))'.*squeeze(move(is,:,iE,ie));
                    trans(:,2) = squeeze(p.trans(is,2,:))'.*squeeze(move(is,:,iE,ie));
                    trans(:,3) = squeeze(p.trans(is,3,:))'.*squeeze(move(is,:,iE,ie));
                    
                    trans(1,1) = 1.0-trans(2,1)-trans(3,1);
                    trans(2,2) = 1.0-trans(3,2)-trans(1,2);
                    trans(3,3) = 1.0-trans(1,3)-trans(2,3);
                    new(is,:,iE,ie) = trans*squeeze(time(t).y2(is,:,iE,ie)');
                end
            end
        end
        
        time(t+1).y0 = new;
        time(t+1).o0 = time(t).o2;

    end
        
end
