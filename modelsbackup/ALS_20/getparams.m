function [e,p] = pmratwkappa(e,p,maxit,eps);

    for kiter = 1:maxit
%        disp([p.kappa,p.phi(2,2)]);
        
% set parameters
        fs = p.yf;%(:,1);
        fm = p.yf;%(:,2);
        fw = p.yf;%(:,3);

% manager productivity
        meff = ones(p.Ns,2);
        meff(1,2) = p.phi(1,2);
        meff(2,2) = p.phi(2,2).*p.psi(2,2);
        meff = meff.*[p.z(:,2),p.z(:,2)];
    
% workers: always forced to work
        weff = ones(p.Ns,2);
        weff(:,2) = p.phi(:,3);
        weff = weff.*[p.z(:,3),p.z(:,3)];

% get Zmwrat
        tmeff = meff(:,1).*(1.0-fm) + meff(:,2).*fm;
        tweff = weff(:,1).*(1.0-fw) + weff(:,2).*fw;
    
        Zmwrat = tmeff./tweff .* p.L0(:,2)./p.L0(:,3);

% set SE productivity using S-M indifference
        p.z(:,1) = p.alph.*p.z(:,2).*Zmwrat.^(p.alph-1.0);

% SE: don't include productivity
        seff = ones(p.Ns,2);
        seff(1,2) = p.phi(1,1);
        seff(2,2) = p.phi(2,1)*p.psi(2,1);
        seff = seff.*[p.z(:,1),p.z(:,1)];

        tseff = seff(:,1).*(1.0-fs) + seff(:,2).*fs;
    
% since we have all effective units, set eqm rel price
% effective units
        teff(:,1) = tseff.*p.L0(:,1);
        teff(:,2) = tmeff.*p.L0(:,2);
        teff(:,3) = tweff.*p.L0(:,3);

% total production
        prod(:,1) = teff(:,1);
        prod(:,2) = teff(:,2).^p.alph.*teff(:,3).^(1.0-p.alph);

% price
        p21 = (1.0-p.thet)/p.thet * sum(prod(1,:))/sum(prod(2,:));
        px(1) = e.P0*p.thet^p.thet;
        px(1) = px(1)/(p21/(1.0-p.thet))^(1.0-p.thet);
        px(2,1) = px(1)*p21;

% set kappa using max psi
		phi = log(1.+px.*p.z.*p.phi) - log(1.+px.*p.z.*p.phi.*p.psi);
		p.kappa = min(phi,[],'all');

% set phi's so everyone is indifferent
		phi = (1.-exp(p.kappa))./(exp(p.kappa)*p.psi-1.);
		phi = phi./(px.*p.z);

		diff = sum((phi-p.phi).^2,'all');
		
		if (diff<=eps)
			break;
		end
		p.phi = phi;
    end
return
% use gdp to normalize infection probability
%    gdp = sum(prod(1,:)).^p.thet.*sum(prod(2,:)).^(1.0-p.thet);
%    p.chi(1) = exp(log(1+gdp)/p.chi(2))-1;
%    p.chi(1) = log(1+gdp);
%    p.chi(2) = log(1+gdp)/log(1+p.chi(1));

% manager and worker wages
    mx = p.alph.*px.*Zmwrat.^(p.alph-1.0);
    wx = (1.0-p.alph).*px.*Zmwrat.^p.alph;

    e.prod0 = prod;
    e.price0 = [px,mx,wx];
    e.earn0 = zeros(p.Ns,p.No,2);
    e.earn0(:,:,1) = [px.*seff(:,1),mx.*meff(:,1),wx.*weff(:,1)];
    e.earn0(:,:,2) = [px.*seff(:,2),mx.*meff(:,2),wx.*weff(:,2)];
    
end
