%% mini functions
classdef f

    methods(Static)
        
        function ff = chi(force,pchi)
%            ff = pchi(2)*log(1+pchi(1).*force);
%            ff = (pchi(1).*force).^pchi(2);
            ff = pchi.*force;
        end
       
%        function ff = mit(time,time0,ld0,ld1)
%%            ff = ld0 + ld1*(time-time0);
%%            ff = 1./(1.0+exp(ff));
%
%%			ff = (time-time0)/365/ld1;
%%			ff = ld0*(1.-ff);
%
        function ff = mit(time,time0,time1,ld)
			ff = max(time-time0,0.)./max(time1-time,0.);
			ff = 1. + ff.^ld;
			ff = 1./ff;
%			ff = max(0.,min(ff,1.));
        end
        
        function ff = getpsi(shift,p,wages,gdp0,gdpd)
            p.psi = shift.*p.psi0 + (1.-shift).*ones(p.Ns,p.No);
            
            gdp1 = (1.-p.lockj.*(1-p.psi)).*p.L0;
            gdp1 = [sum(gdp1(:,1:2),2),gdp1(:,3)];
            gdp1 = gdp1.*wages;
            gdp1 = sum(gdp1,2);

            gdp1 = gdp1./gdp0*1d2;
            ff = gdp1 - gdpd;
        end
        
    end
    
end
