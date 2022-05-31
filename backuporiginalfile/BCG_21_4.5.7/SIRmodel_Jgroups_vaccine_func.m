function [Smat,Vmat,Imat,Rmat,Dmat,Cmat,Nmat]...
         = SIRmodel_Jgroups_vaccine_func(betap,rhop,nperiods,irate, Ninit, gammap, thetap, deltap)


   
Jgroups = length(Ninit);   
                
%initialize groups
Smat = zeros(Jgroups,nperiods);
Vmat = zeros(Jgroups,nperiods);
Imat = zeros(Jgroups,nperiods);
Rmat = zeros(Jgroups,nperiods);
Dmat = zeros(Jgroups,nperiods);
Cmat = zeros(Jgroups,nperiods);
Nmat = zeros(Jgroups,nperiods);

Imat(:,1) = irate*Ninit;
Smat(:,1) = Ninit-Imat(:,1);
Vmat(:,1) = 0*Ninit;
Nmat(:,1) = Ninit;

% evolution before vaccine
for t = 2:1:nperiods
    for j = 1:Jgroups
        Smat(j,t) =  + Smat(j,t-1) ...
                     - betap(j,:,t-1)*(Imat(:,t-1)./Nmat(:,t-1))*Smat(j,t-1) ...
                     - rhop(j,t-1)*Smat(j,t-1);  

        Vmat(j,t) =  + Vmat(j,t-1)...
                     + rhop(j,t-1)*Smat(j,t-1);        
                 
        Imat(j,t) =  + Imat(j,t-1) ...
                     + betap(j,:,t-1)*(Imat(:,t-1)./Nmat(:,t-1))*Smat(j,t-1)...
                     - gammap(j)*Imat(j,t-1);

        Rmat(j,t) =  + Rmat(j,t-1)+gammap(j)*Imat(j,t-1)-thetap(j)*Rmat(j,t-1);

        Dmat(j,t) =  + Dmat(j,t-1) + deltap(j)*thetap(j)*Rmat(j,t-1);

        Cmat(j,t) =  + Cmat(j,t-1) + (1-deltap(j))*thetap(j)*Rmat(j,t-1);

        Nmat(j,t) =  + Nmat(j,1)-Dmat(j,t);   
    end    
end   
