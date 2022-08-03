function [Smat,Imat,Rmat,Dmat,Cmat,Nmat]...
         = SIRmodel_Jgroups_func(betap,nperiods,nperiods_total,irate,Ninit,gammap,thetap,deltap,omegap);
% 
% global nperiods nperiods_total irate Ninit...
%        gammap thetap deltap omegap

   
Jgroups = length(Ninit);   
                
%initialize groups
Smat = zeros(Jgroups,nperiods_total);
Imat = zeros(Jgroups,nperiods_total);
Rmat = zeros(Jgroups,nperiods_total);
Dmat = zeros(Jgroups,nperiods_total);
Cmat = zeros(Jgroups,nperiods_total);
Nmat = zeros(Jgroups,nperiods_total);

Imat(:,1) = irate*Ninit;
Smat(:,1) = Ninit-Imat(:,1);
Nmat(:,1) = Ninit;

% evolution before vaccine
for t = 2:1:nperiods
    for j = 1:Jgroups
        Smat(j,t) =  + Smat(j,t-1) ...
                     - betap(j,:,t-1)*(Imat(:,t-1)./Nmat(:,t-1))*Smat(j,t-1);  

        Imat(j,t) =  + Imat(j,t-1) ...
                     + betap(j,:,t-1)*(Imat(:,t-1)./Nmat(:,t-1))*Smat(j,t-1)...
                     - gammap(j)*Imat(j,t-1);
                 %... 
                 %    - (deltap(j)+omegap(j)*sum(Imat(:,t-1)))*Imat(j,t-1);   

        Rmat(j,t) = Rmat(j,t-1)+gammap(j)*Imat(j,t-1)-thetap(j)*Rmat(j,t-1);

        Dmat(j,t) = Dmat(j,t-1) + deltap(j)*thetap(j)*Rmat(j,t-1);

        Cmat(j,t) =  Cmat(j,t-1) + (1-deltap(j)*thetap(j))*Rmat(j,t-1);

        Nmat(j,t) =  Ninit(j)-Dmat(j,t);   
    end    
end   

% evolution afer vaccine
for t = nperiods+1:1:nperiods_total
    for j = 1:Jgroups
        Smat(j,t) =  + Smat(j,t-1);  
                            
        Imat(j,t) =  + Imat(j,t-1) ...
                     - gammap(j)*Imat(j,t-1);
                 %... 
                 %    - (deltap(j)+omegap(j)*sum(Imat(:,t-1)))*Imat(j,t-1);    

        Rmat(j,t) =  + Rmat(j,t-1)+gammap(j)*Imat(j,t-1)-thetap(j)*Rmat(j,t-1);

        Dmat(j,t) = Dmat(j,t-1) + deltap(j)*thetap(j)*Rmat(j,t-1);

        Cmat(j,t) =  Cmat(j,t-1) + (1-deltap(j)*thetap(j))*Rmat(j,t-1);

        Nmat(j,t) =  Ninit(j)-Dmat(j,t);   
    end
end 