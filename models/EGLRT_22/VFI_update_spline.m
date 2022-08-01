%% Update Value Function using Splines for off-grid
function [Vnew,bprime] = VFI_update_spline(V,Y,util,par,mpar,grid,U,flag,solflag,I,mu,Uroprime,Uryprime,Uioprime,Uiyprime,Usoprime,Usyprime,fminbnd_options)
%when solflag==1, solve the model for I=0 for all t. Otherwise, take the
%solution for I=0 and simulate epidemic backwards using perfect foresight.
if solflag==0
    
    bprime = zeros(size(V.ro(mpar.grid_sim_idx)));
    Vnew   = zeros(size(V.ro(mpar.grid_sim_idx)));
    
    %get continuation values -- abusing notation. These are U.io(bprime,iprime)
    U.ro=Uroprime;U.ry=Uryprime;
    U.io=Uioprime;U.iy=Uiyprime;
    U.so=Usoprime;U.sy=Usyprime;
    
    if     flag==1,V=U.ro;
    elseif flag==2,V=U.ry;
    elseif flag==3,V=U.io;
    elseif flag==4,V=U.iy;
    elseif flag==5,V=U.so;
    elseif flag==6,V=U.sy;
    end
    
    
    % Gridded Interpolants for off-grid continuation values
    V_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},V,'spline','pchip'); %continuation value of respective value function (Uro, Ury etc)
    
    if flag>=2 %ry, io, iy, so, sy
        Uro_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},U.ro,'spline','pchip');
    end
    if flag>=4 %iy, so, sy
        Uio_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},U.io,'spline','pchip');
        Ury_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},U.ry,'spline','pchip');
    end
    if flag==6 %sy
        Uiy_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},U.iy,'spline','pchip');
        Uso_int= griddedInterpolant({grid.b(mpar.grid_sim_idx)},U.so,'spline','pchip');
    end
    
    go_through_grid=mpar.grid_sim_idx;
    
else
    
    bprime = zeros(size(V));
    Vnew   = zeros(size(V));
    go_through_grid=1:mpar.nb;
    
    % Gridded Interpolants for off-grid continuation values
    V_int= griddedInterpolant({grid.b},V,'spline','pchip'); %continuation value of respective value function (Uro, Ury etc)
    
    if flag>=2 %ry, io, iy, so, sy
        Uro_int= griddedInterpolant({grid.b},U.ro,'spline','pchip');
    end
    if flag>=4 %iy, so, sy
        Uio_int= griddedInterpolant({grid.b},U.io,'spline','pchip');
        Ury_int= griddedInterpolant({grid.b},U.ry,'spline','pchip');
    end
    if flag==6 %sy
        Uiy_int= griddedInterpolant({grid.b},U.iy,'spline','pchip');
        Uso_int= griddedInterpolant({grid.b},U.so,'spline','pchip');
    end
    
end


%go through grid
for bb=go_through_grid
    if flag==1 %ro
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.deltao)*V_int({b})^(1-par.alf) + par.deltao*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf)   )^par.kap )^par.gam );
    elseif flag==2 %ry
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.deltay-par.nu)*V_int({b})^(1-par.alf) + par.nu*Uro_int({b})^(1-par.alf) + par.deltay*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf) )^par.kap )^par.gam );
    elseif flag==3 %io
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.piro-par.pido)*(1-par.deltao)*V_int({b})^(1-par.alf) + par.piro*(1-par.deltao)*Uro_int({b})^(1-par.alf) + (par.deltao+(1-par.deltao)*par.pido)*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf)    )^par.kap )^par.gam );
    elseif flag==4 %iy
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.piry-par.pidy)*(1-par.deltay-par.nu)*V_int({b})^(1-par.alf) + (1-par.piry-par.pidy)*par.nu*Uio_int({b})^(1-par.alf) + par.piry*(1-par.deltay-par.nu)*Ury_int({b})^(1-par.alf) + par.piry*par.nu*Uro_int({b})^(1-par.alf)  + (par.deltay+(1-par.deltay)*par.pidy)*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf)   )^par.kap )^par.gam );
    elseif flag==5 %so
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.pi1*((Y(bb)-b))*I-par.pi2*I)*(1-par.deltao)*V_int({b})^(1-par.alf) + (par.pi1*((Y(bb)-b))*I+par.pi2*I)*(1-par.deltao)*Uio_int({b})^(1-par.alf)  + par.deltao*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf)   )^par.kap )^par.gam );
    elseif  flag==6 %sy
        f             = @(b)(-par.z  - ( (1-mu)*util((Y(bb)-b)) + (1-(1-mu)*(1-par.betta))*( (1-par.pi1*((Y(bb)-b))*I-par.pi2*I)*(1-par.deltay-par.nu)*V_int({b})^(1-par.alf) + (1-par.pi1*((Y(bb)-b))*I-par.pi2*I)*par.nu*Uso_int({b})^(1-par.alf) + (par.pi1*((Y(bb)-b))*I+par.pi2*I)*(1-par.deltay-par.nu)*Uiy_int({b})^(1-par.alf) + (par.pi1*((Y(bb)-b))*I+par.pi2*I)*par.nu*Uio_int({b})^(1-par.alf) + par.deltay*(par.xi+par.theta*b^(1-par.rho))^(1-par.alf)    )^par.kap )^par.gam );
    end
    %fminbnd_options=optimset('Display','off','MaxFunEvals',100,'TolX',1e-12);   
    %tic
    [bp,v,exitflag]     = fminbnd(f,0,Y(bb),fminbnd_options);
    %toc
        
    if exitflag~=1
        keyboard
    end
    
    Vnew(bb-(go_through_grid(1)-1))   = -v;
    bprime(bb-(go_through_grid(1)-1)) = bp;
    
end