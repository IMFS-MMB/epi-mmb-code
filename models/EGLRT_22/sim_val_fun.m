%recovered agents
for flag = 1:2
    
    for TT=numel(Ivec)-1:-1:1
        
        I=Ivec(TT);       %infections at inner time loop
        mu=muvec(TT);     %containment at inner time loop
        
        par.pidy=7*0.001/par.days;        %Weekly probability of dying, young
        par.piry=7*1/par.days-par.pidy;   %Weekly probability of recovering, young
        par.pido=7*0.035/par.days;        %Weekly probability of dying, old
        par.piro=7*1/par.days-par.pido;   %Weekly probability of recovering, old
        
        [Unewt,bprimet] = VFI_update_spline(U,Y,util,par,mpar,grid,[],flag,solflag,I,mu,Umat.ro(:,TT+1),Umat.ry(:,TT+1),Umat.io(:,TT+1),Umat.iy(:,TT+1),Umat.so(:,TT+1),Umat.sy(:,TT+1),fminbnd_options);
        
        if     flag==1, Umat.ro(:,TT)=Unewt;bprimemat.ro(:,TT)=bprimet;
        elseif flag==2, Umat.ry(:,TT)=Unewt;bprimemat.ry(:,TT)=bprimet;
        end
        
    end
end

Umat_tmp.io=Umat.io;
Umat_tmp.iy=Umat.iy;

%infected agents
for flag = 3:4
    
    for TT=numel(Ivec)-1:-1:1%, TT   %outer time loop. go backwards in time.
        
        for TTrec=numel(Ivec)-1:-1:TT %inner time loop. Always start at the end and go back to time in the outer time loop
            
            I=Ivec(TTrec);       %infections at inner time loop
            mu=muvec(TTrec);     %containment at inner time loop
            
            par.pidy=par.pidy_vec(TT);  %mortality and recovery rates at outer time loop, i.e. constant rates until end of time; implies random walk forecast of beliefs about mortality rates..
            par.piry=par.piry_vec(TT);
            par.pido=par.pido_vec(TT);
            par.piro=par.piro_vec(TT);
            
            if flag==3
                [Unewt,bprimet] = VFI_update_spline(U,Y,util,par,mpar,grid,[],flag,solflag,I,mu,Umat.ro(:,TTrec+1),Umat.ry(:,TTrec+1),Umat_tmp.io(:,TTrec+1),Umat_tmp.iy(:,TTrec+1),Umat.so(:,TTrec+1),Umat.sy(:,TTrec+1),fminbnd_options);
            elseif flag==4
                [Unewt,bprimet] = VFI_update_spline(U,Y,util,par,mpar,grid,[],flag,solflag,I,mu,Umat.ro(:,TTrec+1),Umat.ry(:,TTrec+1),Umat.io(:,TTrec+1),Umat_tmp.iy(:,TTrec+1),Umat.so(:,TTrec+1),Umat.sy(:,TTrec+1),fminbnd_options);
            end
            
            if     flag==3, Umat_tmp.io(:,TTrec)=Unewt;
            elseif flag==4, Umat_tmp.iy(:,TTrec)=Unewt;
            end
            
        end
        
        if     flag==3,  bprimemat.io(:,TT)=bprimet;Umat.io(:,TT)=Unewt;
        elseif flag==4,  bprimemat.iy(:,TT)=bprimet;Umat.iy(:,TT)=Unewt;
        end
        
    end
end


%susceptible agents
for flag = 5:6
    
    for TT=numel(Ivec)-1:-1:1
        
        I=Ivec(TT);       %infections at inner time loop
        mu=muvec(TT);     %containment at inner time loop
        
        par.pidy=7*0.001/par.days;        %Weekly probability of dying, young
        par.piry=7*1/par.days-par.pidy;   %Weekly probability of recovering, young
        par.pido=7*0.035/par.days;        %Weekly probability of dying, old
        par.piro=7*1/par.days-par.pido;   %Weekly probability of recovering, old
        
        [Unewt,bprimet] = VFI_update_spline(U,Y,util,par,mpar,grid,[],flag,solflag,I,mu,Umat.ro(:,TT+1),Umat.ry(:,TT+1),Umat.io(:,TT+1),Umat.iy(:,TT+1),Umat.so(:,TT+1),Umat.sy(:,TT+1),fminbnd_options);
        
        if     flag==5, Umat.so(:,TT)=Unewt;bprimemat.so(:,TT)=bprimet;
        elseif flag==6, Umat.sy(:,TT)=Unewt;bprimemat.sy(:,TT)=bprimet;
        end
        
    end
end