%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate model thru epidemic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[neglogpost,cmat,comat,cymat,cons_monthly,pidy,pido]=simulate_model(guess,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options)


par.muvec_scale=abs(guess(1));
par.pidy_ini=abs(guess(2));
par.pido_ini=abs(guess(3));
par.gyoung=abs(guess(4));
par.gold=abs(guess(5));

try
    
    
    %scale containment measure
    muvec=muvec*par.muvec_scale;
    
    %constant gain belief evolution
    pidy=Ivec*NaN; %initial cfr
    pido=Ivec*NaN; %initial cfr
    
    pidy(1)=par.pidy_ini;
    pido(1)=par.pido_ini;
    
    for rr=2:1:numel(Ivec)
        pidy(rr)=pidy(rr-1)+par.gyoung*(par.pidy-pidy(rr-1));
        pido(rr)=pido(rr-1)+par.gold*(par.pido-pido(rr-1));
    end
    
    %time varying probabilities of dying and recovering
    par.pidy_vec=pidy;%*0+7*0.001/par.days;
    par.piry_vec=7*1/par.days-par.pidy_vec;
    par.pido_vec=pido;%*0+7*0.035/par.days;
    par.piro_vec=7*1/par.days-par.pido_vec;
    
    if par.piry_vec<0,error('piry <0');end
    if par.piro_vec<0,error('piro <0');end
    
    
    %Store original variables
    Ivec_base=Ivec;
    muvec_base=muvec;
    
    par.pidy_vec_base=par.pidy_vec;
    par.piry_vec_base=par.piry_vec;
    par.pido_vec_base=par.pido_vec;
    par.piro_vec_base=par.piro_vec;
    
    Umat_base=Umat;
    bprimemat_base=bprimemat;
    
    
    if mpar.infect_surprise==0 %perfect foresight for infections and containment
        
        %simulate value functions of agents by age and health types
        sim_val_fun;
        
        %store assets first wave
        bprimemat_foresight=bprimemat;
        
    elseif mpar.infect_surprise==1  %simulate first wave until period 21 (assume steady state in period 22). Then surprise in period 22 that second wave takes place.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%First wave
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %first wave -- take relevant stuff for first wave
        Ivec=Ivec(1:mpar.infect_surprise_week);
        Ivec=[Ivec,0];
        
        muvec=muvec(1:mpar.infect_surprise_week+1);
        
        par.pidy_vec=par.pidy_vec(1:mpar.infect_surprise_week+1);
        par.piry_vec=par.piry_vec(1:mpar.infect_surprise_week+1);
        par.pido_vec=par.pido_vec(1:mpar.infect_surprise_week+1);
        par.piro_vec=par.piro_vec(1:mpar.infect_surprise_week+1);
        
        Umat=Umat_base;
        bprimemat=bprimemat_base;
        
        Umat.ro=Umat.ro(:,end-mpar.infect_surprise_week:end);
        Umat.ry=Umat.ry(:,end-mpar.infect_surprise_week:end);
        Umat.io=Umat.io(:,end-mpar.infect_surprise_week:end);
        Umat.iy=Umat.iy(:,end-mpar.infect_surprise_week:end);
        Umat.so=Umat.so(:,end-mpar.infect_surprise_week:end);
        Umat.sy=Umat.sy(:,end-mpar.infect_surprise_week:end);
        
        bprimemat.ro=bprimemat.ro(:,end-mpar.infect_surprise_week:end);
        bprimemat.ry=bprimemat.ry(:,end-mpar.infect_surprise_week:end);
        bprimemat.io=bprimemat.io(:,end-mpar.infect_surprise_week:end);
        bprimemat.iy=bprimemat.iy(:,end-mpar.infect_surprise_week:end);
        bprimemat.so=bprimemat.so(:,end-mpar.infect_surprise_week:end);
        bprimemat.sy=bprimemat.sy(:,end-mpar.infect_surprise_week:end);
        
        %simulate value functions of agents by age and health types
        sim_val_fun;
        
        %store assets first wave
        bprimemat_first_wave=bprimemat;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Second wave
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %re-set variables to base
        Ivec=Ivec_base;
        muvec=muvec_base;
        
        par.pidy_vec=par.pidy_vec_base;
        par.piry_vec=par.piry_vec_base;
        par.pido_vec=par.pido_vec_base;
        par.piro_vec=par.piro_vec_base;
        
        Umat=Umat_base;
        bprimemat=bprimemat_base;
        
        %second wave -- take relevant stuff for second wave
        Ivec=Ivec(mpar.infect_surprise_week+1:end);
        muvec=muvec(mpar.infect_surprise_week+1:end);
        
        par.pidy_vec=par.pidy_vec(mpar.infect_surprise_week+1:end);
        par.piry_vec=par.piry_vec(mpar.infect_surprise_week+1:end);
        par.pido_vec=par.pido_vec(mpar.infect_surprise_week+1:end);
        par.piro_vec=par.piro_vec(mpar.infect_surprise_week+1:end);
        
        Umat.ro=Umat.ro(:,mpar.infect_surprise_week+1:end);
        Umat.ry=Umat.ry(:,mpar.infect_surprise_week+1:end);
        Umat.io=Umat.io(:,mpar.infect_surprise_week+1:end);
        Umat.iy=Umat.iy(:,mpar.infect_surprise_week+1:end);
        Umat.so=Umat.so(:,mpar.infect_surprise_week+1:end);
        Umat.sy=Umat.sy(:,mpar.infect_surprise_week+1:end);
        
        bprimemat.ro=bprimemat.ro(:,mpar.infect_surprise_week+1:end);
        bprimemat.ry=bprimemat.ry(:,mpar.infect_surprise_week+1:end);
        bprimemat.io=bprimemat.io(:,mpar.infect_surprise_week+1:end);
        bprimemat.iy=bprimemat.iy(:,mpar.infect_surprise_week+1:end);
        bprimemat.so=bprimemat.so(:,mpar.infect_surprise_week+1:end);
        bprimemat.sy=bprimemat.sy(:,mpar.infect_surprise_week+1:end);
        
        %simulate value functions of agents by age and health types
        sim_val_fun;
        
        %re-set variables to base
        Ivec=Ivec_base;
        muvec=muvec_base;
        
        par.pidy_vec=par.pidy_vec_base;
        par.piry_vec=par.piry_vec_base;
        par.pido_vec=par.pido_vec_base;
        par.piro_vec=par.piro_vec_base;
        
        
        %store assets second wave
        bprimemat_second_wave=bprimemat;
        
        %put together assets from both waves
        bprimemat.ro=[bprimemat_first_wave.ro(:,1:mpar.infect_surprise_week) bprimemat_second_wave.ro];
        bprimemat.ry=[bprimemat_first_wave.ry(:,1:mpar.infect_surprise_week) bprimemat_second_wave.ry];
        bprimemat.io=[bprimemat_first_wave.io(:,1:mpar.infect_surprise_week) bprimemat_second_wave.io];
        bprimemat.iy=[bprimemat_first_wave.iy(:,1:mpar.infect_surprise_week) bprimemat_second_wave.iy];
        bprimemat.so=[bprimemat_first_wave.so(:,1:mpar.infect_surprise_week) bprimemat_second_wave.so];
        bprimemat.sy=[bprimemat_first_wave.sy(:,1:mpar.infect_surprise_week) bprimemat_second_wave.sy];
        
    end
    
    
    
    %get consumption
    cmat.ro=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.ro);
    cmat.ry=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.ry);
    cmat.io=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.io);
    cmat.iy=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.iy);
    cmat.so=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.so);
    cmat.sy=(par.w+(1+par.r)*grid.b(mpar.grid_sim_idx)'-bprimemat.sy);
    
    %per capita consumption of young and old
    cymat=(cmat.sy.*Sy'+cmat.iy.*Iy'+cmat.ry.*Ry')./(Sy'+Iy'+Ry');
    comat=(cmat.so.*So'+cmat.io.*Io'+cmat.ro.*Ro')./(So'+Io'+Ro');
    
    %put mpar.timestart periods with pre-epi steady state consumption at
    %the beginning, i.e. before the simulations start
    cymat=[repmat(cymat(:,end),1,mpar.timestart-1) cymat];
    comat=[repmat(comat(:,end),1,mpar.timestart-1) comat];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Epi Simulation results: consumption old and young
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cytmp=(cymat./cmat.sy(:,end)-1)*100;
    cotmp=(comat./cmat.so(:,end)-1)*100;
    
    
    % Put data into timetable
    cons_weekly.Cy=[cytmp(find(par.grid_b_median_assets_idx==mpar.grid_sim_idx),1:numel(Ivec)+mpar.timestart-1)';0];
    cons_weekly.Co=[cotmp(find(par.grid_b_median_assets_idx==mpar.grid_sim_idx),1:numel(Ivec)+mpar.timestart-1)';0];
    
    tty = timetable(newdeaths_weekly.Time,cons_weekly.Cy,'VariableNames',{'Consy'});
    tto = timetable(newdeaths_weekly.Time,cons_weekly.Co,'VariableNames',{'Conso'});
    
    cons_monthly.Cy=retime(tty,'monthly','mean');
    cons_monthly.Co=retime(tto,'monthly','mean');
    
    %March 1 2020 thru April 30 2021
    disty=(dat_young(1:14)-cons_monthly.Cy.Consy(1:14)')/100; %inverse V weighting matrix is in original units (i.e. not in percent)
    disto=(dat_old(1:14)-cons_monthly.Co.Conso(1:14)')/100;
    
    dist=[disty';disto'];
    criterion=dist'*par.invVVmat*dist;
    loglik=numel(dist)/2*log(1/2/pi)-1/2*par.logdetVVmat-1/2*criterion; %in logs
    
    %uniform priors for all parameters except pi
    logprior=0;
    
    for gg=1:1:numel(guess)
        if gg==1 %muscale -- uniform prior
            logprior=logprior+log(unifpdf(abs(guess(gg)),par.LB(gg),par.UB(gg)));
            
        elseif gg==2 %xini young -- uniform prior
            logprior=logprior+log(unifpdf(abs(guess(gg)),par.LB(gg),par.UB(gg)));
            
        elseif gg==3 %xini old -- uniform prior
            logprior=logprior+log(unifpdf(abs(guess(gg)),par.LB(gg),par.UB(gg)));
            
        elseif gg==4 %gyoung -- uniform prior
            logprior=logprior+log(unifpdf(abs(guess(gg)),par.LB(gg),par.UB(gg)));
            
        elseif gg==5 %gold -- uniform prior
            logprior=logprior+log(unifpdf(abs(guess(gg)),par.LB(gg),par.UB(gg)));
        end
    end
    
    logpost=loglik+logprior;
    
    neglogpost=-logpost; %flip sign for maximization using fmincon
    
    %if isinf(logpost)
    %    keyboard
    %end
    
catch
    
    neglogpost=1e5;
    %keyboard
    
end