%this file implements a random walk metropolis mcmc algorithm

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%MCMC%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some mcmc options
n_draws=10000;%number of draws per chain
burn_in_draws=1000;%burn-in draws per chain
n_chains=11;%number of chains

changes=0; %counter for accepted moves in the chain
acceptance_rate=zeros(1,n_chains);
c_val=1.5; %scaling factor

%Scaled inverse hessian
Scaled_invhessian = c_val^2*invhessian;

%initialize matrices to store results in as well as counters
Big_theta=NaN*zeros(numel(para_mode),n_draws,n_chains);
Big_theta_logpost=NaN*zeros(1,n_draws,n_chains);

w = waitbar(0,'MCMC in progress ... please wait ...');

% Create DataQueue and listener
D = parallel.pool.DataQueue;
afterEach(D,@parforWaitbar);
parforWaitbar(w,n_draws*n_chains)

parfor i_chains=1:n_chains
    %Draw initial parameter vector from proposal density
    %note that we assume a multivariate normal density as in Metropolis (1953)
    theta_ini=mvnrnd(para_mode, Scaled_invhessian)';
    
    %get posterior etc at theta_init
    [neglogpost,cmat,comat,cymat,cons_monthly,cfryoung_constgain,cfrold_constgain]=simulate_model(theta_ini,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
    logpost_old = -neglogpost;
    
    theta_ini_mat(:,i_chains)=theta_ini;
    logpost_ini_vec(:,i_chains)=logpost_old;
    
    changes=0; %record changes in log posterior
    
    tmp=NaN*zeros(numel(para_mode),n_draws); %store parameters
    tmp2=NaN*zeros(1,n_draws); %store log post
    
    for i_loop=1:1:n_draws
        send(D,[]); %waitbar
        
        %draw new parameter candidate from proposal density
        cand_theta= mvnrnd(theta_ini,Scaled_invhessian)';
        
        %get new posterior, if possible
        try
            [neglogpost,cmat,comat,cymat,cons_monthly,cfryoung_constgain,cfrold_constgain]=simulate_model(cand_theta,Ivec,muvec,par,U,Y,util,mpar,grid,solflag,Umat,bprimemat,Sy,Ry,Iy,So,Ro,Io,newdeaths_weekly,dat_young,dat_old,dt_cons,fminbnd_options);
        catch
            neglogpost=1e50;
        end
        
        %flip sign
        logpost_new = -neglogpost;
        
        %posterior ratio
        ratio = logpost_new - logpost_old;
        
        %minimum probability of acceptance
        decide_min=min([exp(ratio) 1]);
        
        %draw uniform random variable and decide whether to keep or
        %trash the draw
        ra=rand(1);
        if ra < decide_min
            tmp(1:numel(para_mode),i_loop)=cand_theta;%take new parameter vector
            tmp2(:,i_loop)=logpost_new;
            theta_ini = cand_theta;
            logpost_old=logpost_new;
            changes = changes +1;
        else
            tmp(1:numel(para_mode),i_loop)=theta_ini;%keep old parameter vector
            tmp2(:,i_loop)=logpost_old;
        end
        acceptance_rate(i_chains)=changes/i_loop;
    end
    Big_theta(:,:,i_chains)=tmp;
    Big_theta_logpost(:,:,i_chains)=tmp2;
end

delete(w);

save MCMCdata Big_theta Big_theta_logpost

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Results reporting below%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Big_theta_all=reshape(Big_theta(:,burn_in_draws+1:n_draws,:),numel(para_mode),(n_draws-burn_in_draws)*n_chains);
Big_theta_logpost_all=reshape(Big_theta_logpost(:,burn_in_draws+1:n_draws,:),1,(n_draws-burn_in_draws)*n_chains);

%find posterior_mode
[ww,maxlogpostidx]=max(Big_theta_logpost_all);

%Report Results etc
size_plot_y=ceil(sqrt(size(para_mode,1)));
size_plot_x=round(sqrt(size(para_mode,1)));

%chains
for uu=1:1:n_chains
    figure('Name',['Chains (after burn-in) -- Chain ',num2str(uu)]);
    for mmm=1:1:size(para_mode,1)
        subplot(size_plot_y,size_plot_x,mmm)
        plot(squeeze(Big_theta(mmm,burn_in_draws+1:n_draws,uu)));
        title(par.guess_name(mmm),'Interpreter','none');
        drawnow
        %axis tight
    end
end
orient landscape

%some useful information
disp(' ');
disp(['n_chains: ',num2str(n_chains)]);
disp(['n_draws per chain: ',num2str(n_draws)]);
disp(['burn_in_draws per chain: ',num2str(burn_in_draws)]);
acceptance_rate
disp(' ');

%prior posterior plots
figure('name','Prior-Posterior');
ia=2;ib=3;
trunc = 1e-12; steps = 2^15-1;
for zz=1:1:numel(sol)
    
    if zz==1 %muscale -- uniform prior
        gg_prior=linspace(par.LB(zz),par.UB(zz),steps); 
        prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
        
    elseif zz==2 %pid_ini young -- uniform prior        
        gg_prior=linspace(par.LB(zz),par.UB(zz),steps); 
        prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
        
    elseif zz==3 %pid_xini old -- uniform prior 
        gg_prior=linspace(par.LB(zz),par.UB(zz),steps); 
        prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));
        
    elseif zz==4 %gyoung -- uniform prior 
        gg_prior=linspace(par.LB(zz),par.UB(zz),steps);
        prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));   
                
    elseif zz==5 %gold -- uniform prior 
        gg_prior=linspace(par.LB(zz),par.UB(zz),steps);
        prior_dens=unifpdf(gg_prior,par.LB(zz),par.UB(zz));       
        
    end
    
    [F,XI]=ksdensity(Big_theta_all(zz,:));
    subplot(ia,ib,zz)
    plot(gg_prior,prior_dens,'k--','linewidth',1.5); hold on
    plot(XI,F,'b-','linewidth',2.5);
    vline(Big_theta_all(zz,maxlogpostidx));
    axis tight;
    title(char(par.guess_name(zz)),'Interpreter','none');
end
suptitle('Priors vs. Posteriors');
legend1=legend('Prior','Posterior (MCMC)');
set(legend1,...
    'Position',[0.418594533212181 0.906613971512973 0.174041297935103 0.0217391304347826],...
    'Orientation','horizontal');
legend boxoff;orient landscape;
print -dpdf -fillpage prior_posterior_mcmc

%print para estimates
for qq=1:1:numel(para_mode)
    par.guess_name(qq)
    [post_mean, post_median, post_var, hpd_interval] = posterior_moments(Big_theta_all(qq,:),1,0.95)
end

%print out mode
disp(' ');
disp(['         Max. Log Posterior (Optimization)     ','Max. Log Posterior (MCMC)']);
fprintf('                      %6.6f                     %6.6f \n', logpost_mode,Big_theta_logpost_all(maxlogpostidx))
disp(' ');

for ii=1:1:size(para_mode)
    fprintf('%s              %6.4f                       %6.4f \n',char(par.guess_name(ii,:)), para_mode(ii),Big_theta_all(ii,maxlogpostidx))
end

%print posterior mode parameters
post_mode_mcmc=Big_theta_all(:,maxlogpostidx)

%time for MCMC
time=toc;
disp(' ');
disp(['Total time for MCMC: ',num2str(time/3600),' hours']);

save mcmc_results
