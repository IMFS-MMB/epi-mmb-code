function g = plots(data,p,r,name)

legnom3 = {'Model True','Model Observed','Data'};
legnom6 = {'Low-Skill SE','High-Skill SE','Mgr','Mgr','Wkr','Wkr'};

%% SIR PLOTS
g.YSO = zeros(p.Nt,3);
g.SIR = zeros(p.Nt,3);
g.D01 = zeros(p.Nt,2);

g.ucr = g.SIR;
g.d01 = g.SIR;

% true and observed infection rates, total
for t = 1:p.Nt
    g.YSO(t,:) = [sum(data(t).y0(1,:,:,:),'all'), ...
                    sum(data(t).y0(2,:,:,:),'all'), ...
                    sum(data(t).o0,'all')];
    
    g.SIR(t,1) = sum(data(t).y0(:,:,1,:),'all')+sum(data(t).o0(1,:));
    g.SIR(t,2) = sum(data(t).y0(:,:,2,:),'all')+sum(data(t).o0(2,:));
    g.SIR(t,3) = sum(data(t).y0(:,:,3,:),'all')+sum(data(t).o0(3,:));

    g.ucr(t,1) = sum(data(t).y0(:,:,:,1:2),'all')+sum(data(t).o0(:,1:2),'all');
    g.ucr(t,2) = sum(data(t).y0(:,:,:,3:4),'all')+sum(data(t).o0(:,3:4),'all');
    g.ucr(t,3) = sum(data(t).y0(:,:,:,5:6),'all')+sum(data(t).o0(:,5:6),'all');

    if (t==p.Nt)
        break
    end
    
% die with virus
    g.D01(t+1,1) = sum(p.ym.*data(t).y0(:,:,2,:),'all');
    g.d01(t+1,1) = sum(p.ym.*data(t).y0(:,:,:,3:4),'all');
    g.D01(t+1,1) = g.D01(t+1,1) + p.om*sum(data(t).o0(2,:),'all');
    g.d01(t+1,1) = g.d01(t+1,1) + p.om*sum(data(t).o0(:,3:4),'all');

% die because of virus
    g.D01(t+1,2) = sum((1.0-p.ydelt).*p.ym.*data(t).y0(:,:,2,:),'all');
    g.d01(t+1,2) = sum((1.0-p.ydelt).*p.ym.*data(t).y0(:,:,:,3:4),'all');
    g.D01(t+1,2) = g.D01(t+1,2) + (1.0-p.odelt)*p.om*sum(data(t).o0(2,:),'all');
    g.d01(t+1,2) = g.d01(t+1,2) + (1.0-p.odelt)*p.om*sum(data(t).o0(:,3:4),'all');
end
g.N = sum(g.YSO,2);
g.SIR = g.SIR./g.N(1);
g.ucr = g.ucr./g.N(1);

g.D01 = cumsum(g.D01)./g.N(1);
g.d01 = cumsum(g.d01)./g.N(1);
g.d01(:,3) = 1.0 - sum(g.SIR,2);

%% SIR
%h1 = figure('color','w');

g.ucr(:,2) = sum(g.ucr(:,2:3),2);
temp = g.ucr(2:end,2)-g.ucr(1:end-1,2);
g.ucr(:,2) = [0;temp];
%p1 = semilogy(r.stime,g.ucr(:,2)*p.pop); hold on;

g.SIR(:,2) = sum(g.SIR(:,2:3),2);
temp = g.SIR(2:end,2)-g.SIR(1:end-1,2);
g.SIR(:,2) = [0;temp];
%p2 = semilogy(r.stime,g.SIR(:,2)*p.pop); hold on;

temp = r.crd(2:end,1)-r.crd(1:end-1,1);
temp = [0;temp];
for i = 1:size(temp)
	if (temp(i)>0)
		break
	end
end
temp(i:end) = movmean(temp(i:end),7);
% p3 = plot(r.xaxis,temp); hold on;
% 
% ylim([1 Inf])
% 
% vals = {3,3,3;'-','--',':';[.7 .7 .7],[.4 .4 .4],'k'}';
% h = [p2;p1;p3];
% set(h,r.style,vals);
% legend(h,legnom3,'NumColumns',3,'location','southoutside');
% 
% grid on;
% xlim([1 p.Nt])
% xticks(r.xtnum);xticklabels(r.xtlab); 
% set(gca,'fontsize',12'); 


%% deaths
% h2 = figure('color','w');
% 
% %p1 = plot(r.stime,g.d01(r.stime,1:2)*p.pop); hold on;
% %vals = {3,3;'-','-';'b','r'}';
% %set(p1,r.style,vals)
% 
% %p2 = plot(r.stime,g.D01(r.stime,:)*p.pop); hold on;
% %vals = {2,2;'--','--';'b','r'}';
% %set(p2,r.style,vals)
% 
% p1 = semilogy(r.stime,g.d01(r.stime,1)*p.pop); hold on;
% p2 = semilogy(r.stime,g.D01(r.stime,1)*p.pop); hold on;
% 
% %title('- Observed  -- True','fontsize',16);	%...
% %   ,'units','normalized','position',[0.3,0.9,0]);
% p3 = semilogy(r.xaxis,r.crd(:,3),':'); hold on;
% 
% h = [p2;p1;p3];
% set(h,r.style,vals);
% %if (p.uk==0)
% %    title('- Observed  -- True  x Korea','fontsize',16);
% 	legend(h,legnom3,'NumColumns',3,'location','southoutside');
% %    ylim([1 10^6])
% %else
% %    title('- Observed  -- True  x UK','fontsize',16);
% %    legend({'Model Observed','Model True','Data'},'location','southeast','fontsize',16);
% %    ylim([1 10^7])
% %end
% ylim([1 Inf])
% 
% grid on
% xlim([1 p.Nt])
% xticks(r.xtnum);xticklabels(r.xtlab);
% set(gca,'fontsize',12'); 

%% ECONOMIC PLOTS
% employmep.Nt and incomes
g.prod = zeros(p.Nt,p.Ns);
g.emps = zeros(p.Nt,p.Ns,p.No);
g.earn = g.emps;
g.util = g.emps;
g.home = g.emps;
for io = 1:p.No
    for is = 1:p.Ns
        for t = 1:p.Nt
            g.prod(t,:) = sum(data(t).prod,2);

            temp = squeeze(data(t).y0(is,io,:,:));
        	g.emps(t,is,io) = sum(temp,'all');
        	g.earn(t,is,io) = sum(squeeze(data(t).earn(is,io,:,:)).*temp,'all');
        	g.earn(t,is,io) = g.earn(t,is,io)/g.emps(t,is,io);
        	g.home(t,is,io) = sum(squeeze(data(t).choice(is,io,:,:)).*temp,'all');
        	g.home(t,is,io) = g.home(t,is,io)/g.emps(t,is,io);
        	g.util(t,is,io) = sum(squeeze(data(t).util(is,io,:,:)).*temp,'all');
        	g.util(t,is,io) = g.util(t,is,io)/g.emps(t,is,io);
        end
    end
end
g.earn(isnan(g.earn))=0.0;
g.emps = g.emps./sum(g.YSO(:,1:2),2);

%% GDP
%h3 = figure('color','w');

g.prod(:,3) = g.prod(:,1).^p.thet.*g.prod(:,2).^(1.0-p.thet);
g.prod(:,1:2) = g.prod(:,1:2);  %./g.YSO(:,1:2);
g.prod(:,3) = g.prod(:,3);  %./sum(g.YSO(:,1:2),2);
g.prod = log(g.prod)-log(g.prod(1,:));

% p1 = plot(r.stime,g.prod(r.stime,:)); hold on;
% 
% % vals = {2,2,2;'--',':','-';'b','r','k'}';
% vals = {3,3,3;'--',':','-';[.4 .4 .4],[.7 .7 .7],'k'}';
% set(p1,r.style,vals)
% 
% xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
% ylim([-0.30 0.00]);
% 
% legend(p1,{'Low-Skill','High-Skill','Total'},'location','southoutside','NumColumns',3);
% set(gca,'fontsize',12); 
% grid on

%% Earnings
%h4 = figure('color','w');
g.earn = log(g.earn) - log(g.earn(1,:,:));

% p1 = plot(r.stime,g.earn(r.stime,:,1)); hold on;
% p2 = plot(r.stime,g.earn(r.stime,:,2)); hold on;
% p3 = plot(r.stime,g.earn(r.stime,:,3)); hold on;
% xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
% 
% vals = {2,2,2,2,2,2;'-','-',':',':','--','--';'k',[.7 .7 .7],'k',[.7 .7 .7],'k',[.7 .7 .7]}';
% h = [p1;p2;p3];
% set(h,r.style,vals);
% legend(h,legnom6,'location','southoutside','NumColumns',3);
% set(gca,'fontsize',12); 
% grid on

%% Emp Shares
%h5 = figure('color','w');
g.emps = g.emps - g.emps(1,:,:);

% p1 = plot(r.stime,g.emps(r.stime,:,1)); hold on;
% p2 = plot(r.stime,g.emps(r.stime,:,2)); hold on;
% p3 = plot(r.stime,g.emps(r.stime,:,3)); hold on;
% xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
% 
% h = [p1;p2;p3];
% set(h,r.style,vals);
% legend(h,legnom6,'location','southoutside','NumColumns',3);
% set(gca,'fontsize',12); 
% grid on


%% SAVE AND CLOSE
% saveas(h1,append(name,'_sir.pdf'));
% saveas(h2,append(name,'_death.pdf'));
% saveas(h3,append(name,'_gdp.pdf'));
% saveas(h4,append(name,'_earn.pdf'));
% saveas(h5,append(name,'_emp.pdf'));
% close all 
