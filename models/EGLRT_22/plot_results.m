%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weekstr=[];
for ww=1:1:numel(Ivec)
    weekstr=[weekstr {['Week ',num2str(ww)]}];
end


%%Plot value functions etc
figure
subplot(2,2,1)
plot(grid.b,bprime.ro,'-*')
hold on
plot(grid.b,bprime.ry,'-o')
hold on
plot(grid.b,bprime.io,'-+')
hold on
plot(grid.b,bprime.iy,'-d')
hold on
plot(grid.b,bprime.so,'-s')
hold on
plot(grid.b,bprime.sy,'->')
hold on
plot(grid.b,grid.b,'b--','LineWidth',2)
axis tight;
legend({'Recovered old','Recovered young','Infected old','Infected young','Susceptible old','Susceptible young','45 degree line'},'Location','southeast');
legend boxoff;
xlabel('b');
ylabel('bprime');
title('Assets');

subplot(2,2,2)
plot(grid.b,c.ro,'-*')
hold on
plot(grid.b,c.ry,'-o')
hold on
plot(grid.b,c.io,'-+')
hold on
plot(grid.b,c.iy,'-d')
hold on
plot(grid.b,c.so,'-s')
hold on
plot(grid.b,c.sy,'->')
hold on
plot(grid.b,par.w*ones(mpar.nb,1),'b-')
text(grid.b(2),par.w+15,'w');
axis tight;
legend({'Recovered old','Recovered young','Infected old','Infected young','Susceptible old','Susceptible young'},'Location','southeast');
legend boxoff;
xlabel('b');
ylabel('c');
title('Consumption');

subplot(2,2,3)
plot(grid.b,U.ro,'-*')
hold on
plot(grid.b,U.ry,'-o')
hold on
plot(grid.b,U.io,'-+')
hold on
plot(grid.b,U.iy,'-+')
hold on
plot(grid.b,U.so,'-s')
hold on
plot(grid.b,U.sy,'->')
axis tight;
legend({'Recovered old','Recovered young','Infected old','Infected young','Susceptible old','Susceptible young'},'Location','southeast');
legend boxoff;
xlabel('b');
ylabel('U');
title('Value functions');

subplot(2,2,4)
semilogy((distF.ro(2:end)))
hold on
semilogy((distF.ry(2:end)))
hold on
semilogy((distF.io(2:end)))
hold on
semilogy((distF.iy(2:end)))
hold on
semilogy((distF.so(2:end)))
hold on
semilogy((distF.sy(2:end)))
axis tight;
title('Distance between updates of V -logscale')
legend({'Recovered old','Recovered young','Infected old','Infected young','Susceptible old','Susceptible young'})
legend boxoff;
xlabel('b');
ylabel('Dist');

suptitle('Value and Policy Functions (i=0 for all t)');
drawnow;
orient landscape;
print -dpdf -fillpage value_and_policy_functions




%simulation
% %plot consumption per capita by types and age
time=1:1:numel(muvec);
figure;
subplot(2,3,1)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.sy./cmat.sy(:,end)-1)*100);
ylabel('b');
xlabel('weeks');
zlabel('% of ini. cons.');
title('Consumption, sy')

subplot(2,3,4)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.so./cmat.so(:,end)-1)*100);
ylabel('b');
xlabel('weeks');
zlabel('% of ini. cons.');
title('Consumption, so')

subplot(2,3,2)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.iy./cmat.iy(:,end)-1)*100);
ylabel('b');
xlabel('weeks');
zlabel('% of ini. cons.');
title('Consumption, iy')

subplot(2,3,5)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.io./cmat.io(:,end)-1)*100);
xlabel('b');
ylabel('time');
zlabel('% of ini. cons.');
title('Consumption, io')

subplot(2,3,3)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.ry./cmat.ry(:,end)-1)*100);
ylabel('b');
xlabel('weeks');
zlabel('% of ini. cons.');
title('Consumption, ry')

subplot(2,3,6)
mesh(time,grid.b(mpar.grid_sim_idx),(cmat.ro./cmat.ro(:,end)-1)*100);
ylabel('b');
xlabel('weeks');
zlabel('% of ini. cons.');
title('Consumption, ro')

suptitle('Consumption by types and age in the epidemic')

orient landscape
print -dpdf -fillpage value_functions_epi_simulation






%simulation
% %plot consumption per capita by types and age (median assets)
time=1:1:numel(muvec);
figure;
subplot(2,3,1)
plot(time,(cmat.sy(2,:)./cmat.sy(2,end)-1)*100);
xlabel('weeks');
ylabel('% of ini. cons.');
title('Consumption, sy')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)


subplot(2,3,4)
plot(time,(cmat.so(2,:)./cmat.so(2,end)-1)*100);
xlabel('weeks');
ylabel('% of ini. cons.');
title('Consumption, so')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)

subplot(2,3,2)
plot(time,(cmat.iy(2,:)./cmat.iy(2,end)-1)*100);
xlabel('weeks');
ylabel('% of ini. cons.');
title('Consumption, iy')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)

subplot(2,3,5)
plot(time,(cmat.io(2,:)./cmat.io(2,end)-1)*100);
ylabel('time');
ylabel('% of ini. cons.');
title('Consumption, io')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)

subplot(2,3,3)
plot(time,(cmat.ry(2,:)./cmat.ry(2,end)-1)*100);
xlabel('weeks');
ylabel('% of ini. cons.');
title('Consumption, ry')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)

subplot(2,3,6)
plot(time,(cmat.ro(2,:)./cmat.ro(2,end)-1)*100);
xlabel('weeks');
ylabel('% of ini. cons.');
title('Consumption, ro')
axis([ 0 62 -70 10])
vline(mpar.infect_surprise_week+1)

suptitle('Consumption by types and age in the epidemic (at median assets)')

orient landscape
print -dpdf -fillpage value_functions_epi_simulation_2







%Probs of getting infected
tauymat=par.pi1*cmat.sy.*Ivec+par.pi2.*Ivec;
tauomat=par.pi1*cmat.so.*Ivec+par.pi2.*Ivec;
    
if isempty(find(tauymat>1))==0
    error('tauy > 1');
end
if isempty(find(tauomat>1))==0
    error('tauyo > 1');
end

figure;
subplot(2,2,1)
mesh(time,grid.b(mpar.grid_sim_idx),tauymat);
ylabel('b');
xlabel('weeks');
%zlabel('% of ini. cons.');
title('tau, y')
view([-23 14]);

subplot(2,2,2)
mesh(time,grid.b(mpar.grid_sim_idx),tauomat);
ylabel('b');
xlabel('weeks');
%zlabel('% of ini. cons.');
title('tau, o')
view([-23 14]); 


%constant gain learning for cfr's
subplot(2,2,3)
plot(cfryoung_constgain); hold on
plot(cfryoung_constgain*0+par.pidy);
legend('Const. Gain','Steady State');
title('Weekly Prob. of Dying, Constant Gain Learning, Young');
axis tight
vline(mpar.infect_surprise_week+1)

subplot(2,2,4)
plot(cfrold_constgain); hold on
plot(cfryoung_constgain*0+par.pido);
title('Weekly Prob. of Dying, Constant Gain Learning, Old');
legend('Const. Gain','Steady State');
axis tight
vline(mpar.infect_surprise_week+1)

orient landscape
print -dpdf -fillpage tauo_tauy_and_cons_gain_paths





%epi dynamics
figure;
subplot(2,4,1:2)
plot(time,Ivec,'b-');
title('Infections (Data), I');
xlabel('weeks');
vline(mpar.infect_surprise_week+1)

subplot(2,4,3:4)
plot(time,muvec,'b-');
title('Containment (Data), \mu');
xlabel('weeks');
vline(mpar.infect_surprise_week+1)

subplot(2,4,5)
plot(time,Iy);hold on;
plot(time,Io);hold on;
title('Infected (SIR implied)');
xlabel('weeks');
legend('young','old');
legend boxoff;
vline(mpar.infect_surprise_week+1)

subplot(2,4,6)
plot(time,Sy);hold on;
plot(time,So);hold on;
title('Susceptibles (SIR implied)');
xlabel('weeks');
legend('young','old');
legend boxoff;
vline(mpar.infect_surprise_week+1)

subplot(2,4,7)
plot(time,Ry);hold on;
plot(time,Ro);hold on;
title('Recovered (SIR implied)');
xlabel('weeks');
legend('young','old');
legend boxoff;
vline(mpar.infect_surprise_week+1)

subplot(2,4,8)
plot(time,Dy);hold on;
plot(time,Do);hold on;
title('Deceased (SIR implied)');
xlabel('weeks');
legend('young','old');
legend boxoff;
vline(mpar.infect_surprise_week+1)

orient landscape
print -dpdf -fillpage infections_containment_used_in_epi_simulation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%KEY FIGURE: EFFECTS OF EPI+CONTAINMENT ON YOUNG AND OLD
%Shades
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conflevel=0.95;
factor=-norminv((1-conflevel)/2);

figure;
sub1=subplot(1,1,1);

plot((1:1:14),tty_data.Consy_data,'b--','LineWidth',2); hold on %data mean young

grpyat = [ (1:1:14)' (tty_data.Consy_data(1:end)-factor*100*par.dat_stderr_young(1:end)'); (14:-1:1)' (tty_data.Consy_data(end:-1:1))+factor*100*par.dat_stderr_young(end:-1:1)'];
patch(grpyat(:,1),grpyat(:,2),[0.8 0.8 1],'edgecolor',[0.8 0.8 1]); hold on %data 95% young
%use uisetcolor to get RBG codes

plot((1:1:14),cons_monthly.Cy.Consy(1:end-1),'b-','LineWidth',2); hold on %model young

plot((1:1:14),tto_data.Conso_data,'r--','LineWidth',2);%data mean old

grpyat = [ (1:1:14)' (tto_data.Conso_data(1:end)-factor*100*par.dat_stderr_old(1:end)'); (14:-1:1)' (tto_data.Conso_data(end:-1:1))+factor*100*par.dat_stderr_old(end:-1:1)'];
patch(grpyat(:,1),grpyat(:,2),[1 0.8 0.8],'edgecolor',[1 0.8 0.8]); hold on%data 95% old

plot((1:1:14),cons_monthly.Co.Conso(1:end-1),'r-','LineWidth',2); hold on%model old

plot((1:1:14),tty_data.Consy_data,'b--','LineWidth',2); hold on %data mean young
plot((1:1:14),cons_monthly.Cy.Consy(1:end-1),'b-','LineWidth',2); hold on %model young
plot((1:1:14),tto_data.Conso_data,'r--','LineWidth',2);%data mean old

plot((1:1:14),0*tto_data.Conso_data,'k:','LineWidth',1.5);
titl=title('Consumption of Young and Old in the Pandemic','FontSize',12);
legend1=legend('Data: Young (Mean)','Data: Young (95%)', 'Model: Young', 'Data: Old (Mean)','Data: Old (95%)', 'Model: Old','Location','Southwest','FontSize',12);
box off;
legend box off;
ylabel('% Deviations from January 2020','FontSize',12);
set(legend1,...
    'Position',[0.317548452641451 0.250788781770377 0.264285714285714 0.229761904761905],...
    'Orientation','vertical');
sub1.XLim(1)=1;%'01-Mar-2020';
sub1.XLim(2)=14;%'01-Mar-2021';
sub1.YLim(1)=-50;
sub1.YLim(2)=15;

datstring=['Mar 2020'
    'Apr 2020'
    'May 2020'
    'Jun 2020'
    'Jul 2020'
    'Aug 2020'
    'Sep 2020'
    'Oct 2020'
    'Nov 2020'
    'Dec 2020'
    'Jan 2021'
    'Feb 2021'
    'Mar 2021'
    'Apr 2021'];

%xticks(cons_monthly.Cy.Time(1:end-1));
xticks(1:14);
%xticklabels({cons_monthly.Cy.Time(1:end-1)});
xticklabels(datstring)
xtickangle(-35);
set(gca,'FontSize',12);
set(titl,'FontSize',18);

orient landscape
print -dpdf -fillpage Figure4shades
