load all_results;

%plotting
ia=4;
ib=4;
fsize=11;
horz=100;

time=0:1:horz-1;

figure;

subplot(ia,ib,1)
plot(time(1:end-1),0*100*(C1(2:horz)-C1bar)/C1bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(C1(2:horz)-C1bar)/C1bar,'b-','LineWidth',2);hold off
box off;
title('Consumption Social, C(s)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,2)
plot(time(1:end-1),0*100*(C(2:horz)-Cbar)/Cbar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(C(2:horz)-Cbar)/Cbar,'b-','LineWidth',2);hold off
box off;
title('Consumption, C','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time(1:end-1),0*100*(GDP(2:horz)-GDPbar)/GDPbar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(GDP(2:horz)-GDPbar)/GDPbar,'b-','LineWidth',2);hold off
box off;
title('GDP, GDP','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,4)
plot(time(1:end-1),0*100*(w(2:horz)-wbar)/wbar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(w(2:horz)-wbar)/wbar,'b-','LineWidth',2);hold off
box off;
title('Wage, w','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time(1:end-1),0*100*(no1(2:horz)-no1bar)/no1bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(no1(2:horz)-no1bar)/no1bar,'b-','LineWidth',2);hold off
box off;
title('Operative Social, No(s)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,6)
plot(time(1:end-1),0*100*(N1(2:horz)-N1bar)/N1bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(N1(2:horz)-N1bar)/N1bar,'b-','LineWidth',2);hold off
box off;
title('Firms Social, N(s)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,7)
plot(time,0*100*(Ne1(1:horz)-Ne1bar)/Ne1bar,'m:','LineWidth',1.5);hold on
plot(time,100*(Ne1(1:horz)-Ne1bar)/Ne1bar,'b-','LineWidth',2);hold off
box off;
title('Entrants Social, Ne(s)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,8)
plot(time(1:end-1),0*100*(zc1(2:horz)-zc1bar)/zc1bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(zc1(2:horz)-zc1bar)/zc1bar,'b-','LineWidth',2);hold off
box off;
title('Cut-off Social, zc(s)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,9)
plot(time(1:end-1),0*100*(no2(2:horz)-no2bar)/no2bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(no2(2:horz)-no2bar)/no2bar,'b-','LineWidth',2);hold off
box off;
title('Operative Non-Social, No(ns)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,10)
plot(time(1:end-1),0*100*(N2(2:horz)-N2bar)/N2bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(N2(2:horz)-N2bar)/N2bar,'b-','LineWidth',2);hold off
box off;
title('Firms Non-Social, N(ns)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,11)
plot(time,0*100*(Ne2(1:horz)-Ne2bar)/Ne2bar,'m:','LineWidth',1.5);hold on
plot(time,100*(Ne2(1:horz)-Ne2bar)/Ne2bar,'b-','LineWidth',2);hold off
box off;
title('Entrants Non-Social, Ne(ns)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,12)
plot(time(1:end-1),0*100*(zc2(2:horz)-zc2bar)/zc2bar,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(zc2(2:horz)-zc2bar)/zc2bar,'b-','LineWidth',2);hold off
box off;
title('Cut-off Non-Social, zc(ns)','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,13)
plot(time,328000000*ii(1:horz),'b-','LineWidth',2);hold off
box off;
title('Infected, ii','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,14)
plot(time,328000000*ss(1:horz),'b-','LineWidth',2);hold off
box off;
title('Susceptibles, ss','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,15)
plot(time,328000000*rr(1:horz),'b-','LineWidth',2);hold off
box off;
title('Recovered, rr','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,16)
plot(time,328000000*dd(1:horz),'b-','LineWidth',2);hold off
box off;
title('Deaths, dd','FontSize',fsize);
set(gca,'FontSize',fsize);

xlabel('Weeks','FontSize',fsize);

%suptitle('Figure 1: Simulation Results Benchmark');
orient landscape
print -dpdf -fillpage fig_Bench_paper


 
