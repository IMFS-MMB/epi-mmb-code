%plotting
ia=4;
ib=3;
fsize=12;
horz=100;

time=0:1:horz-1;

figure;

subplot(ia,ib,1)
plot(time(1:end-1),0*100*(y(2:horz)-y_ss)/y_ss,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(y(2:horz)-y_ss)/y_ss,'b-','LineWidth',2);hold off
box off;
title('GDP, Y','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,2)
plot(time(1:end-1),0*100*(c(2:horz)-c_ss)/c_ss,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(c(2:horz)-c_ss)/c_ss,'b-','LineWidth',2);hold off
box off;
title('Consumption, C','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time(1:end-1),0*100*(x(2:horz)-x_ss)/x_ss,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(x(2:horz)-x_ss)/x_ss,'b-','LineWidth',2);hold off
box off;
title('Investment, X','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,4)
plot(time,0*100*(k(1:horz)-k_ss)/k_ss,'m:','LineWidth',1.5);hold on
plot(time,100*(k(1:horz)-k_ss)/k_ss,'b-','LineWidth',2);hold off
box off;
title('Capital, K','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time(1:end-1),0*100*(c(2:horz)-c_ss)/c_ss,'m:','LineWidth',1.5);hold on
plot(time(1:end-1),100*(n(2:horz)-n_ss)/n_ss,'b-','LineWidth',2);hold off
box off;
title('Hours, N','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,6)
plot(time(1:end-1),100*((rr(2:horz)).^52-1),'b-','LineWidth',2);hold off
box off;
title('Real Interest Rate, rr','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,7)
plot(time,100*i(1:horz),'b-','LineWidth',2);hold off
box off;
title('Infected, I','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,8)
plot(time,100*s(1:horz),'b-','LineWidth',2);hold off
box off;
title('Susceptibles, S','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,9)
plot(time,100*r(1:horz),'b-','LineWidth',2);hold off
box off;
title('Recovered, R','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,10)
plot(time,100*dd(1:horz),'b-','LineWidth',2);hold off
box off;
title('Deaths, D','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,11)
plot(time(1:end-1),100*(Rb(2:horz).^52-1),'b-','LineWidth',2);hold off
box off;
title('Policy Rate, Rb','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,12)
plot(time(1:end-1),100*(pie(2:horz).^52-1),'b-','LineWidth',2);hold off
box off;
title('Inflation Rate, pie','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);


suptitle('Figure 1: Simulation Results');
orient landscape
print -dpdf -fillpage fig1


 
