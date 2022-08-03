%plotting
ia=2;
ib=3;
fsize=12;
horz=100;
time=0:1:horz-1;


figure;
subplot(1,2,1)
plot(time(1:end-1),100*(cs(2:horz)-c_ss)/c_ss,'r--','LineWidth',2); hold on
plot(time(1:end-1),100*(ci(2:horz)-c_ss)/c_ss,'k-.','LineWidth',2); hold on
plot(time(1:end-1),100*(cr(2:horz)-c_ss)/c_ss,'m:','LineWidth',2); hold on
plot(time(1:end-1),0*time(1:end-1),'b-','LineWidth',1); hold off
box off;
title('Consumption by Type','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Susceptibles','Infected','Recovered','Location','SouthEast');
legend boxoff

 
subplot(1,2,2)
plot(time(1:end-1),100*(ns(2:horz)-n_ss)/n_ss,'r--','LineWidth',2); hold on
plot(time(1:end-1),100*(ni(2:horz)-n_ss)/n_ss,'k-.','LineWidth',2); hold on
plot(time(1:end-1),100*(nr(2:horz)-n_ss)/n_ss,'m:','LineWidth',2); hold on
plot(time(1:end-1),0*time(1:end-1),'b-','LineWidth',1); hold off
box off;
title('Hours by Type','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Susceptibles','Infected','Recovered','Location','best');
legend boxoff
 
suptitle('Figure 2: Simulation Results');
orient landscape
print -dpdf -fillpage fig2


 
