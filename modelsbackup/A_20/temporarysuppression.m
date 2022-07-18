%This m file implements a simulation with a temporary intense suppression
%of transmission of the virus In this simulation, we assume that extreme
%suppression efforts get R_t well below 1 for over 1 month, but then these
%are gradually relaxed. We run the simulation for a year and a half

y0 = zeros(11,1);
y0(3) = 1/10000000;
y0(2) = 3*y0(3);
y0(4) = 0;
y0(1) = 1 - y0(2) - y0(3) - y0(4);
y0(6) = 10;
y0(7) = -4;
y0(5) = (y0(6)+y0(7))/2;
y0(8) = -4;
y0(9) = 10;
y0(10) = 1/35;
y0(11) = 1/100;
Rinfty = (y0(10)+y0(11))/2;

t0 = 0;
tfinal = 548;   
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[t,y] = ode113(@disease,[t0 tfinal],y0,opts);

plot(t,y(:,5))
title('The value of $R_t$ assumed over time','interpreter','latex')
xlabel('days')
ylabel('$R_t$','interpreter','latex')
print('Rtscenario3','-dpdf','-bestfit');

plot(t,y(:,3)+y(:,4))
title('The ratio of cumulative cases to the population','interpreter','latex')
xlabel('days')
ylabel('cumulative cases over population','interpreter','latex')
print('Caseslongscenario3','-dpdf','-bestfit');

plot(t(1:150),(y(1:150,3)+y(1:150,4))*330000000)
title('The number of cumulative cases at the start','interpreter','latex')
xlabel('days')
ylabel('cumulative cases','interpreter','latex')
print('Casesshortscenario3','-dpdf','-bestfit');

plot(t,y(:,3))
title('The fraction of the population currently infected over time','interpreter','latex')
xlabel('days')
ylabel('Fraction infected','interpreter','latex')
print('infectedscenario3','-dpdf','-bestfit');

plot(t(1:150),y(1:150,3)*330000000)
title('The numbers currently infected over time','interpreter','latex')
xlabel('days')
ylabel('Numbers infected','interpreter','latex')
print('infectedshortscenario3','-dpdf','-bestfit');
