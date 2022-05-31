%This m file conducts a base-case simulation with various constant values of R_t

%To do so, we run the model 5 times and store the results an then plot them

%The cases are R_t = 2.8, 2.5, 2.2, 2.0, 1.8

%We take the population of the US to be 330 million

%the initial number of infected is 33 and exposed is 4 times 33

Rbase = [3.0,2.8,2.5,2.2,2.0,1.8,1.6];

t0 = 0;
tfinal = 548; 
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

for i = 1:7;

y0 = zeros(11,1);
y0(3) = 1/10000000;
y0(2) = 3*y0(3);
y0(4) = 0;
y0(1) = 1 - y0(2) - y0(3) - y0(4);
y0(6) = Rbase(i);
y0(7) = Rbase(i);
y0(5) = (y0(6)+y0(7))/2;
y0(8) = Rbase(i);
y0(9) = Rbase(i);
y0(10) = 1/10;
y0(11) = 1/10;
Rinfty = (y0(10)+y0(11))/2;

if i == 1;
[t1,y1] = ode113(@disease,[t0 tfinal],y0,opts);
cases1 = (y1(:,3)+y1(:,4));
slope = (log(cases1(50))-log(cases1(40)))./(t1(50)-t1(40));
doubling(1) = log(2)./slope;
end;

if i == 2;
[t2,y2] = ode113(@disease,[t0 tfinal],y0,opts);
cases2 = (y2(:,3)+y2(:,4));
slope = (log(cases2(50))-log(cases2(40)))./(t2(50)-t2(40));
doubling(2) = log(2)./slope;
end;

if i == 3;
[t3,y3] = ode113(@disease,[t0 tfinal],y0,opts);
cases3 = (y3(:,3)+y3(:,4));
slope = (log(cases2(50))-log(cases3(40)))./(t3(50)-t3(40));
doubling(3) = log(2)./slope;
end;

if i == 4;
[t4,y4] = ode113(@disease,[t0 tfinal],y0,opts);
cases4 = (y4(:,3)+y4(:,4));
slope = (log(cases4(50))-log(cases4(40)))./(t4(50)-t4(40));
doubling(4) = log(2)./slope;
end;

if i == 5;
[t5,y5] = ode113(@disease,[t0 tfinal],y0,opts);
cases5 = (y5(:,3)+y5(:,4));
slope = (log(cases5(50))-log(cases5(40)))./(t5(50)-t5(40));
doubling(5) = log(2)./slope;
end;

if i == 6;
[t6,y6] = ode113(@disease,[t0 tfinal],y0,opts);
cases6 = (y6(:,3)+y6(:,4));
slope = (log(cases6(50))-log(cases6(40)))./(t6(50)-t6(40));
doubling(6) = log(2)./slope;
end;

if i == 7;
[t7,y7] = ode113(@disease,[t0 tfinal],y0,opts);
cases7 = (y7(:,3)+y7(:,4));
slope = (log(cases7(50))-log(cases7(40)))./(t7(50)-t7(40));
doubling(7) = log(2)./slope;
end;

end;

%Now we plot cumulative cases relative to the population

plot(t1,cases1,t2,cases2,t3,cases3,t4,cases4,t5,cases5,t6,cases6,t7,cases7)
title('The ratio of cumulative cases to the population','interpreter','latex')
xlabel('days')
ylabel('cumulative cases over population','interpreter','latex')
legend({'Rt = 3.0','Rt = 2.8','Rt = 2.5','Rt = 2.2','Rt = 2.0','Rt = 1.8','Rt = 1.6'},'Location','northwest')
print('Caseslongscenario1','-dpdf','-bestfit');

N = 330000000;

plot(t1(1:80),cases1(1:80)*N,t2(1:80),cases2(1:80)*N,t3(1:80),cases3(1:80)*N,t4(1:80),cases4(1:80)*N,t5(1:80),cases5(1:80)*N,t6(1:80),cases6(1:80)*N,t7(1:80),cases7(1:80)*N)
title('The number of cumulative cases for the first 2+ months','interpreter','latex')
xlabel('days')
ylabel('cumulative cases','interpreter','latex')
legend({'Rt = 3.0','Rt = 2.8','Rt = 2.5','Rt = 2.2','Rt = 2.0','Rt = 1.8','Rt = 1.6'},'Location','northwest')
print('Casesshortscenario1','-dpdf','-bestfit');

plot(t1,y1(:,3),t2,y2(:,3),t3,y3(:,3),t4,y4(:,3),t5,y5(:,3),t6,y6(:,3),t7,y7(:,3))
title('The fraction of the population currently infected over time','interpreter','latex')
xlabel('days')
ylabel('Fraction infected','interpreter','latex')
legend({'Rt = 3.0','Rt = 2.8','Rt = 2.5','Rt = 2.2','Rt = 2.0','Rt = 1.8','Rt = 1.6'},'Location','northwest')
print('Infectedscenario1','-dpdf','-bestfit');



