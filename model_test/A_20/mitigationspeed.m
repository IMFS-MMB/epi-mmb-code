%This m file conducts a simulations of mitigation speed maintained forever with various values
%of \eta

%To do so, we run the model 5 times and store the results an then plot them

%The cases are

%We take the population of the US to be 330 million

%the initial number of infected is 33 and exposed is 4 times 33

etacase = [1/5,1/10,1/20,1/50,1/100];

t0 = 0;
tfinal = 548;  

for i = 1:5;

y0 = zeros(11,1);
y0(3) = 1/10000000;
y0(2) = 3*y0(3);
y0(4) = 0;
y0(1) = 1 - y0(2) - y0(3) - y0(4);
y0(6) = 3.0;
y0(7) = 3.0;
y0(5) = (y0(6)+y0(7))/2;
y0(8) = 1.6;
y0(9) = 1.6;
y0(10) = etacase(i);
y0(11) = etacase(i);
Rinfty = (y0(10)+y0(11))/2;

if i == 1;
[t1,y1] = ode45(@disease,[t0 tfinal],y0);
cases1 = (y1(:,3)+y1(:,4));
slope = (log(cases1(12))-log(cases1(1)))./(t1(12)-t1(1));
doubling(1) = log(2)./slope;
end;

if i == 2;
[t2,y2] = ode45(@disease,[t0 tfinal],y0);
cases2 = (y2(:,3)+y2(:,4));
slope = (log(cases2(12))-log(cases2(1)))./(t2(12)-t2(1));
doubling(2) = log(2)./slope;
end;

if i == 3;
[t3,y3] = ode45(@disease,[t0 tfinal],y0);
cases3 = (y3(:,3)+y3(:,4));
slope = (log(cases2(12))-log(cases3(1)))./(t3(12)-t3(1));
doubling(3) = log(2)./slope;
end;

if i == 4;
[t4,y4] = ode45(@disease,[t0 tfinal],y0);
cases4 = (y4(:,3)+y4(:,4));
slope = (log(cases4(12))-log(cases4(1)))./(t4(12)-t4(1));
doubling(4) = log(2)./slope;
end;

if i == 5;
[t5,y5] = ode45(@disease,[t0 tfinal],y0);
cases5 = (y5(:,3)+y5(:,4));
slope = (log(cases5(12))-log(cases5(1)))./(t5(12)-t5(1));
doubling(5) = log(2)./slope;
end;

end;

%First we plot the path of R_t

plot(t1,y1(:,5),t2,y2(:,5),t3,y3(:,5),t4,y4(:,5),t5,y5(:,5))
title('The value of $R_t$ assumed over time','interpreter','latex')
xlabel('days')
ylabel('$R_t$','interpreter','latex')
legend({'very fast','fast','moderate','slow','very slow'},'Location','northeast')
print('Rtscenario2','-dpdf','-bestfit');



%Now we plot cumulative cases relative to the population

plot(t1,cases1,t2,cases2,t3,cases3,t4,cases4,t5,cases5)
title('The ratio of cumulative cases to the population','interpreter','latex')
xlabel('days')
ylabel('cumulative cases over population','interpreter','latex')
legend({'very fast','fast','moderate','slow','very slow'},'Location','northwest')
print('Caseslongscenario2','-dpdf','-bestfit');

N = 330000000;

plot(t1(1:33),cases1(1:33)*N,t2(1:25),cases2(1:25)*N,t3(1:25),cases3(1:25)*N,t4(1:25),cases4(1:25)*N,t5(1:25),cases5(1:25)*N)
title('The number of cumulative cases for the first 2 months','interpreter','latex')
xlabel('days')
ylabel('cumulative cases','interpreter','latex')
legend({'very fast','fast','moderate','slow','very slow'},'Location','northwest')
print('Casesshortscenario2','-dpdf','-bestfit');

plot(t1,y1(:,3),t2,y2(:,3),t3,y3(:,3),t4,y4(:,3),t5,y5(:,3))
title('The fraction of the population currently infected over time','interpreter','latex')
xlabel('days')
ylabel('Fraction infected','interpreter','latex')
legend({'very fast','fast','moderate','slow','very slow'},'Location','northwest')
print('Infectedscenario2','-dpdf','-bestfit');