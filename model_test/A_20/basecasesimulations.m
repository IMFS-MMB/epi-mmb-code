%This m file conducts a base-case simulation with various constant values of R_t

%To do so, we run the model 5 times and store the results an then plot them

%The cases are R_t = 2.8, 2.5, 2.2, 2.0, 1.8

%We take the population of the US to be 330 million

%the initial number of infected is 33 and exposed is 4 times 33

Rbase = [3.0,2.8,2.5,2.2,2.0,1.8,1.6];

t0 = 0;
tfinal = 548; 
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

helper=load('inf_ini.mat');

y0 = zeros(11,1);
%y0(3) = 1/10000000;
y0(3)=helper.helper;
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


[t,y] = ode113(@disease,[t0 tfinal],y0,opts);
%cases1 = (y(:,3)+y1(:,4));
%slope = (log(cases1(50))-log(cases1(40)))./(t1(50)-t1(40));
%doubling(1) = log(2)./slope;

