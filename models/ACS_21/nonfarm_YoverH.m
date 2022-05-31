clear all

load us_YoverH_nnfarm

labor_prod = table2array(nonfarmyoverh);


log_labor_prod=log(labor_prod);

% compute trend and cyclical compont of the logarithm, in this case the
% cyclical component will be in percentage

[labor_prod_trend,Cyc_labor_prod] = hpfilter(log_labor_prod,1600);



startDate = datenum('01-1-2019');
endDate = datenum('01-01-2021');
xData = linspace(startDate,endDate,9);

save Cyc_labor_prod Cyc_labor_prod;

figure(1)
plot(xData,100*Cyc_labor_prod(289:297),'b','LineWidth',3)
ax = gca; 
ax.XTick = xData(1:1:9);
datetick('x','QQ-YY','Keepticks')
xticks
xlim([startDate endDate]);
title('Labor productivity')

