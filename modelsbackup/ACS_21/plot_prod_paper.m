load all_results
load Cyc_labor_prod


for i=1:4
   
    Z(i,1)=100*(Ztilde(i*10)-Ztildebar)/Ztildebar;
    
end

startDate = datenum('01-04-2020');
endDate = datenum('01-01-2021');
xData = linspace(startDate,endDate,4);


figure(1)



subplot(2,2,2)
plot(100*(Ztilde(6:42)-Ztildebar)/Ztildebar,'b-','LineWidth',2);
title('Aggregate Productivity Model');
subplot(2,2,1)
plot(xData,100*Cyc_labor_prod(294:297),'b-','LineWidth',2)
ax = gca; 
ax.XTick = xData(1:1:4);
datetick('x','QQ-YY','Keepticks')
xticks
xlim([startDate endDate]);
title('Aggregate Productivity Data');
subplot(2,2,3)
plot(100*(harm_av(6:42)-harm_avbar)/harm_avbar,'b-','LineWidth',2);
title('Aggregate Productivity - No Reallocation');
subplot(2,2,4)
plot(100*(omega(6:42)-omegabar)/omegabar,'b-','LineWidth',2);
title('Relative Size Social-Sector');
% subplot(1,3,3)
% plot(Z,'m--','LineWidth',2);hold on
% plot(100*Cyc_labor_prod(294:297),'b-','LineWidth',2);
% legend('Model','Data');
% title('Aggregate Productivity Data vs. Model');
