contain_daily=[0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.75
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
0.29
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
1.00
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.50
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00
0.00];
%daily containment from March 1, 2020 thru May 15, 2021
    

%containment
dt_contain = datetime('01/Mar/2020') : datetime('15/May/2021'); 
% Put data into timetable
tt = timetable(dt_contain',contain_daily,'VariableNames',{'contain'});
% Calculate weekly containment
contain_weekly = retime(tt,'Weekly','mean');

% figure(1);
% subplot(1,2,2)
% plot(contain_weekly.Time,contain_weekly.contain);
% title('Containment (Data)');
% 
% xticks(contain_weekly.Time);
% sub1.XLim(1)='01-Mar-2020';
% sub1.XLim(2)='09-May-2021';
% xtickangle(-90);
% 
% 
% drawnow;