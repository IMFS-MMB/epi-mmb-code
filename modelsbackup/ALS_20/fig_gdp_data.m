h3 = figure('color','w');

r.style = {'linewidth','linestyle','color'};
sk=[0	2.391208835	-1.456913459	-7.443068274	-7.728551219	-4.997261516	-2.140111654	-1.383797859	-4.611778098	-4.348492443
    0	2.93565184	0.407590743	2.327795717	-1.948810863	-4.692745864	-0.136287708	-0.244227055	-0.929699791	1.865508595
0	2.497687327	-1.017576318	-1.110083256	-3.977798335	-5.087881591	-1.202590194	-1.110083256	-1.942645698	0.277520814];
uk=[0	0.436092177	0.190897045	-8.908301119	-29.85404401	-27.12308932	-18.95993553	-12.89750297	-10.50019124
    0	-0.054985166	-0.265566651	-5.048574129	-19.06838751	-18.05151291	-12.77083119	-8.504684268	-7.328469681
0	0.295857988	0	-7.297830375	-25.34516765	-23.37278107	-16.37080868	-11.04536489	-9.171597633];

% sk=sk/100;
% uk=uk/100;

X=sk;
p1 = plot(0:9,X); hold on;

% vals = {2,2,2;'--',':','-';'b','r','k'}';
vals = {3,3,3;'--',':','-';[.4 .4 .4],[.7 .7 .7],'k'}';
set(p1,r.style,vals)

xlim([0 9]);
xticks(0:1:9);
xticklabels({'2019','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'}); 
ylim([-30 5]);

legend(p1,{'Low-Skill','High-Skill','Total'},'location','southoutside','NumColumns',3);
set(gca,'fontsize',12); 
grid on
saveas(h3,'../data/sk_gdp_data.pdf');

h4 = figure('color','w');

X=uk;
p1 = plot(0:8,X); hold on;

% vals = {2,2,2;'--',':','-';'b','r','k'}';
vals = {3,3,3;'--',':','-';[.4 .4 .4],[.7 .7 .7],'k'}';
set(p1,r.style,vals)

xlim([0 8]);
xticks(0:1:8);
xticklabels({'2019','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'}); 
ylim([-30 5]);

legend(p1,{'Low-Skill','High-Skill','Total'},'location','southoutside','NumColumns',3);
set(gca,'fontsize',12); 
grid on
saveas(h4,'../data/uk_gdp_data.pdf');

close all
