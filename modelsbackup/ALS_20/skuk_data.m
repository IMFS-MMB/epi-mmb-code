function r = skuk_data(country,filename);

cd sir_data
r.time = readmatrix("skukus123.csv",'Range','A:A','OutputType','char');
%%%%%%%%%%%%%% change by KaiLong
% due to the fact that some versions of matlab import the csv header in the
% character vector
tmp_dateheader = sum(contains(r.time,'date'));
if tmp_dateheader > 0
    r.time = r.time(2:end);
end
%%%%%%%%end of change%%%%%%%%%%%
r.time = datestr(r.time);
r.crd = readmatrix("skukus123.csv",'Range','B2');
cd ..

for i = 1:size(r.crd,1)
    if r.crd(i,1)~=1
        break
    end
end
if country==2
    r.crd(1:i-1,:) = [];
    r.time(1:i-1,:) = [];
    for i = 1:size(r.crd,1)
        if r.crd(i,1)~=2
            break
        end
    end
end
r.crd(i:end,:) = [];
r.time(i:end,:) = [];

r.tested = max(r.crd(:,6));
r.pop = r.crd(1,7);
r.ofrac = r.crd(1,8);
r.crd(:,[1 3 5:end]) = [];

% owindata starts on 31 Dec; change to jan 20 %21
r.time = r.time(21:end,:);

% release data from coronaboard.kr begins on 21 Jan; change to 20
rcum = released(country);
rcum = [0;rcum];
rcum = rcum(1:size(r.crd(21:end,1)));
r.crd = [r.crd(21:end,1) rcum r.crd(21:end,2)];

% format date
r.time(:,3) = blanks(size(r.time,1));
for j = 1:length(r.time)
    if r.time(j,1)=="0"
        r.time(j,1) = blanks(1);
    end
end

%% PLOTS
if 0
    %h1 = figure('color','w');
    noms = {'linewidth','linestyle','color'};
    left_color = [0 0 0];
    right_color = [1 1 1];
    set(h1,'defaultAxesColorOrder',[left_color; right_color]);

    %p1 = semilogy(r.crd(:,1:2)); hold on;
    % vals = {3,2;':','--';'k','b'}';
    vals = {3,2;':','--';[.4 .4 .4],[.7 .7 .7]}';
    set(p1,noms,vals)
    ylim([1 Inf])

    yyaxis right
    %p2 = semilogy(r.crd(:,3)); hold on;
    % vals = {3;'-';'r'}';
    vals = {3;'-';'k'}';
    set(p2,noms,vals)

    xlim([1 Inf])
    xticks([1:60:j]);
    xticklabels(r.time(1:60:end,1:6));
    legend({'Confirmed','Released','Deaths (right)'},'NumColumns',3,'location','southoutside','fontsize',12);
    set(gca,'fontsize',12);
    grid on

    %% SAVE AND CLOSE
    saveas(h1,filename);
    close(h1);
end

r.xaxis = 30 + [1:size(r.time,1)]';
%r.crd(1:30,2);
end
