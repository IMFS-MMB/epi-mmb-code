%clear; clc;close all

%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE TO SELF:
% 0/ RUN PARAM AND SOLVESAFE IF CHANGING SS
% 1/ PERHAPS MODIFY HOME-CHOICE PROBABILITY WHEN MAKING OCC CHOICE
%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------
%% COUNTRY SPECIFIC SHIT
%---------------------------
% demographics: input from sk and uk separately
%---------------------------
folder = 'figures_uk/';

%https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/bulletins/annualmidyearpopulationestimates/latest#ageing
%0-14: 11960981;
p.pop = 43155239;
p.o0 = 11680587;
p.o0 = p.o0/(p.pop-p.o0);

%---------------------------
% technology: input from sk and uk separately
%---------------------------
p.L0 = [2046.07,	412.445+1590.73,	14022.76; ...
		1825.955,	257.389+1132.154,	9982.274];
p.L0 = p.L0/sum(p.L0,'all');

% low(m,w); high(m,w)
wages = [	4464.300	2012.320;
			5181.220	2660.520];

share = [sum(p.L0(:,1:2),2),p.L0(:,3)];
share = wages.*share;
share = sum(share,2);

p.thet = share(2)/share(1)+1;
p.thet = 1/p.thet;

share = wages.*p.L0(:,2:3);
p.alph = share(:,2)./share(:,1)+1;
p.alph = 1./p.alph;

%---------------------------
%% SIR params
%---------------------------
gpara = 1/14;

% toggle to match initial rise and death rate by occ
vpara = gpara * 3.9;	%2.177;

% cfr of old
mpara = 0.152495752154401;

%------------
% get sick params
%------------
common_params;

%------------
%% hetero params
%------------
p.yv0 = [	0.040725529193878	-0.520682871341705	0.166166186332703;
			-0.335949033498764	-0.647982895374298	-0.065895065665245];
p.psi0=[	0.168915167450905	0.153531476855278	0.077645488083363;
			0.149883821606636	0.177869692444801	0.121577709913254];
p.lockj = [	0.736981272697449	0.74560022354126	0.730324387550354;
			0.959773659706116	0.813492894172668	0.781773507595062];
gdpd = [74.6876437134882;80.9316124926297];

% force of infection
% b1+b2 by sector to 3:4 death rate (scotland)
% b1:b2 to match death rate by occ
%p.yb = zeros(p.Ns,p.No-1,p.No-1);
%p.yb(1,:,:) = [0.4,0.2;0.2,0.4];
%p.yb(2,:,:) = [0.5,0.1;0.1,0.5];

p.yv = p.yv0 - sum(p.yv0.*p.L0,'all');
p.yv = -p.yv * vpara/min(p.yv,[],'all');
p.yv = p.yv + vpara;

p.yq = vpara/7;
p.oq = vpara;

% wfh prod
gdp0 = [sum(p.L0(:,1:2),2),p.L0(:,3)];
gdp0= gdp0.*wages;
gdp0 = sum(gdp0,2);

shift = [0.387424095277876;0.267806503904904];
%shift = fsolve(@(x) f.getpsi(x,p,wages,gdp0,gdpd),shift, ...
%		optimset('display','off','MaxIter',maxit,'MaxFunEvals',maxit ...
%		));
p.psi = shift.*p.psi0 + (1.-shift).*ones(p.Ns,p.No);

% switch logit extreme values
% only need to run once!
p.qloc = 0.0345355698720605;
p.loc = [	0	-0.0141446806644398	2.14423696009744;
			0	-0.265715026850898	1.91157301210788];
runsafe = 0;
common_safe;
if runsafe==1
	format long g
	p.qloc
	p.loc
	format short
	return
end
%------------
%------------

%% base: 1 infection, do nothing
% initialize
Ia = Ia0;
base = safe;
base(1).y0(:,:,1,:) = (1d0-Ia)*safe(1).y0(:,:,1,:);
base(1).y0(:,:,2,:) = Ia*safe(1).y0(:,:,1,:);
base(1).o0(1,:) = (1d0-Ia)*safe(1).o0(1,:);
base(1).o0(2,:) = Ia*safe(1).o0(1,:);
base(2) = base(1);

% policies
for t = 1:p.Nt
	%t0 = 0.0001; %the shock is here
	t0 = I0in;
    t1 = 0.03;
	t2 = 0.3;
	q1 = 0.55;
    if (t<42)
        base(t).tau = p.omeg * [0.0;0.0];
    elseif (t<51)	% Feb 1	first detection	(Jan 31 in Korea)
    	base(t).tau = p.omeg * [0.00;t0];
    elseif (t<65)	% Feb 10 first quarantine
    	base(t).tau = p.omeg * [0.00;t0];
    elseif (t<93)	% Feb 24 testing restarts
    	base(t).tau = p.omeg * [0.00;t0+(t1-t0)*(t-64)/(92-64)];
%    elseif (t<124)	% Mar 23 lockdown
%    	base(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];
    else	% Apr 23 a-testing starts / May 30: mass track and trace in place
        base(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];%t0*(t-123)/(161-123)
    end
end
base = solvemodel(base,Ia,p);

if 0 %commented out for baseline only, by KaiLong
    %----------------------
    %% CALIBRATE
    % INITIAL PEAK
    % DEATHS BY SECTOR
    % DEATHS BY OCCUPATION
    % TESTING RATES
    %----------------------
    %% lockdown
    %----------------------
    % initialize
    Ia = Ia0;
    lock = base;

    % policies: 4.23: start asymp 6.1: reach max.
    %https://en.wikipedia.org/wiki/Timeline_of_the_COVID-19_pandemic_in_the_United_Kingdom_(January%E2%80%93June_2020)
    % Mar 15 lockdown imminent (t=85)
    for t = 1:p.Nt
        lock(t).j = f.mit(t,93,93+ld0,ld1);
        t0 = 0.0001;
        t1 = 0.03;
        t2 = 0.3;
        q1 = 0.55;
        if (t<42)
            lock(t).tau = p.omeg * [0.0;0.0];
        elseif (t<51)	% Feb 1	first detection	(Jan 31 in Korea)
            lock(t).tau = p.omeg * [0.00;t0];
        elseif (t<65)	% Feb 10 first quarantine
            lock(t).tau = p.omeg * [0.00;t0];
            lock(t).Q = 1;
        elseif (t<93)	% Feb 24 testing restarts
            lock(t).tau = p.omeg * [0.00;t0+(t1-t0)*(t-64)/(92-64)];
            lock(t).Q = 1;
        else	% Apr 23 a-testing starts (124) / May 30: mass track and trace in place
            lock(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];%t0*(t-123)/(161-123)
            lock(t).Q = 2;
            lock(t).q = q1;
        end
    end
    lock = solvemodel(lock,Ia,p);

    %% high tests and track
    % initialize
    Ia = Ia0;
    track = base;

    % policies
    for t = 1:p.Nt
        t0 = 0.03;
        t1 = 0.8;
        q0 = 0.1;
        q1 = 0.94;
        q2 = 0.61;
        q3 = 0.9;
        q4 = 0.78;
        if (t<30)
            track(t).tau = p.omeg * [0.0;0.0];
        elseif (t<61)	% Jan 20 first detection (21 in Korea)
            track(t).tau = p.omeg * t0*ones(2,1);
            track(t).Q = 1;
            track(t).q = q0;
        elseif (t<117)	% Feb 20 SCJ outbreak (19 in Korea)
            track(t).tau = p.omeg * ones(2,1)* (t0+(t1-t0)*(t-60)/(117-60));
            track(t).Q = 1;
            track(t).q = q1;
    % Apr 18 (117) lifting of restrictions/mass testing.
    % May 1 (130) Nightclub incident.
        elseif (t<236)
            track(t).tau = p.omeg * ones(2,1)* t1;
            track(t).Q = 1;
            track(t).q = q2 + f.mit(t,117,236,3)* (q1-q2);
    % Aug 15 (236) Seoul restrictions.
        elseif (t<265)
            track(t).tau = p.omeg * ones(2,1)* t1;
            track(t).Q = 1;
            track(t).q = q3 + f.mit(t,236,265,2)* (q2-q3);
    % Sep 13 (265) Seoul eases restrictions.
        else
            track(t).tau = p.omeg * ones(2,1)* t1;
            track(t).Q = 1;
            track(t).q = q4 + f.mit(t,265,265+58,2)* (q3-q4);
    % Oct 12 lifting again (Korea Oct 12).
        end
    end
    track = solvemodel(track,Ia,p);

end
%% PLOTS 
% data
r = skuk_data(2,append(folder,'uk_data.pdf'));
r.style = {'linewidth','linestyle','color'};
r.stime = [1:p.Nt];
if p.Nt==365
    r.xtnum = [31:60:p.Nt];
    r.xtlab = {'Jan 21','Mar 21','May 20','Jul 19','Sep 17','Nov 16'};
else
    r.xtnum = [31:120:p.Nt];
    r.xtlab = {'Mar 21','Jul 19','Nov 16','Mar 16','Jul 14','Nov 11'};
end    
legnom = {'Model','No policy','Track'};
legnom2 = {'Model','Early','No policy','Long','Track'};
vals = {3,3,3;'-',':','--';'k',[.7 .7 .7],[.4 .4 .4]}';
vals2 = {3,3,3,3,3;'-',':','--','--',':';'k',[.8 .8 .8],[.5 .5 .5],'k',[.2 .2 .2]}';
% 
baseg = plots(base,p,r,append(folder,'uk_base'));
%lockg = plots(lock,p,r,append(folder,'uk_lock'));
%trackg = plots(track,p,r,append(folder,'uk_track'));

% by KaiLong, extract results
uk_sir(:,1) = baseg.SIR(:,2);%-baseg.SIR(:,3);
uk_death(:,1) = baseg.D01(:,1);
uk_gdp(:,1) = baseg.prod(:,3);
%uk_emp(:,1) = baseg.emp;


% uk_sir(:,1) = lockg.SIR(:,2);%-lockg.SIR(:,3);
% uk_sir(:,2) = baseg.SIR(:,2);%-baseg.SIR(:,3);
% uk_sir(:,3) = trackg.SIR(:,2);%-trackg.SIR(:,3);

% uk_death(:,1) = lockg.D01(:,1);
%uk_death(:,2) = baseg.D01(:,1);
% uk_death(:,3) = trackg.D01(:,1);

% uk_gdp(:,1) = lockg.prod(:,3);
%uk_gdp(:,2) = baseg.prod(:,3);
% uk_gdp(:,3) = trackg.prod(:,3);

% %-------------------    
% if 0   %commented out for baseline only, by KaiLong
%     h2 = figure('color','w');
%     p2 = semilogy(r.stime,uk_death*p.pop); hold on;
%     set(p2,r.style,vals)
%     legend(p2,legnom,'location','southoutside','NumColumns',3);
%     xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
%     ylim([1 Inf]);
%     set(gca,'fontsize',12');     
%     grid on
%     saveas(h2,append(folder,'uk_cf_death.pdf'));
% end
% %-------------------    
% if 0 %commented out for baseline only, by KaiLong
%     h3 = figure('color','w');
%     p3 = plot(r.stime,uk_gdp); hold on;
%     set(p3,r.style,vals)
%     legend(p3,legnom,'location','southoutside','NumColumns',3);
%     xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
%     ylim([-Inf, 0]);
%     set(gca,'fontsize',12');     
%     grid on
%     saveas(h3,append(folder,'uk_cf_gdp.pdf'));
% 
%     close all
%     addpath('.');
% 
%     cd figures_uk
%     %cd figures_nocc
%     %cd figures_nome
% 
%     !find -name "*base*.pdf" -exec pdfcrop {} \;
%     !find -name "*lock*.pdf" -exec pdfcrop {} \;
%     !find -name "*track*.pdf" -exec pdfcrop {} \;
%     !find -name "*data*.pdf" -exec pdfcrop {} \;
%     !find -name "*cf_*.pdf" -exec pdfcrop {} \;
%     !rename -f 's/-crop.pdf/.pdf/' *-crop.pdf;
%     %!rename -f 's/base/uk_base/' base*;
%     %!rename -f 's/lock/uk_lock/' lock*;
%     %!rename -f 's/track/uk_track/' track*;
%     cd ..
% 
%     %----------------------
%     %% early lockdown
%     %----------------------
%     % initialize
%     Ia = Ia0;
%     early = base;
% 
%     % policies: 4.23: start asymp 6.1: reach max.
%     %https://en.wikipedia.org/wiki/Timeline_of_the_COVID-19_pandemic_in_the_United_Kingdom_(January%E2%80%93June_2020)
%     % Mar 15 lockdown imminent (t=85)
%     for t = 1:p.Nt
%         early(t).j = f.mit(t,65,65+ld0,ld1);
%         t0 = 0.0001;
%         t1 = 0.03;
%         t2 = 0.3;
%         q1 = 0.55;
%         if (t<42)
%             early(t).tau = p.omeg * [0.0;0.0];
%         elseif (t<51)	% Feb 1	first detection	(Jan 31 in Korea)
%             early(t).tau = p.omeg * [0.00;t0];
%         elseif (t<65)	% Feb 10 first quarantine
%             early(t).tau = p.omeg * [0.00;t0];
%             early(t).Q = 1;
%         elseif (t<93)	% Feb 24 testing restarts
%             early(t).tau = p.omeg * [0.00;t0+(t1-t0)*(t-64)/(92-64)];
%             early(t).Q = 2;
%             early(t).q = q1;%	q1*(t-92)/(124-92);
%     %    elseif (t<124)	% Mar 23 lockdown
%     %    	early(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];
%     %        early(t).Q = 2;
%     %        early(t).q = q1;%	q1*(t-92)/(124-92);
%         else	% Apr 23 a-testing starts / May 30: mass track and trace in place
%             early(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];%t0*(t-123)/(161-123)
%             early(t).Q = 2;
%             early(t).q = q1;
%         end
%     end
%     early = solvemodel(early,Ia,p);
%     % 
% 
%     %----------------------
%     %% long lockdown
%     %----------------------
%     % initialize
%     Ia = Ia0;
%     long = base;
% 
%     % policies: 4.23: start asymp 6.1: reach max.
%     %https://en.wikipedia.org/wiki/Timeline_of_the_COVID-19_pandemic_in_the_United_Kingdom_(January%E2%80%93June_2020)
%     % Mar 15 lockdown imminent (t=85)
%     for t = 1:p.Nt
%         long(t).j = f.mit(t,93,93+360,ld1);
%         t0 = 0.0001;
%         t1 = 0.03;
%         t2 = 0.3;
%         q1 = 0.55;
%         if (t<42)
%             long(t).tau = p.omeg * [0.0;0.0];
%         elseif (t<51)	% Feb 1	first detection	(Jan 31 in Korea)
%             long(t).tau = p.omeg * [0.00;t0];
%         elseif (t<65)	% Feb 10 first quarantine
%             long(t).tau = p.omeg * [0.00;t0];
%             long(t).Q = 1;
%         elseif (t<93)	% Feb 24 testing restarts
%             long(t).tau = p.omeg * [0.00;t0+(t1-t0)*(t-64)/(92-64)];
%             long(t).Q = 1;
%     %    elseif (t<124)	% Mar 23 lockdown
%     %    	long(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];
%     %        long(t).Q = 2;
%     %        long(t).q = q1;%	q1*(t-92)/(124-92);
%         else	% Apr 23 a-testing starts / May 30: mass track and trace in place
%             long(t).tau = p.omeg * [0.00;t1+(t2-t1)*(t-92)/(161-92)];%t0*(t-123)/(161-123)
%             long(t).Q = 2;
%             long(t).q = q1;
%         end
%     end
%     long = solvemodel(long,Ia,p);% 
% 
%     earlyg = plots(early,p,r,append(folder,'uk_early'));
%     longg = plots(long,p,r,append(folder,'uk_long'));
% 
%     uk_sir(:,4) = longg.SIR(:,2);%-trackg.SIR(:,3);
%     uk_sir(:,5) = uk_sir(:,3);
%     uk_sir(:,3) = uk_sir(:,2);
%     uk_sir(:,2) = earlyg.SIR(:,2);%-baseg.SIR(:,3);
% 
%     uk_death(:,4) = longg.D01(:,1);
%     uk_death(:,5) = uk_death(:,3);
%     uk_death(:,3) = uk_death(:,2);
%     uk_death(:,2) = earlyg.D01(:,1);
% 
%     uk_gdp(:,4) = longg.prod(:,3);
%     uk_gdp(:,5) = uk_gdp(:,3);
%     uk_gdp(:,3) = uk_gdp(:,2);
%     uk_gdp(:,2) = earlyg.prod(:,3);
% 
%     %-------------------    
%     h2 = figure('color','w');
%     p2 = semilogy(r.stime,uk_death*p.pop); hold on;
%     set(p2,r.style,vals2)
%     legend(p2,legnom2,'location','southoutside','NumColumns',3);
%     xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
%     ylim([1 Inf]);
%     set(gca,'fontsize',12');     
%     grid on
%     saveas(h2,append(folder,'uk_cf1_death.pdf'));
% 
%     %-------------------    
%     h3 = figure('color','w');
%     p3 = plot(r.stime,uk_gdp); hold on;
%     set(p3,r.style,vals2)
%     legend(p3,legnom2,'location','southoutside','NumColumns',3);
%     xlim([1 p.Nt]);xticks(r.xtnum);xticklabels(r.xtlab); 
%     ylim([-Inf, 0]);
%     set(gca,'fontsize',12');     
%     grid on
%     saveas(h3,append(folder,'uk_cf1_gdp.pdf'));
% 
%     close all
%     addpath('.');
% 
%     cd figures_uk
%     %cd figures_nocc
%     %cd figures_nome
% 
%     !find -name "*early*.pdf" -exec pdfcrop {} \;
%     !find -name "*long*.pdf" -exec pdfcrop {} \;
%     !find -name "*cf1*.pdf" -exec pdfcrop {} \;
%     !rename -f 's/-crop.pdf/.pdf/' *-crop.pdf;
%     %!rename -f 's/early/uk_early/' early*;
%     %!rename -f 's/long/uk_long/' long*;
%     cd ..
% 
%     nmax = r.xaxis(end);
%     uk_tab(1,:) = floor(uk_death(nmax,:)*p.pop+0.5);
%     uk_tab(2,:) = sum(uk_gdp(31:nmax,:),1)/(nmax-30);
% 
%     save uk_results
% 
%     return
% end

