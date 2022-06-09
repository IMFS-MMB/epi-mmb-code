%%% This file computes everything else: (i) baseline experiments, (ii)
%%% policy-by-policy analysis, (iii) US package multipliers
clear
clc
close all

load params.mat
load data_shocks.mat

options                = struct;
options.shock_a        = 0;
options.shock_mp       = 0;
options.shock_util_a   = 1;
options.shock_policies = 1;
    
if options.shock_policies
    options.crisis_G         = 1;
    options.crisis_taul      = 1;
    options.crisis_varsigma  = 1;
    options.crisis_transfer  = 1;
    options.crisis_govt_wage = 1;
else
    options.crisis_G         = 0;
    options.crisis_taul      = 0;
    options.crisis_varsigma  = 0;
    options.crisis_transfer  = 0;
    options.crisis_govt_wage = 0;
end

% Tables and figures in the paper
options.table_multipliers = 1;
options.bar_plots         = 1;
options.us_package        = 1;

policy_size = 0.2/(21.7/4); % size of the fiscal impulse --> USD 200bn 

% Variables to plot
vars      = {'GDP_fix'; 'Cb'; 'Cs'; 'N_a'; 'N_n'; 'w_a'; 'w_n'; 'Pi'; 'R'; 'V'; 'Fb'; 'sprqb'};
var_names = {'GDP (Fix)'; 'Cons. Borrower'; 'Cons. Saver'; 'Labor Sector A'; 'Labor Sector B'; 'Wage Sector A'; 'Wage Sector B'; 'Inflation'; 'Policy Rate'; 'Bank Stock'; 'Default Rate'; 'Credit Spread'};

vars_paper      = {'util_a'; 'unemp'; 'Cb'; 'N_a'; 'N_n'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
var_names_paper = {'Shock'; 'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};

vars_slides      = {'util_a'; 'GDP_fix'; 'R'; 'Cb'; 'N_a'; 'N_n'};
var_names_slides = {'Shock'; 'Real GDP'; 'Policy Rate'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'};
    
Tirf = 9;
time = (1:Tirf)';

% Run Dynare once
dynare model.mod noclearall;

% Indices for variables of interest
ind_unemp  = strmatch('unemp', M_.endo_names, 'exact');
ind_ua_var = strmatch('util_a', M_.endo_names, 'exact');

ind_tfp       = strmatch('ee_a', M_.exo_names, 'exact');
ind_mp        = strmatch('ee_mp', M_.exo_names, 'exact');
ind_ua        = strmatch('ee_util_a', M_.exo_names, 'exact');
ind_G         = strmatch('ee_G', M_.exo_names, 'exact');
ind_varsigma  = strmatch('ee_varsigma', M_.exo_names, 'exact');
ind_govt_wage = strmatch('ee_govt_wage', M_.exo_names, 'exact');
ind_transfer  = strmatch('ee_transfer', M_.exo_names, 'exact');
ind_tau_l     = strmatch('ee_taul', M_.exo_names, 'exact');

oo_.exo_simul(:,:) = 0; % set all shocks to zero

%% (1) TFP Shock
if options.shock_a == 1
    oo_.exo_simul(2,ind_tfp) = 0.01;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf.crisis = oo_.endo_simul(:,1:end);

    % Plots
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf.crisis(currind,1:Tirf);
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R')
            crisis = exp(irf.crisis(currind,1:Tirf));
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf.crisis(currind,1:Tirf))/detss.(currvar)-1);
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
end

%% (2) MP Shock
if options.shock_a == 1
    oo_.exo_simul(2,ind_mp) = 0.01;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf.crisis = oo_.endo_simul(:,1:end);

    % Plots
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf.crisis(currind,1:Tirf);
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R')
            crisis = exp(irf.crisis(currind,1:Tirf));
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf.crisis(currind,1:Tirf))/detss.(currvar)-1);
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
end

%% (3) Crisis shock

irf_store = struct;

if options.shock_util_a == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.nopol = oo_.endo_simul(:,1:end);
    
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            plot(time, crisis, 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            plot(time, [crisis], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    dt_start = datetime(2020,03,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;
    
    f = figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    f.Position = [962 42 958 954];
    print -depsc figures/crisis
    
    % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_slides
    
    % Inflation Figure
    ind_infl = strmatch('Pi', M_.endo_names, 'exact');
    ind_pa   = strmatch('p_a', M_.endo_names, 'exact');
    infl_data = [1.0210637^(1/4); 1.004440387^(1/4); 1.0125337^(1/4); 1.01182235^(1/4)];
    infl_data = [infl_data; NaN(Tirf-length(infl_data),1)];
    figure
    subplot(2,1,1)
    var_nopol = (1.02^(1/4))*exp(irf_store.nopol(ind_infl,1:Tirf));
    plot(dateplot, var_nopol', 'Linewidth', 3), hold on
    plot(dateplot, infl_data, '--', 'Linewidth', 3)
    title('Inflation, n-sector')
    ylabel('Gross Inflation Rate')
    axis tight
    grid minor
    set(gca, 'XTick', dateplot);
    datetick('x','QQ-YY','keeplimits','keepticks');
    legend('Model', 'Data, CPI', 'Location', 'Southeast')
    subplot(2,1,2)
    var_nopol = 100*(exp(irf_store.nopol(ind_pa,1:Tirf))/detss.p_a-1);
    plot(dateplot, var_nopol, 'Linewidth', 3)
    title('Price, a-sector')
    ylabel('% dev. from SS')
    axis tight
    grid minor
    set(gca, 'XTick', dateplot);
    datetick('x','QQ-YY','keeplimits','keepticks');
    print -depsc figures/inflation_crisis
    
    
    ind_rot = strmatch('rot', M_.endo_names, 'exact');
    rot_costs = irf_store.nopol(ind_rot,1:Tirf);
end



%% (4) Crisis shock + Govt Spending

if options.crisis_G == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(2,ind_G)    = policy_size*detss.GDP/detss.G;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.G = oo_.endo_simul(:,1:end);
    
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.G(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.G(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.G(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'unemp'; 'Cb'; 'N_a'; 'N_n'; 'Pi'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Inflation Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.G(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.G(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.G(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.G(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'G', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_pol_G
    
    % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.G(currind,1:Tirf)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.G(currind,1:Tirf))';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.G(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.G(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'G', 'Location', 'Southeast')
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_polG_slides 
    
    ind_infl = strmatch('Pi', M_.endo_names, 'exact');
    ind_pa   = strmatch('p_a', M_.endo_names, 'exact');

end

%% (5) Crisis shock + Payroll Tax Cut

if options.crisis_taul == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(2,ind_tau_l) = policy_size*detss.GDP/(detss.tau_l*(detss.w_n*detss.N_n + detss.w_a*detss.N_a));
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.tau_l = oo_.endo_simul(:,1:end);
    
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.tau_l(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.tau_l(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.tau_l(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'unemp'; 'Cb'; 'N_a'; 'N_n'; 'Pi'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Inflation Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.tau_l(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.tau_l(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.tau_l(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.tau_l(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'tau_l', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_pol_taul
    
        % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.tau_l(currind,1:Tirf)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.tau_l(currind,1:Tirf))';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.tau_l(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.tau_l(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'tau_l', 'Location', 'Southeast')
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_poltaul_slides
    
end

%% (6) Crisis shock + Unemployment Insurance

if options.crisis_varsigma == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(2,ind_varsigma) = 0.50*policy_size*detss.GDP/detss.varsigma/(1-detss.N_n-detss.N_a);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.varsigma = oo_.endo_simul(:,1:end);
    
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.varsigma(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.varsigma(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.varsigma(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'unemp'; 'Cb'; 'N_a'; 'N_n'; 'Pi'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Inflation Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.varsigma(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.varsigma(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.varsigma(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.varsigma(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'ui', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_pol_varsigma
    
        % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.varsigma(currind,1:Tirf)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.varsigma(currind,1:Tirf))';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.varsigma(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.varsigma(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'ui', 'Location', 'Southeast')
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_polui_slides
end

%% (7) Crisis shock + Transfer

if options.crisis_transfer == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(2,ind_transfer) = policy_size*detss.GDP;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.transfer = oo_.endo_simul(:,1:end);
    
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.transfer(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.transfer(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.transfer(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'unemp'; 'Cb'; 'N_a'; 'N_n'; 'Pi'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Inflation Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.transfer(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.transfer(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.transfer(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.transfer(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'Tb', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_pol_transfer
    
        % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.transfer(currind,1:Tirf)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.transfer(currind,1:Tirf))';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.transfer(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.transfer(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'Tb', 'Location', 'Southeast')
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_poltb_slides
end


%% (8) Crisis shock + Furloughed workers

if options.crisis_govt_wage == 1
    
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(2,ind_govt_wage) = policy_size*detss.GDP/detss.N_a/detss.w_a ;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.govt_wage = oo_.endo_simul(:,1:end);

    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.govt_wage(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.govt_wage(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.govt_wage(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'unemp'; 'Cb'; 'N_a'; 'N_n'; 'Pi'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Inflation Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m') || strcmp(currvar,'ci')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.govt_wage(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.govt_wage(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.govt_wage(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.govt_wage(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'Ta', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_pol_govtwage
    %     print -depsc figures/crisis_pol_govtwage_losigmaa
    
        % Figure for slides
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot = dt_start:calmonths(3):dt_end;

    f = figure;
    for ii = 1:length(vars_slides)
        currvar = vars_slides{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_slides{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.govt_wage(currind,1:Tirf)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.govt_wage(currind,1:Tirf))';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.govt_wage(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.govt_wage(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis, crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'Ta', 'Location', 'Southeast')
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_polta_slides
end

%% (9) Table with multipliers

if options.table_multipliers == 1
    
    pols = {'G'; 'tau_l'; 'varsigma'; 'transfer'; 'govt_wage'};
    pol_names = {'Govt. Cons.'; 'Income Tax'; 'UI'; 'Transfer'; 'Firm Assistance'};
    
    mult = struct;
    mult.gdp = struct;
    mult.cb  = struct;
    mult.cs  = struct;
    mult.y  = struct;
    mult.income = struct;
    
    ind_gdp   = strmatch('GDP', M_.endo_names, 'exact');
    ind_y   = strmatch('Y', M_.endo_names, 'exact');
    ind_cb    = strmatch('Cb', M_.endo_names, 'exact');
    ind_cs    = strmatch('Cs', M_.endo_names, 'exact');
    ind_R     = strmatch('R', M_.endo_names, 'exact');
    ind_spend_G = strmatch('spend_G', M_.endo_names, 'exact');
    ind_spend_tau_l = strmatch('spend_tau_l', M_.endo_names, 'exact');
    ind_spend_varsigma = strmatch('spend_varsigma', M_.endo_names, 'exact');
    ind_spend_transfer = strmatch('spend_transfer', M_.endo_names, 'exact');
    ind_spend_govt_wage = strmatch('spend_govt_wage', M_.endo_names, 'exact');
    ind_spend = [ind_spend_G; ind_spend_tau_l; ind_spend_varsigma; ind_spend_transfer; ind_spend_govt_wage];
    ind_income    = strmatch('income', M_.endo_names, 'exact');
    ind_na = strmatch('N_a', M_.endo_names, 'exact');
    ind_nn = strmatch('N_n', M_.endo_names, 'exact');
    
    horz = 20;
    
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        GDP_nopol = exp(irf_store.nopol(ind_gdp, 1:horz))';
        GDP_pol = exp(irf_store.(currpol)(ind_gdp, 1:horz))';
        
        cb_nopol = exp(irf_store.nopol(ind_cb, 1:horz))';
        cb_pol = exp(irf_store.(currpol)(ind_cb, 1:horz))';
        
        cs_nopol = exp(irf_store.nopol(ind_cs, 1:horz))';
        cs_pol = exp(irf_store.(currpol)(ind_cs, 1:horz))';
        
        y_nopol = exp(irf_store.nopol(ind_y, 1:horz))';
        y_pol = exp(irf_store.(currpol)(ind_y, 1:horz))';
        
        spend_nopol = irf_store.nopol(ind_spend(jj), 1:horz)';
        spend_pol = irf_store.(currpol)(ind_spend(jj), 1:horz)';
        
        income_nopol = irf_store.nopol(ind_income, 1:horz)';
        income_pol = irf_store.(currpol)(ind_income, 1:horz)';
        
        employ_nopol = exp(irf_store.nopol(ind_na, 1:horz))' + exp(irf_store.nopol(ind_nn, 1:horz))';
        employ_pol = exp(irf_store.(currpol)(ind_na, 1:horz))' + exp(irf_store.(currpol)(ind_nn, 1:horz))';
        
        R = [1; exp(irf_store.nopol(ind_R, 1:horz-1))'];
        cumR = 1./exp(cumsum(log(R)));
        
        mult.gdp.(currpol) = sum(cumR.*(GDP_pol-GDP_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.cb.(currpol) = sum(cumR.*(cb_pol-cb_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.cs.(currpol) = sum(cumR.*(cs_pol-cs_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.y.(currpol) = sum(cumR.*(y_pol-y_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.income.(currpol) = sum(cumR.*(income_pol-income_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.employ.(currpol) = sum(cumR.*(employ_pol-employ_nopol))./sum(cumR.*(spend_pol-spend_nopol));
    end
    
    % Table with multipliers
%     fileID = fopen('tables/multipliers.txt','w');
%     fprintf(fileID, ' $G$         & Govt. Consumption & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.G, mult.income.G, mult.cb.G, mult.cs.G, mult.gdp.G, '\\');
%     fprintf(fileID, ' $\\tau_t^l$  & Income Tax        & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.tau_l, mult.income.tau_l, mult.cb.tau_l, mult.cs.tau_l, mult.gdp.tau_l, '\\');
%     fprintf(fileID, ' $\\varsigma$ & UI                & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n',mult.employ.varsigma, mult.income.varsigma, mult.cb.varsigma, mult.cs.varsigma, mult.gdp.varsigma, '\\');
%     fprintf(fileID, ' $T_t^b$     & Uncond. Transfer  & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.transfer, mult.income.transfer, mult.cb.transfer, mult.cs.transfer, mult.gdp.transfer, '\\');
%     fprintf(fileID, ' $T_t^a$     & Liquidity Assist. & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.govt_wage, mult.income.govt_wage, mult.cb.govt_wage, mult.cs.govt_wage, mult.gdp.govt_wage, '\\');
%     fclose(fileID);
end

%% (10) Bar plot with consumption and income responses

if options.bar_plots==1
    
    pols = {'G'; 'tau_l'; 'varsigma'; 'transfer'; 'govt_wage'};
    pol_names = {'Govt. Consumption'; 'Income Tax'; 'UI'; 'Transfer'; 'Firm Assistance'};
    
    resp = struct;
    resp.gdp = struct;
    resp.cb  = struct;
    resp.cs  = struct;
    resp.y  = struct;
    resp.income = struct;
    
    ind_gdp             = strmatch('GDP', M_.endo_names, 'exact');
    ind_y               = strmatch('Y', M_.endo_names, 'exact');
    ind_cb              = strmatch('Cb', M_.endo_names, 'exact');
    ind_cs              = strmatch('Cs', M_.endo_names, 'exact');
    ind_R               = strmatch('R', M_.endo_names, 'exact');
    ind_spend_G         = strmatch('spend_G', M_.endo_names, 'exact');
    ind_spend_tau_l     = strmatch('spend_tau_l', M_.endo_names, 'exact');
    ind_spend_varsigma  = strmatch('spend_varsigma', M_.endo_names, 'exact');
    ind_spend_transfer  = strmatch('spend_transfer', M_.endo_names, 'exact');
    ind_spend_govt_wage = strmatch('spend_govt_wage', M_.endo_names, 'exact');
    ind_spend           = [ind_spend_G; ind_spend_tau_l; ind_spend_varsigma; ind_spend_transfer; ind_spend_govt_wage];
    ind_income          = strmatch('income', M_.endo_names, 'exact');
    ind_wa              = strmatch('w_a', M_.endo_names, 'exact');
    ind_na              = strmatch('N_a', M_.endo_names, 'exact');
    ind_wn              = strmatch('w_n', M_.endo_names, 'exact');
    ind_nn              = strmatch('N_n', M_.endo_names, 'exact');
    ind_vsigma          = strmatch('varsigma', M_.endo_names, 'exact');
    ind_taul            = strmatch('tau_l', M_.endo_names, 'exact');
    
    % (1) Total Change in Income
    d_income = [];
    d_spend  = [];
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        d_income_nopol = 100*(irf_store.nopol(ind_income, 2)-irf_store.nopol(ind_income, 1))/detss.income;
        d_income_pol   = 100*(irf_store.(currpol)(ind_income, 2)-irf_store.(currpol)(ind_income, 1))/detss.income;
        
        ind_spend_pol = ind_spend(jj);
        spend_name    = strcat('spend_',currpol);
        ss_pol        = detss.(spend_name);
        d_spend_nopol = 100*(irf_store.nopol(ind_spend_pol, 2)-irf_store.nopol(ind_spend_pol, 1));
        d_spend_pol   = 100*(irf_store.(currpol)(ind_spend_pol, 2)-irf_store.(currpol)(ind_spend_pol, 1));
        
        d_income = [d_income; d_income_pol-d_income_nopol];
        d_spend  = [d_spend; d_spend_pol-d_spend_nopol];
    end
    
    figure
    X = categorical(pol_names);
    X = reordercats(X,pol_names);
    y = d_income;
    b = bar(X, y, 'FaceColor', 'flat'); title('Policy Impact on Total Income'), ylabel('% Change on Impact')
    b(1).CData = [0 0 0];  
    print -depsc figures/bar_total_income
    
    % (2) Change on income per worker
    f=figure;
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        d_income_nopol = (irf_store.nopol(ind_income, 2)-irf_store.nopol(ind_income, 1))/detss.income;
        d_income_pol   = (irf_store.(currpol)(ind_income, 2)-irf_store.(currpol)(ind_income, 1))/detss.income;
        
        income_a_nopol = exp(irf_store.nopol(ind_wa, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        income_n_nopol = exp(irf_store.nopol(ind_wn, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        income_u_nopol = exp(irf_store.nopol(ind_vsigma, :));
        
        income_a_pol = exp(irf_store.(currpol)(ind_wa, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        income_n_pol = exp(irf_store.(currpol)(ind_wn, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        income_u_pol = exp(irf_store.(currpol)(ind_vsigma, :));
        
        d_income_a_nopol = (income_a_nopol(2) - income_a_nopol(1))/income_a_nopol(1);
        d_income_a_pol   = (income_a_pol(2) - income_a_pol(1))/income_a_pol(1);
        
        d_income_n_nopol = (income_n_nopol(2) - income_n_nopol(1))/income_n_nopol(1);
        d_income_n_pol   = (income_n_pol(2) - income_n_pol(1))/income_n_pol(1);
        
        d_income_u_nopol = (income_u_nopol(2) - income_u_nopol(1))/income_u_nopol(1);
        d_income_u_pol   = (income_u_pol(2) - income_u_pol(1))/income_u_pol(1);
        
        %     d_na_nopol = (exp(irf_store.nopol(ind_na, 2)) - exp(irf_store.nopol(ind_na,1)))/exp(irf_store.nopol(ind_na, 1));
        %     d_nn_nopol = (exp(irf_store.nopol(ind_nn, 2)) - exp(irf_store.nopol(ind_nn,1)))/exp(irf_store.nopol(ind_nn, 1));
        %     d_nu_nopol = ((1-exp(irf_store.nopol(ind_na, 2))-exp(irf_store.nopol(ind_nn, 2))) - (1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1))))/(1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1)));
        %
        %     d_na_pol = (exp(irf_store.(currpol)(ind_na, 2)) - exp(irf_store.(currpol)(ind_na,1)))/exp(irf_store.(currpol)(ind_na, 1));
        %     d_nn_pol = (exp(irf_store.(currpol)(ind_nn, 2)) - exp(irf_store.(currpol)(ind_nn,1)))/exp(irf_store.(currpol)(ind_nn, 1));
        %     d_nu_pol = ((1-exp(irf_store.(currpol)(ind_na, 2))-exp(irf_store.(currpol)(ind_nn, 2))) - (1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1))))/(1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1)));
        %
        X = categorical({'Sector A','Sector N','Unemployed'});
        X = reordercats(X,{'Sector A','Sector N','Unemployed'});
        y = -100*[d_income_a_nopol-d_income_a_pol; d_income_n_nopol-d_income_n_pol; d_income_u_nopol-d_income_u_pol];
        subplot(3,2,jj)
        b = bar(X, y, 'FaceColor', 'flat'); title(pol_names{jj}), ylabel('% Change')
        b(1).CData = [0 0 0];  
    end
    suptitle('Net Income per Worker')
    f.Position = [962 42 958 954];
    print -depsc figures/bar_income_per_worker
    
    % (3) Change on income per worker, E and U only
    f=figure;
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        d_income_nopol = (irf_store.nopol(ind_income, 2)-irf_store.nopol(ind_income, 1))/detss.income;
        d_income_pol   = (irf_store.(currpol)(ind_income, 2)-irf_store.(currpol)(ind_income, 1))/detss.income;
        
        income_a_nopol = exp(irf_store.nopol(ind_wa, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        income_n_nopol = exp(irf_store.nopol(ind_wn, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        income_u_nopol = exp(irf_store.nopol(ind_vsigma, :));
        
        income_a_pol = exp(irf_store.(currpol)(ind_wa, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        income_n_pol = exp(irf_store.(currpol)(ind_wn, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        income_u_pol = exp(irf_store.(currpol)(ind_vsigma, :));
        
        d_income_a_nopol = (income_a_nopol(2) - income_a_nopol(1))/income_a_nopol(1);
        d_income_a_pol   = (income_a_pol(2) - income_a_pol(1))/income_a_pol(1);
        
        d_income_n_nopol = (income_n_nopol(2) - income_n_nopol(1))/income_n_nopol(1);
        d_income_n_pol   = (income_n_pol(2) - income_n_pol(1))/income_n_pol(1);
        
        d_income_u_nopol = (income_u_nopol(2) - income_u_nopol(1))/income_u_nopol(1);
        d_income_u_pol   = (income_u_pol(2) - income_u_pol(1))/income_u_pol(1);
        
        %     d_na_nopol = (exp(irf_store.nopol(ind_na, 2)) - exp(irf_store.nopol(ind_na,1)))/exp(irf_store.nopol(ind_na, 1));
        %     d_nn_nopol = (exp(irf_store.nopol(ind_nn, 2)) - exp(irf_store.nopol(ind_nn,1)))/exp(irf_store.nopol(ind_nn, 1));
        %     d_nu_nopol = ((1-exp(irf_store.nopol(ind_na, 2))-exp(irf_store.nopol(ind_nn, 2))) - (1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1))))/(1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1)));
        %
        %     d_na_pol = (exp(irf_store.(currpol)(ind_na, 2)) - exp(irf_store.(currpol)(ind_na,1)))/exp(irf_store.(currpol)(ind_na, 1));
        %     d_nn_pol = (exp(irf_store.(currpol)(ind_nn, 2)) - exp(irf_store.(currpol)(ind_nn,1)))/exp(irf_store.(currpol)(ind_nn, 1));
        %     d_nu_pol = ((1-exp(irf_store.(currpol)(ind_na, 2))-exp(irf_store.(currpol)(ind_nn, 2))) - (1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1))))/(1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1)));
        %
        X = categorical({'Employed','Unemployed'});
        X = reordercats(X,{'Employed','Unemployed'});
        y = -100*[d_income_a_nopol-d_income_a_pol; d_income_u_nopol-d_income_u_pol];
        subplot(3,2,jj)
        b = bar(X, y, 'FaceColor', 'flat'); title(pol_names{jj}), ylabel('% Change')
        b(1).CData = [0 0 0];  
    end
    suptitle('Net Income per Worker')
    f.Position = [962 42 958 954];
    print -depsc figures/bar_income_per_worker_eu
    
    % (4) Number of workers
    f=figure;
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        d_na_nopol = (exp(irf_store.nopol(ind_na, 2)) - exp(irf_store.nopol(ind_na,1)));
        d_nn_nopol = (exp(irf_store.nopol(ind_nn, 2)) - exp(irf_store.nopol(ind_nn,1)));
        d_nu_nopol = ((1-exp(irf_store.nopol(ind_na, 2))-exp(irf_store.nopol(ind_nn, 2))) - (1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1))));
        
        d_na_pol = (exp(irf_store.(currpol)(ind_na, 2)) - exp(irf_store.(currpol)(ind_na,1)));
        d_nn_pol = (exp(irf_store.(currpol)(ind_nn, 2)) - exp(irf_store.(currpol)(ind_nn,1)));
        d_nu_pol = ((1-exp(irf_store.(currpol)(ind_na, 2))-exp(irf_store.(currpol)(ind_nn, 2))) - (1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1))));
        
        X = categorical({'Sector A','Sector N','Unemployed'});
        X = reordercats(X,{'Sector A','Sector N','Unemployed'});
        y = [d_na_nopol d_na_pol; d_nn_nopol d_nn_pol; d_nu_nopol d_nu_pol];
        subplot(3,2,jj)
        b = bar(X, y, 'FaceColor', 'flat'); title(pol_names{jj}), ylabel('Absolute Change')
        b(1).CData = [0 0 0];  
        b(2).CData = [1 1 1];       
    end
    legend('No Policy', 'Policy', 'Location', 'Northwest')
    suptitle('Mass of Workers')
    f.Position = [962 42 958 954];
    print -depsc figures/bar_number_workers
    
    % (4) Number of workers for the slides
    f=figure;
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        %     d_income_nopol = (irf_store.nopol(ind_income, 2)-irf_store.nopol(ind_income, 1))/detss.income;
        %     d_income_pol   = (irf_store.(currpol)(ind_income, 2)-irf_store.(currpol)(ind_income, 1))/detss.income;
        %
        %     income_a_nopol = exp(irf_store.nopol(ind_wa, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        %     income_n_nopol = exp(irf_store.nopol(ind_wn, :)).*(1-exp(irf_store.nopol(ind_taul, :)));
        %     income_u_nopol = exp(irf_store.nopol(ind_vsigma, :));
        %
        %     income_a_pol = exp(irf_store.(currpol)(ind_wa, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        %     income_n_pol = exp(irf_store.(currpol)(ind_wn, :)).*(1-exp(irf_store.(currpol)(ind_taul, :)));
        %     income_u_pol = exp(irf_store.(currpol)(ind_vsigma, :));
        %
        %     d_income_a_nopol = (income_a_nopol(2) - income_a_nopol(1))/income_a_nopol(1);
        %     d_income_a_pol   = (income_a_pol(2) - income_a_pol(1))/income_a_pol(1);
        %
        %     d_income_n_nopol = (income_n_nopol(2) - income_n_nopol(1))/income_n_nopol(1);
        %     d_income_n_pol   = (income_n_pol(2) - income_n_pol(1))/income_n_pol(1);
        %
        %     d_income_u_nopol = (income_u_nopol(2) - income_u_nopol(1))/income_u_nopol(1);
        %     d_income_u_pol   = (income_u_pol(2) - income_u_pol(1))/income_u_pol(1);
        
        d_na_nopol = (exp(irf_store.nopol(ind_na, 2)) - exp(irf_store.nopol(ind_na,1)));
        d_nn_nopol = (exp(irf_store.nopol(ind_nn, 2)) - exp(irf_store.nopol(ind_nn,1)));
        d_nu_nopol = ((1-exp(irf_store.nopol(ind_na, 2))-exp(irf_store.nopol(ind_nn, 2))) - (1-exp(irf_store.nopol(ind_na, 1))-exp(irf_store.nopol(ind_nn, 1))));
        
        d_na_pol = (exp(irf_store.(currpol)(ind_na, 2)) - exp(irf_store.(currpol)(ind_na,1)));
        d_nn_pol = (exp(irf_store.(currpol)(ind_nn, 2)) - exp(irf_store.(currpol)(ind_nn,1)));
        d_nu_pol = ((1-exp(irf_store.(currpol)(ind_na, 2))-exp(irf_store.(currpol)(ind_nn, 2))) - (1-exp(irf_store.(currpol)(ind_na, 1))-exp(irf_store.(currpol)(ind_nn, 1))));
        
        X = categorical({'Sector A','Sector N','Unemployed'});
        X = reordercats(X,{'Sector A','Sector N','Unemployed'});
        y = [d_na_nopol d_na_pol; d_nn_nopol d_nn_pol; d_nu_nopol d_nu_pol];
        subplot(2,3,jj)
        bar(X, y), title(pol_names{jj}), ylabel('Absolute Change')
    end
    legend('No Policy', 'Policy', 'Location', 'Northwest')
%     suptitle('Mass of Workers')
    f.Position = [962 42 958 954];
    print -depsc figures/bar_number_workers_slides
end

%% (11) US Fiscal Package

if options.us_package == 1
    
    % Run all policies
    oo_.exo_simul      = data_shocks;
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.all = oo_.endo_simul(:,1:end);
   
    % Figures
    figure
    for ii = 1:length(vars)
        currvar = vars{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names{ii};
        
        subplot(3,4,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.all(currind,1:Tirf)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.all(currind,1:Tirf))';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('level')
        else
            crisis    = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.all(currind,1:Tirf))/detss.(currvar)-1)';
            plot(time, [crisis, crisis_pol], 'Linewidth', 2), hold on
            plot(time, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            ylabel('% from SS')
        end
    end
    
    % Figure for paper
    vars_paper      = {'Bg'; 'unemp'; 'Cb'; 'N_a'; 'N_n'; 'R'; 'Fb'; 'sprqb'; 'f'; 'm'};
    var_names_paper = {'Public Debt'; 'Unemployment'; 'Cons. Borrower'; 'Labor Sector A'; 'Labor Sector N'; 'Policy Rate'; 'Default Rate'; 'Credit Spread'; 'Mass of Firms'; 'Entrants'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(5,2,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m') || strcmp(currvar,'ci')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.all(currind,1:Tirf)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.all(currind,1:Tirf))';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.all(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.all(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, crisis, 'Linewidth', 2), hold on
            plot(dateplot, crisis_pol, '--', 'Linewidth', 2), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    legend('No Policy', 'Ta', 'Location', 'Southeast')
    f.Position = [962 42 958 954];
    print -depsc figures/crisis_uspackage
    
    % Figure for slides
    vars_paper      = {'Bg'; 'G'; 'tau_l'; 'varsigma'; 'transfer'; 'govt_wage'};
    var_names_paper = {'Public Debt'; 'Govt. Consumption'; 'Income Tax'; 'Unemployment Insurance'; 'Transfer'; 'Liq. Assistance'};
    
    dt_start = datetime(2020,3,31);
    dt_end   = dt_start+calmonths(3*Tirf-1);
    dateplot    = dt_start:calmonths(3):dt_end;
    
    f=figure;
    for ii = 1:length(vars_paper)
        currvar = vars_paper{ii};
        currind = strmatch(currvar, M_.endo_names, 'exact');
        currname = var_names_paper{ii};
        subplot(2,3,ii)
        if strcmp(currvar,'T') || strcmp(currvar,'rot') || strcmp(currvar,'m') || strcmp(currvar,'ci') || strcmp(currvar,'transfer') || strcmp(currvar,'govt_wage')
            crisis = irf_store.nopol(currind,1:Tirf)';
            crisis_pol = irf_store.all(currind,1:Tirf)';
            plot(dateplot, [crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'R') || strcmp(currvar,'sprqb') || strcmp(currvar,'Fb')  || strcmp(currvar,'unemp')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.all(currind,1:Tirf))';
            plot(dateplot, [crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*ones(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('level')
        elseif strcmp(currvar,'Bg')
            crisis = exp(irf_store.nopol(currind,1:Tirf))';
            crisis_pol = exp(irf_store.all(currind,1:Tirf))';
            crisis_pol = 4*100*(crisis_pol-crisis)./detss.(currvar);
            plot(dateplot,crisis_pol, 'Linewidth', 2, 'color', [0.8500 0.3250 0.0980]), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname)
            axis tight
            grid minor
            ylabel('% difference, ann.')
        else
            crisis = 100*(exp(irf_store.nopol(currind,1:Tirf))/detss.(currvar)-1)';
            crisis_pol = 100*(exp(irf_store.all(currind,1:Tirf))/detss.(currvar)-1)';
            plot(dateplot, [crisis_pol], 'Linewidth', 3), hold on
            plot(dateplot, detss.(currvar)*zeros(size(time)), 'k--')
            title(currname, 'FontSize', 18)
            axis tight
            grid minor
            ylabel('% from SS')
        end
        set(gca, 'XTick', dateplot);
        datetick('x','QQ-YY','keeplimits','keepticks');
    end
    f.Position = [1 41 1920 963];
    print -depsc figures/crisis_uspackage_policies
    
    % G only
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(1:6,ind_G)  = data_shocks(1:6,ind_G);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.G = oo_.endo_simul(:,1:end);
    
    % UI only
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(1:6,ind_varsigma) = data_shocks(1:6,ind_varsigma);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.varsigma = oo_.endo_simul(:,1:end);
    
    % Transfer only
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(1:6,ind_transfer) = data_shocks(1:6,ind_transfer);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.transfer = oo_.endo_simul(:,1:end);
    
    % Govt wage only
    oo_.exo_simul(1:6,ind_ua) = data_shocks(1:6,ind_ua);
    oo_.exo_simul(1:6,ind_govt_wage) = data_shocks(1:6,ind_govt_wage);
    perfect_foresight_solver;
    oo_.exo_simul(:,:) = 0;

    irf_store.govt_wage = oo_.endo_simul(:,1:end);
    
    
    mult = struct;
    mult.gdp = struct;
    mult.cb  = struct;
    mult.cs  = struct;
    mult.y  = struct;
    mult.income = struct;
    pols = {'all'; 'G'; 'varsigma'; 'transfer'; 'govt_wage'};
    pol_names = {'All'; 'Govt. Cons.'; 'UI'; 'Transfer'; 'Firm Assistance'};
    
    ind_gdp   = strmatch('GDP', M_.endo_names, 'exact');
    ind_y   = strmatch('Y', M_.endo_names, 'exact');
    ind_cb    = strmatch('Cb', M_.endo_names, 'exact');
    ind_cs    = strmatch('Cs', M_.endo_names, 'exact');
    ind_R     = strmatch('R', M_.endo_names, 'exact');
    ind_spend_G = strmatch('spend_G', M_.endo_names, 'exact');
    ind_spend_tau_l = strmatch('spend_tau_l', M_.endo_names, 'exact');
    ind_spend_varsigma = strmatch('spend_varsigma', M_.endo_names, 'exact');
    ind_spend_transfer = strmatch('spend_transfer', M_.endo_names, 'exact');
    ind_spend_govt_wage = strmatch('spend_govt_wage', M_.endo_names, 'exact');
    ind_spend = [0; ind_spend_G; ind_spend_varsigma; ind_spend_transfer; ind_spend_govt_wage];
    ind_income    = strmatch('income', M_.endo_names, 'exact');
    ind_na = strmatch('N_a', M_.endo_names, 'exact');
    ind_nn = strmatch('N_n', M_.endo_names, 'exact');
    
    horz = 20;
    
    for jj = 1:length(pols)
        currpol = pols{jj};
        
        GDP_nopol = exp(irf_store.nopol(ind_gdp, 1:horz))';
        GDP_pol = exp(irf_store.(currpol)(ind_gdp, 1:horz))';
        
        cb_nopol = exp(irf_store.nopol(ind_cb, 1:horz))';
        cb_pol = exp(irf_store.(currpol)(ind_cb, 1:horz))';
        
        cs_nopol = exp(irf_store.nopol(ind_cs, 1:horz))';
        cs_pol = exp(irf_store.(currpol)(ind_cs, 1:horz))';
        
        y_nopol = exp(irf_store.nopol(ind_y, 1:horz))';
        y_pol = exp(irf_store.(currpol)(ind_y, 1:horz))';
        
        if strcmp(currpol,'all')
            spend_nopol = irf_store.nopol(ind_spend_G, 1:horz)'+irf_store.nopol(ind_spend_varsigma, 1:horz)'+irf_store.nopol(ind_spend_transfer, 1:horz)'+irf_store.nopol(ind_spend_govt_wage, 1:horz)';
            spend_pol   = irf_store.all(ind_spend_G, 1:horz)'+irf_store.all(ind_spend_varsigma, 1:horz)'+irf_store.all(ind_spend_transfer, 1:horz)'+irf_store.all(ind_spend_govt_wage, 1:horz)';
        else
            spend_nopol = irf_store.nopol(ind_spend(jj), 1:horz)';
            spend_pol   = irf_store.(currpol)(ind_spend(jj), 1:horz)';
        end
        
        income_nopol = irf_store.nopol(ind_income, 1:horz)';
        income_pol = irf_store.(currpol)(ind_income, 1:horz)';
        
        employ_nopol = exp(irf_store.nopol(ind_na, 1:horz))' + exp(irf_store.nopol(ind_nn, 1:horz))';
        employ_pol = exp(irf_store.(currpol)(ind_na, 1:horz))' + exp(irf_store.(currpol)(ind_nn, 1:horz))';
        
        R = [1; exp(irf_store.nopol(ind_R, 1:horz-1))'];
        cumR = 1./exp(cumsum(log(R)));
        
        mult.gdp.(currpol) = sum(cumR.*(GDP_pol-GDP_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.cb.(currpol) = sum(cumR.*(cb_pol-cb_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.cs.(currpol) = sum(cumR.*(cs_pol-cs_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.y.(currpol) = sum(cumR.*(y_pol-y_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.income.(currpol) = sum(cumR.*(income_pol-income_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult.employ.(currpol) = sum(cumR.*(employ_pol-employ_nopol))./sum(cumR.*(spend_pol-spend_nopol));
    end
    
    % Table with Income
%     fileID = fopen('tables/multipliers_us.txt','w');
%     fprintf(fileID, 'US Package \n');
%     fprintf(fileID,' All Policies &                 & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.all, mult.income.all, mult.cb.all, mult.cs.all, mult.gdp.all, '\\');
%     fprintf(fileID,' $G$         & Govt. Consumption & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.G, mult.income.G, mult.cb.G, mult.cs.G, mult.gdp.G, '\\');
% %     fprintf(' $\\tau_t^l$  & Income Tax        & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f%s \n', mult.employ.tau_l,mult.income.tau_l, mult.cb.tau_l, mult.cs.tau_l, mult.gdp.tau_l, '\\');
%     fprintf(fileID,' $\\varsigma$ & UI                & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.varsigma, mult.income.varsigma, mult.cb.varsigma, mult.cs.varsigma, mult.gdp.varsigma, '\\');
%     fprintf(fileID,' $T_t^b$     & Uncond. Transfer  & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.transfer, mult.income.transfer, mult.cb.transfer, mult.cs.transfer, mult.gdp.transfer, '\\');
%     fprintf(fileID,' $T_t^a$     & Liquidity Assist. & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.govt_wage, mult.income.govt_wage, mult.cb.govt_wage, mult.cs.govt_wage, mult.gdp.govt_wage, '\\');    
%     fclose(fileID);
end

%% (12) Plots with Multiplier "IRFs" for the US fiscal package

mult_ts = struct;
horz = linspace(1, 40, 40)';
for jj = 1:length(pols)
    for hh = 1:length(horz)
        currpol = pols{jj};

        GDP_nopol = exp(irf_store.nopol(ind_gdp, 1:hh))';
        GDP_pol = exp(irf_store.(currpol)(ind_gdp, 1:hh))';

        cb_nopol = exp(irf_store.nopol(ind_cb, 1:hh))';
        cb_pol = exp(irf_store.(currpol)(ind_cb, 1:hh))';

        cs_nopol = exp(irf_store.nopol(ind_cs, 1:hh))';
        cs_pol = exp(irf_store.(currpol)(ind_cs, 1:hh))';

        y_nopol = exp(irf_store.nopol(ind_y, 1:hh))';
        y_pol = exp(irf_store.(currpol)(ind_y, 1:hh))';
        
        if strcmp(currpol,'all')
            spend_nopol = irf_store.nopol(ind_spend_G, 1:hh)'+irf_store.nopol(ind_spend_varsigma, 1:hh)'+irf_store.nopol(ind_spend_transfer, 1:hh)'+irf_store.nopol(ind_spend_govt_wage, 1:hh)';
            spend_pol   = irf_store.all(ind_spend_G, 1:hh)'+irf_store.all(ind_spend_varsigma, 1:hh)'+irf_store.all(ind_spend_transfer, 1:hh)'+irf_store.all(ind_spend_govt_wage, 1:hh)';
        else
            spend_nopol = irf_store.nopol(ind_spend(jj), 1:hh)';
            spend_pol   = irf_store.(currpol)(ind_spend(jj), 1:hh)';
        end

        income_nopol = irf_store.nopol(ind_income, 1:hh)';
        income_pol = irf_store.(currpol)(ind_income, 1:hh)';

        employ_nopol = exp(irf_store.nopol(ind_na, 1:hh))' + exp(irf_store.nopol(ind_nn, 1:hh))';
        employ_pol = exp(irf_store.(currpol)(ind_na, 1:hh))' + exp(irf_store.(currpol)(ind_nn, 1:hh))';

        R = [1; exp(irf_store.nopol(ind_R, 1:hh-1))'];
        cumR = 1./exp(cumsum(log(R)));

        mult_ts.gdp.(currpol)(hh) = sum(cumR.*(GDP_pol-GDP_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult_ts.cb.(currpol)(hh) = sum(cumR.*(cb_pol-cb_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult_ts.cs.(currpol)(hh) = sum(cumR.*(cs_pol-cs_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult_ts.y.(currpol)(hh) = sum(cumR.*(y_pol-y_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult_ts.income.(currpol)(hh) = sum(cumR.*(income_pol-income_nopol))./sum(cumR.*(spend_pol-spend_nopol));
        mult_ts.employ.(currpol)(hh) = sum(cumR.*(employ_pol-employ_nopol))./sum(cumR.*(spend_pol-spend_nopol));
    end
end

for jj = 1:length(pols)
    currpol  = pols{jj};
    currname = pol_names{jj};
    h=figure
    subplot(2,2,1)
    plot(mult_ts.employ.(currpol), 'Linewidth', 2)
    title('Multiplier: Employment')
    subplot(2,2,2)
    plot(mult_ts.gdp.(currpol), 'Linewidth', 2)
    title('Multiplier: GDP')
    subplot(2,2,3)
    plot(mult_ts.cs.(currpol), 'Linewidth', 2)
    title('Multiplier: Cs')
    subplot(2,2,4)
    plot(mult_ts.cb.(currpol), 'Linewidth', 2)
    title('Multiplier: Cb')
    suptitle(sprintf('%s, CARES Act',currname))
    saveas(h,sprintf('figures/mult_cares_%s.png',currpol))
end

delete *.log *_dynamic.m *_results.mat *_variables.m *_static.m
rmdir('model', 's')
rmdir('+model', 's')