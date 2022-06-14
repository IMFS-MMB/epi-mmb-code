%%% This file computes fiscal multipliers without the Covid shock
% 
% clear
% clc
% close all

load params.mat

% Run Dynare once
dynare model.mod noclearall;

% Indices for variables of interest
ind_unemp  = strmatch('unemp', M_.endo_names, 'exact');
ind_ua_var = strmatch('util_a', M_.endo_names, 'exact');
ind_ua     = strmatch('ee_util_a', M_.exo_names, 'exact');

ind_G         = strmatch('ee_G', M_.exo_names, 'exact');
ind_varsigma  = strmatch('ee_varsigma', M_.exo_names, 'exact');
ind_govt_wage = strmatch('ee_govt_wage', M_.exo_names, 'exact');
ind_transfer  = strmatch('ee_transfer', M_.exo_names, 'exact');
ind_tau_l     = strmatch('ee_taul', M_.exo_names, 'exact');

%% (1) Baseline, no policy
oo_.exo_simul(:,ind_ua) = 0;
perfect_foresight_solver;

irf_store.nopol = oo_.endo_simul(:,1:end);


% %% (2) Govt Spending
% oo_.exo_simul(2,ind_G) = 0.01;
% perfect_foresight_solver;
% 
% irf_store.G = oo_.endo_simul(:,1:end);
% 
% oo_.exo_simul(2,ind_G) = 0;
% 
% %% (3) Tax cut
% oo_.exo_simul(2,ind_tau_l) = 0.01;
% perfect_foresight_solver;
% 
% irf_store.tau_l = oo_.endo_simul(:,1:end);
% 
% oo_.exo_simul(2,ind_tau_l) = 0;
% 
% %% (4) Unemployment insurance
% oo_.exo_simul(2,ind_varsigma) = 0.01;
% perfect_foresight_solver;
% 
% irf_store.varsigma = oo_.endo_simul(:,1:end);
% 
% oo_.exo_simul(2,ind_varsigma) = 0;
% 
% %% (5) Transfer
% oo_.exo_simul(2,ind_transfer) = 0.01;
% perfect_foresight_solver;
% 
% irf_store.transfer = oo_.endo_simul(:,1:end);
% 
% oo_.exo_simul(2,ind_transfer) = 0;
% 
% %% (6) Payroll subsidies
% oo_.exo_simul(2,ind_govt_wage) = 0.01;
% perfect_foresight_solver;
% 
% irf_store.govt_wage = oo_.endo_simul(:,1:end);
% 
% oo_.exo_simul(2,ind_govt_wage) = 0;

%% (7) Table with multipliers

% pols = {'G'; 'tau_l'; 'varsigma'; 'transfer'; 'govt_wage'};
% pol_names = {'Govt. Cons.'; 'Income Tax'; 'UI'; 'Transfer'; 'Firm Assistance'};
% 
% mult = struct;
% mult.gdp = struct;
% mult.cb  = struct;
% mult.cs  = struct;
% mult.y  = struct;
% mult.income = struct;
% 
% ind_gdp   = strmatch('GDP', M_.endo_names, 'exact');
% ind_y     = strmatch('Y', M_.endo_names, 'exact');
% ind_cb    = strmatch('Cb', M_.endo_names, 'exact');
% ind_cs    = strmatch('Cs', M_.endo_names, 'exact');
% ind_R     = strmatch('R', M_.endo_names, 'exact');
% ind_spend_G = strmatch('spend_G', M_.endo_names, 'exact');
% ind_spend_tau_l = strmatch('spend_tau_l', M_.endo_names, 'exact');
% ind_spend_varsigma = strmatch('spend_varsigma', M_.endo_names, 'exact');
% ind_spend_transfer = strmatch('spend_transfer', M_.endo_names, 'exact');
% ind_spend_govt_wage = strmatch('spend_govt_wage', M_.endo_names, 'exact');
% ind_spend = [ind_spend_G; ind_spend_tau_l; ind_spend_varsigma; ind_spend_transfer; ind_spend_govt_wage];
% ind_income    = strmatch('income', M_.endo_names, 'exact');
% ind_na = strmatch('N_a', M_.endo_names, 'exact');
% ind_nn = strmatch('N_n', M_.endo_names, 'exact');

horz = 20;

for jj = 1:length(pols)
    currpol = pols{jj};

    GDP_nopol = exp(irf_store.nopol(ind_gdp, 1:horz))';
%     GDP_pol = exp(irf_store.(currpol)(ind_gdp, 1:horz))';

    cb_nopol = exp(irf_store.nopol(ind_cb, 1:horz))';
%     cb_pol = exp(irf_store.(currpol)(ind_cb, 1:horz))';

    cs_nopol = exp(irf_store.nopol(ind_cs, 1:horz))';
%     cs_pol = exp(irf_store.(currpol)(ind_cs, 1:horz))';

    y_nopol = exp(irf_store.nopol(ind_y, 1:horz))';
%     y_pol = exp(irf_store.(currpol)(ind_y, 1:horz))';

    spend_nopol = irf_store.nopol(ind_spend(jj), 1:horz)';
%     spend_pol = irf_store.(currpol)(ind_spend(jj), 1:horz)';

    income_nopol = irf_store.nopol(ind_income, 1:horz)';
%     income_pol = irf_store.(currpol)(ind_income, 1:horz)';

    employ_nopol = exp(irf_store.nopol(ind_na, 1:horz))' + exp(irf_store.nopol(ind_nn, 1:horz))';
%     employ_pol = exp(irf_store.(currpol)(ind_na, 1:horz))' + exp(irf_store.(currpol)(ind_nn, 1:horz))';

    R = [1; exp(irf_store.nopol(ind_R, 1:horz-1))'];
    cumR = 1./exp(cumsum(log(R)));
% 
%     mult.gdp.(currpol) = sum(cumR.*(GDP_pol-GDP_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     mult.cb.(currpol) = sum(cumR.*(cb_pol-cb_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     mult.cs.(currpol) = sum(cumR.*(cs_pol-cs_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     mult.y.(currpol) = sum(cumR.*(y_pol-y_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     mult.income.(currpol) = sum(cumR.*(income_pol-income_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     mult.employ.(currpol) = sum(cumR.*(employ_pol-employ_nopol))./sum(cumR.*(spend_pol-spend_nopol));
end

% Table with multipliers
%fileID = fopen('tables/multipliers_normal.txt','w');
% fprintf(fileID, ' $G$         & Govt. Consumption & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.G, mult.income.G, mult.cb.G, mult.cs.G, mult.gdp.G, '\\');
% fprintf(fileID, ' $\\tau_t^l$  & Income Tax        & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.tau_l, mult.income.tau_l, mult.cb.tau_l, mult.cs.tau_l, mult.gdp.tau_l, '\\');
% fprintf(fileID, ' $\\varsigma$ & UI                & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n',mult.employ.varsigma, mult.income.varsigma, mult.cb.varsigma, mult.cs.varsigma, mult.gdp.varsigma, '\\');
% fprintf(fileID, ' $T_t^b$     & Uncond. Transfer  & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.transfer, mult.income.transfer, mult.cb.transfer, mult.cs.transfer, mult.gdp.transfer, '\\');
% fprintf(fileID, ' $T_t^a$     & Liquidity Assist. & %4.4f & %4.4f & %4.4f & %4.4f & %4.4f %s \n', mult.employ.govt_wage, mult.income.govt_wage, mult.cb.govt_wage, mult.cs.govt_wage, mult.gdp.govt_wage, '\\');
% fclose(fileID);

%% (8) Plots with Multiplier "IRFs"

% mult_ts = struct;
% horz = linspace(1, 40, 40)';
% for jj = 1:length(pols)
%     for hh = 1:length(horz)
%         currpol = pols{jj};
% 
%         GDP_nopol = exp(irf_store.nopol(ind_gdp, 1:hh))';
%         GDP_pol = exp(irf_store.(currpol)(ind_gdp, 1:hh))';
% 
%         cb_nopol = exp(irf_store.nopol(ind_cb, 1:hh))';
%         cb_pol = exp(irf_store.(currpol)(ind_cb, 1:hh))';
% 
%         cs_nopol = exp(irf_store.nopol(ind_cs, 1:hh))';
%         cs_pol = exp(irf_store.(currpol)(ind_cs, 1:hh))';
% 
%         y_nopol = exp(irf_store.nopol(ind_y, 1:hh))';
%         y_pol = exp(irf_store.(currpol)(ind_y, 1:hh))';
% 
%         spend_nopol = irf_store.nopol(ind_spend(jj), 1:hh)';
%         spend_pol = irf_store.(currpol)(ind_spend(jj), 1:hh)';
% 
%         income_nopol = irf_store.nopol(ind_income, 1:hh)';
%         income_pol = irf_store.(currpol)(ind_income, 1:hh)';
% 
%         employ_nopol = exp(irf_store.nopol(ind_na, 1:hh))' + exp(irf_store.nopol(ind_nn, 1:hh))';
%         employ_pol = exp(irf_store.(currpol)(ind_na, 1:hh))' + exp(irf_store.(currpol)(ind_nn, 1:hh))';
% 
%         R = [1; exp(irf_store.nopol(ind_R, 1:hh-1))'];
%         cumR = 1./exp(cumsum(log(R)));
% 
%         mult_ts.gdp.(currpol)(hh) = sum(cumR.*(GDP_pol-GDP_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%         mult_ts.cb.(currpol)(hh) = sum(cumR.*(cb_pol-cb_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%         mult_ts.cs.(currpol)(hh) = sum(cumR.*(cs_pol-cs_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%         mult_ts.y.(currpol)(hh) = sum(cumR.*(y_pol-y_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%         mult_ts.income.(currpol)(hh) = sum(cumR.*(income_pol-income_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%         mult_ts.employ.(currpol)(hh) = sum(cumR.*(employ_pol-employ_nopol))./sum(cumR.*(spend_pol-spend_nopol));
%     end
% end
%
% 
% 
% for jj = 1:length(pols)
%     currpol  = pols{jj};
%     currname = pol_names{jj};
%     h=figure
%     subplot(2,2,1)
%     plot(mult_ts.employ.(currpol), 'Linewidth', 2)
%     title('Multiplier: Employment')
%     subplot(2,2,2)
%     plot(mult_ts.gdp.(currpol), 'Linewidth', 2)
%     title('Multiplier: GDP')
%     subplot(2,2,3)
%     plot(mult_ts.cs.(currpol), 'Linewidth', 2)
%     title('Multiplier: Cs')
%     subplot(2,2,4)
%     plot(mult_ts.cb.(currpol), 'Linewidth', 2)
%     title('Multiplier: Cb')
%     suptitle(sprintf('%s, Normal Times',currname))
%     saveas(h,sprintf('figures/mult_normal_%s.png',currpol))
% end
% 
% delete *.log *_dynamic.m *_results.mat *_variables.m *_static.m
% rmdir('model', 's')
% rmdir('+model', 's')
