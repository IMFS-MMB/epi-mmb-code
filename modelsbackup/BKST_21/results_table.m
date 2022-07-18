 
% This code uses string variable filename_results_table
%
% If flag_append is one, this code will append the results to
% results_table_matrix. Otherwise, it will create variable
% results_table_matrix and set flag_append = 1

%-------------------------------------------------------------------------%
%                           Computes statistics                           %
%-------------------------------------------------------------------------%

% starts vector that will store information
X = [];

% weeks to peak
X = [X; t_peak_I];

% seriously ill peak
[~, t_peak_s] = max(M_s');
X = [X; t_peak_s'];

% seriously ill peak (by age)
for i_age = 1:n_age
    x = M_s(i_age, t_peak_s(i_age));
    X = [X; x];
end

% fraction dead (one year) deceased
for i_age = 1:n_age
    x = M_d(i_age, 52);
    X = [X; x];
end

x = frac_young * M_d(i_young, 52) + frac_old * M_d(i_old, 52);
X = [X; x];

% fraction dead (LR)
for i_age = 1:n_age
    x = M_d(i_age, T);
    X = [X; x];
end

x = frac_young * M_d(i_young, T) + frac_old * M_d(i_old, T);
X = [X; x];

% Fraction dead covid-19 LR - old/all
x = M_d(i_old, T) * frac_old / ...
    (M_d(i_old, T) * frac_old + M_d(i_young, T) * frac_young);
X = [X; x];

% lives saved rel to BM
if exist('S_table', 'var')==1
    for i = 1:n_age
        M = S_table.M_d(i_age, T);
        x = (1 - M_d(i_age, T))/(1 - M);
        X = [X; x];
    end
    
    % all agents
    M_y = S_table.M_d(i_young, T);
    M_o = S_table.M_d(i_old, T);
    M0 = frac_young * M_y + frac_old * M_o;
    M1 = frac_young * M_d(i_young, T) + frac_old * M_d(i_old, T);
    
    x = (1-M1)/(1-M0);
    X = [X; x];
else
    for i = 1:n_age
        x = 1;
        X = [X; x];
    end
    X = [X; 1];
end

% immune in LR
for i_age = 1:n_age
    x = M_r(i_age, T);
    X = [X; x];
end

x = frac_young * M_r(i_young, T) + frac_old * M_r(i_old, T);
X = [X; x];

% number of tests (peak, LR, max)
M = frac_young * M_t(i_young, :) + frac_old * M_t(i_old, :);
x = [M(t_peak_I); M(T); max(M)];
X = [X; x];

% total number of tests (1 year, 1.5 years, 2 years)
n_weeks_vec = [52, 78, 104];

for i_n_weeks = 1:length(n_weeks_vec)
    n_weeks = n_weeks_vec(i_n_weeks);
    if T >= n_weeks
        x = sum(M(1:n_weeks));
    else
        x = NaN;
    end
    X = [X; x];
end

% GDP gain / number of tests (1 year, 1.5 years)
if M(T)>0 && exist('S_table', 'var')==1
    US_gdp = 20.54*10^12;
    US_pop = 328.2*10^6;
    
    factor_vec = [1, 1.5]; % 1 year, 1.5 years
    
    for i_factor = 1:length(factor_vec)
        factor = factor_vec(i_factor);
        n_weeks = factor * 52;
        
        gdp_ND = income_ND * n_weeks;
        gdp_benchmark = sum(S_table.gdp(1:n_weeks));
        gdp_current_sim_USD = sum(gdp(1:n_weeks));
        
        gdp_ND_USD = factor * US_gdp;
        gdp_benchmark_USD = (gdp_ND_USD / gdp_ND) * gdp_benchmark;
        gdp_current_sim_USD = (gdp_ND_USD / gdp_ND) * gdp_current_sim_USD;
        
        nr_tests = sum(M(1:n_weeks)) * US_pop;
        
        if gdp_current_sim_USD - gdp_benchmark_USD > 0
            x = (gdp_current_sim_USD - gdp_benchmark_USD) / nr_tests;
        else
            x = NaN;
        end
        
        X = [X; x];
    end
else
    x = NaN(2, 1);
    X = [X; x];
end

% GDP at peak
x = gdp(t_peak_I);
X = [X; x];

% GDP at peak - rel to ND
x = gdp(t_peak_I)/income_ND;
X = [X; x];

% GDP at peak - rel to BM
if exist('S_table', 'var')==1
    x = gdp(t_peak_I)/S_table.gdp(S_table.t_peak_I);
    X = [X; x];
else
    x = 1;
    X = [X; x];
end

% GDP in one year (cumulative)
x = sum(gdp(1:52));
X = [X; x];

x = sum(gdp(1:52))/(52 * income_ND);
X = [X; x];

if exist('S_table', 'var')==1
    x = sum(gdp(1:52))/sum(S_table.gdp(1:52));
    X = [X; x];
else
    x = 1;
    X = [X; x];
end

% GDP per capita at peak
x = gdp_pc(t_peak_I);
X = [X; x];

% GDP pc at peak - rel. to ND
x = gdp_pc(t_peak_I)/income_ND;
X = [X; x];

% GDP pc at peak - rel. to BM
if exist('S_table', 'var')==1
    x = gdp_pc(t_peak_I)/S_table.gdp_pc(S_table.t_peak_I);
    X = [X; x];
else
    x = 1;
    X = [X; x];
end

% GDP per capita in one year (cumulative)
x = sum(gdp_pc(1:52));
X = [X; x];

x = sum(gdp_pc(1:52))/(52 * income_ND);
X = [X; x];

if exist('S_table', 'var')==1
    x = sum(gdp_pc(1:52))/sum(S_table.gdp_pc(1:52));
    X = [X; x];
else
    x = 1;
    X = [X; x];
end

% cost of saved lives
if exist('S_table', 'var')==1
    t0 = 78;
    
    US_gdp = 20.54*10^12;
    US_pop = 328.2*10^6;

    gdp_1y_BM = sum(S_table.gdp(1:t0));
    gdp_1y = sum(gdp(1:t0));
    
    frac_dead_1y_BM = frac_young * S_table.M_d(i_young, t0) + frac_old * S_table.M_d(i_old, t0);
    frac_dead_1y = frac_young * M_d(i_young, t0) + frac_old * M_d(i_old, t0);
    
    var_gdp = (gdp_1y - gdp_1y_BM)/gdp_1y_BM * US_gdp;
    lives_saved = -(frac_dead_1y - frac_dead_1y_BM) * US_pop;
    
    if lives_saved > 0
        x = -(var_gdp/lives_saved)/10^6;
    else
        x = NaN;
    end
    
    X = [X; x];
else
    x = NaN;
    X = [X; x];
end

% time at home rel. to ND
for i = 1:3
    if i==1
        t = t_peak_I;
    elseif i==2
        t = 26; % 6 months = 26 weeks
    elseif i==3
        t = 52;
    end
    
    for i_age = 1:n_age
        x = (d_h(i_age, t) + v_h(i_age, t))/(d_h(i_age, T) + v_h(i_age, T));
        X = [X; x];
    end
    
    x0 = frac_young * (d_h(i_young, t) + v_h(i_young, t)) + frac_old * (d_h(i_old, t) + v_h(i_old, t));
    x1 = frac_young * (d_h(i_young, T) + v_h(i_young, T)) + frac_old * (d_h(i_old, T) + v_h(i_old, T));
    x = x0/x1;
    X = [X; x];
end

% weekly time at home
for i = 1:3
    if i==1
        t = t_peak_I;
    elseif i==2
        t = 26; % 6 months = 26 weeks
    elseif i==3
        t = 52;
    end
    
    for i_age = 1:n_age
        x = (d_h(i_age, t) + v_h(i_age, t)) * 112;
        X = [X; x];
    end
    
    x = frac_young * (d_h(i_young, t) + v_h(i_young, t)) + frac_old * (d_h(i_old, t) + v_h(i_old, t));
    x = x * 112;
    X = [X; x];
end

% value of healthy (without lambda_p)
for i_age = 1:n_age
    x = V_h_private(i_age, 1);
    X = [X; x];
end

x = frac_young * V_h_private(i_young, 1) + frac_old * V_h_private(i_old, 1);
X = [X; x];

% value of healthy (includes lambda_p)
for i_age = 1:n_age
    x = V_h(i_age, 1);
    X = [X; x];
end

x = frac_young * V_h(i_young, 1) + frac_old * V_h(i_old, 1);
X = [X; x];

% R0
R0a = zeros(n_age, 1);

for i_age = 1:n_age
    T_i = (phi(i_nosymptoms, i_age) + ...
        (1-phi(i_nosymptoms, i_age))*alpha(i_age)) * Delta(i_age);
    T_i = T_i + (1-Delta(i_age));
    T_i = T_i^-1;
    
    T_s = (phi(i_symptoms, i_age) + ...
        (1-phi(i_symptoms, i_age)) * delta_1(i_age)) * Delta(i_age);
    T_s = T_s + (1-Delta(i_age));
    T_s = T_s^-1;
    
    P_s = (1-phi(i_nosymptoms, i_age))*alpha(i_age)*Delta(i_age);
    P_s = P_s/(1-(1-phi(i_nosymptoms, i_age))*(1-alpha(i_age))*Delta(i_age));
    
    n_a = n_h(i_age, T) + l_h(i_age, T);  
    n_bar = frac_young * (n_h(i_young, T) + l_h(i_young, T)) + ...
            frac_old   * (n_h(i_old  , T) + l_h(i_old  , T));
    
    R0a(i_age) = (n_a * T_i + n_plus_l_s * P_s * T_s) * n_bar * Pi0;
end

R0 = frac_young * R0a(i_young) + frac_old * R0a(i_old);

x = R0;
X = [X; x];

% R
for t = [t_peak_I, max(t_peak_I-4, 1), min(t_peak_I+4, T)]
    
    susceptible = frac_young * (M_h(i_young, t) + M_fh(i_young, t)) + ...
                  frac_old   * (M_h(i_old  , t) + M_fh(i_old  , t));
    
    alive = frac_young * (1 - M_d(i_young, t) - M_dn(i_young, t)) + ...
            frac_old   * (1 - M_d(i_old  , t) - M_dn(i_old  , t));
    
    x = R0 * susceptible/alive;
    X = [X; x];
end

% CEV
if exist('S_c_eq', 'var')==1
    for i_age = 1:n_age
        factor_c_eq = c_eq(V_h_private, i_age);
        x = (factor_c_eq - 1) * 100;
        X = [X; x];
    end
else
    for i_age = 1:n_age
        X = [X; NaN];
    end
end

%-------------------------------------------------------------------------%
%                             Exports results                             %
%-------------------------------------------------------------------------%
%{
% appends results
if exist('flag_append', 'var')==1 && flag_append==1
    results_table_matrix = [results_table_matrix, X];
else
    results_table_matrix = X;
end

flag_append = 1;

% saves results
filename = ['tables/', filename_results_table, '.csv'];
% csvwrite(filename, results_table_matrix)
dlmwrite(filename, results_table_matrix, 'delimiter', ',', 'precision', 9);
%}










