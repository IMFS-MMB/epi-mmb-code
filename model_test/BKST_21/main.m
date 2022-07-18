
% h: healthy (susceptible)
% s: serious symptoms (hospitalized)

%clear
%clc

%tic

%-------------------------------------------------------------------------%
%                                 Moments                                 %
%-------------------------------------------------------------------------%

total_work_time_data = 0.357;
v_data = 0.08 * total_work_time_data;
n_data = total_work_time_data - v_data;
l_data = 0.154; % leisure time in the no-disease case
x_over_income_data = 0.125; % x/income in the no-disease case
Z_data = 0.00032607; % hospital beds
zeta_data = 0.42; % fraction of activities that can be separated

%-------------------------------------------------------------------------%
%                           Technical parameters                          %
%-------------------------------------------------------------------------%

n_age = 2; % number of age groups
i_young = 1;
i_old = 2;
i_nosymptoms = 1;
i_symptoms = 2;
T = 100;
% T = 4 * 52;
tol = 1e-5;

%-------------------------------------------------------------------------%
%                    Internally calibrated parameters                     %
%-------------------------------------------------------------------------%

%----------------%
% No teleworking %
%----------------%

% theta = 0.03327394;
% gamma = 0.63589003;
% lambda_d = 1.5654178;
% lambda_i = 4.502474 - lambda_d;
% alpha1 = 0;
% alpha2 = 0;
% Pi_star = 0.106669183;
% Pi0 = 12.127;
% b = 5.6; % calibrated to fit Swedish data
% % b = 11.0121; % calibrated to get VSL of 9.3M

%------------------%
% With teleworking %
%------------------%

theta = 0.0334618228741673;
gamma = 0.634789816062878;
lambda_d = 1.56193260834488;
alpha1 = 1.05481946598541; % teleworking parameter
alpha2 = 0.959714181749383; % teleworking parameter: output is v * (alpha1 - alpha2 * v)
lambda_i = 1.0685; % 50% increase time at home if infected
Pi_star = 0.112617662;
Pi0 = 13.425;
b = 11.0049; % calibrated to get VSL of 9.3M
% b = 4.63; % calibrated to fit Swedish data

%-------------------------------------------------------------------------%
%                    Externally calibrated parameters                     %
%-------------------------------------------------------------------------%

% parameters - demography
frac_old = 0.214;
frac_young = 1 - frac_old;
frac_age = [frac_young, frac_old];

% parameters - utility
rho = -1.72; % 1/(1-rho) is the elasticity of subst bw goods and leisure time
time_discount = 0.96^(1/52); % discount factor
survival_young = 1;
survival_old = 0.95^(1/52);

% survival vector
Delta = zeros(n_age, 1);
Delta(i_young) = survival_young;
Delta(i_old) = survival_old;

% effective discount (beta)
beta = time_discount * Delta;

% parameters - infection
helper=load('inf_ini.mat');
initial_sick=helper.helper;
%initial_sick = 0.0001; % initial fraction of people of each age that is infected
n_plus_l_s = 0.158159722; % n+l por sick people (affects aggregate infection)

% parameters - wages
w = 1; % basic hourly wage
income_old_over_income_young = 0.6;

% parameters - disease
delta_1 = zeros(n_age, 1); % prob. of dying cond. on developing covid symptoms
delta_1(i_young) = 0.06553201;
delta_1(i_old) = 0.73799969;

delta_2 = zeros(n_age, 1); % prob. dying if symptoms and no hospital bed
delta_2(:) = 1;

alpha = zeros(n_age, 1); % prob. of developing symptoms
alpha(:) = 1;

phi = zeros(2, n_age); % prob. of recovering (2 is for whether you have symptoms or not)
phi(i_nosymptoms, i_young) = 0.96673931;
phi(i_nosymptoms, i_old) = 0.90896255;
phi(i_symptoms, i_young) = 0.28409091;
phi(i_symptoms, i_old) = 0.28409091;

Z = 1; % use this if I don't use hospital bed constraints
% Z = Z_data;

% t_vaccine = T + 1; % in case we don't want a vaccine
t_vaccine = 1.5 * 52;

% parameters - selective mixing
zeta = 0; % fraction of activities that are separated
vartheta = [frac_young, frac_old]; % fraction of space dedicated to given age group

%-------------------------------------------------------------------------%
%                                  Flags                                  %
%-------------------------------------------------------------------------%

% if 1, uses behavior in the case of no disease
flag_epidemiological = 0;

% if 1, value functions of young have the following variables of old 
% agents: delta_1, phi, beta, alpha
flag_fake_young = 0;

% if 1, doesn't iterate on Pi and uses Pi stored in struct S_Pi
flag_Pi = 0;

% if 1, uses Pi and M_s stored in memory as initial guesses for new
% equilibrium
flag_Pi_guess = 0;

%-------------------------------------------------------------------------%
%                             Runs something                              %
%-------------------------------------------------------------------------%

do = 2; %14

if do==1  % to calibrate b fitting Swedish data (15.7% increase in d+v in week 9)
    
%     b = 5.6;
    
    equilibrium
    
    % time at home (no disease)
    time_home_ND = 0;
    
    for i_age = 1:n_age
        time_home_ND = time_home_ND + ...
            frac_age(i_age) * (d_h(i_age, T) + v_h(i_age, T));
    end
    
    % time at home (path during disease)
    time_home = zeros(1, T);
    
    for i_age = 1:n_age
        time_home = time_home + frac_age(i_age) * ...
            (M_h(i_age, :) .* (d_h(i_age, :) + v_h(i_age, :)) + ...
             M_f(i_age, :) .* (d_f(i_age, :) + v_f(i_age, :)) + ...
             M_r(i_age, :) .* (d_r(i_age, :) + v_r(i_age, :)) + ...
             M_i(i_age, :) .* (d_i(i_age, :) + v_i(i_age, :))) ...
             ./ (M_h(i_age, :) + M_f(i_age, :) + M_r(i_age, :) + M_i(i_age, :));
    end
    
    % time at home (percentage change)
    time_home_pctchange = 100 * (time_home/time_home_ND - 1);
    
    % prints info
    fprintf('increase in d+v in week 9 (pct. change): %.8f\n', time_home_pctchange(9))
    
    % only young healthy
%     x = d_h + v_h;
%     X = 100 * (x(i_young, :)/x(i_young, T) - 1);
    
    % prints info
%     fprintf('increase in d+v in week 9 (young, pct. change): %.8f\n', X(9))

elseif do==2
    
    % selected policies (tables, figures) for 4 calibration strategies
    % (with/without teleworking, high/low b)
    
%     for i_simulation = 1:4
    for i_simulation = 2:2
        if i_simulation==1 % with telworking, low b
            filename_results_table = 'tw1_bL';
            theta = 0.0334618228741673;
            gamma = 0.634789816062878;
            lambda_d = 1.56193260834488;
            alpha1 = 1.05481946598541;
            alpha2 = 0.959714181749383;
            lambda_i = 1.0685; % 50% increase time at home if infected
            Pi_star = 0.112617662;
            Pi0 = 13.425;
            b = 4.63; % calibrated to fit Swedish data
        elseif i_simulation==2 % with telworking, high b
            filename_results_table = 'tw1_bH';
            theta = 0.0334618228741673;
            gamma = 0.634789816062878;
            lambda_d = 1.56193260834488;
            alpha1 = 1.05481946598541;
            alpha2 = 0.959714181749383;
            lambda_i = 1.0685; % 50% increase time at home if infected
            Pi_star = 0.112617662;
            Pi0 = 13.425;
            b = 11.0049; % calibrated to get VSL of 9.3M
        elseif i_simulation==3 % without telworking, low b
            filename_results_table = 'tw0_bL';
            theta = 0.03327394;
            gamma = 0.63589003;
            lambda_d = 1.5654178;
            lambda_i = 4.502474 - lambda_d;
            alpha1 = 0;
            alpha2 = 0;
            Pi_star = 0.106669183;
            Pi0 = 12.127;
            b = 5.6; % calibrated to fit Swedish data
        else % without telworking, high b
            filename_results_table = 'tw0_bH';
            theta = 0.03327394;
            gamma = 0.63589003;
            lambda_d = 1.5654178;
            lambda_i = 4.502474 - lambda_d;
            alpha1 = 0;
            alpha2 = 0;
            Pi_star = 0.106669183;
            Pi0 = 12.127;
            b = 11.0121; % calibrated to get VSL of 9.3M
        end
        
        policies
    end

elseif do==3

    % Some numbers for the "no disease" column in the benchmark results table

    % computes benchmark equilibrium (COVID world)
    equilibrium

    % initiates X
    X = [];

    % GDP at peak - rel. to BM
    x = income_ND / gdp(t_peak_I);
    X = [X; x];

    % GDP 1 year - rel. to BM
    x = 52 * income_ND / sum(gdp(1:52));
    X = [X; x];

    % hrs. home (yng) - peak
    x = 112 * (d_h_term(i_young) + v_h_term(i_young));
    X = [X; x];

    % hrs. home (old) - peak
    x = 112 * (d_h_term(i_old) + v_h_term(i_old));
    X = [X; x];

    % hrs. home (yng) - 6m
    x = 112 * (d_h_term(i_young) + v_h_term(i_young));
    X = [X; x];

    % hrs. home (old) - 6m
    x = 112 * (d_h_term(i_old) + v_h_term(i_old));
    X = [X; x];

    % value - healthy (yng)
    x = V_h_term(i_young);
    X = [X; x];

    % value - healthy (old)
    x = V_h_term(i_old);
    X = [X; x];

    % value - healthy (all)
    x = frac_young * V_h_term(i_young) + frac_old * V_h_term(i_old);
    X = [X; x];

    disp(X)

elseif do==4
    
    planner
    
elseif do==5 % simulates some optimal lockdowns (4 weeks)
    
    filename_results_table = 'optimal_4week';
    
    % parameters: tw1, bH
    theta = 0.0334618228741673;
    gamma = 0.634789816062878;
    lambda_d = 1.56193260834488;
    alpha1 = 1.05481946598541;
    alpha2 = 0.959714181749383;
    lambda_i = 1.0685; % 50% increase time at home if infected
    Pi_star = 0.112617662;
    Pi0 = 13.425;
    b = 11.0049; % calibrated to get VSL of 9.3M
    
    % 1: benchmark
    equilibrium
    results_table
    
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    % 2: universal
    lambda_p = [0.952714326, 0.952714326, 0.952714326, 0.952714326, 0.74890537, 0.74890537, 0.74890537, 0.74890537, 0.746958584, 0.746958584, 0.746958584, 0.746958584, 0.733525518, 0.733525518, 0.733525518, 0.733525518, 0.722812661, 0.722812661, 0.722812661, 0.722812661, 0.713482993, 0.713482993, 0.713482993, 0.713482993, 0.704964309, 0.704964309, 0.704964309, 0.704964309, 0.697182255, 0.697182255, 0.697182255, 0.697182255, 0.689659299, 0.689659299, 0.689659299, 0.689659299, 0.682115039, 0.682115039, 0.682115039, 0.682115039, 0.674397618, 0.674397618, 0.674397618, 0.674397618, 0.666077673, 0.666077673, 0.666077673, 0.666077673, 0.656537981, 0.656537981, 0.656537981, 0.656537981, 0.645325248, 0.645325248, 0.645325248, 0.645325248, 0.631278323, 0.631278323, 0.631278323, 0.631278323, 0.612068875, 0.612068875, 0.612068875, 0.612068875, 0.58360761, 0.58360761, 0.58360761, 0.58360761, 0.539718682, 0.539718682, 0.539718682, 0.539718682, 0.230422777, 0.230422777, 0.230422777, 0.230422777, 0.230422777, 0.230422777];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp0_univ';
    figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp0_univ';
    figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp0_univ';
    figures_report
    
    % 3: age specific
    lambda_p = [0.970508218, 0.970508218, 0.970508218, 0.970508218, 0.758040022, 0.758040022, 0.758040022, 0.758040022, 0.75806864, 0.75806864, 0.75806864, 0.75806864, 0.745269026, 0.745269026, 0.745269026, 0.745269026, 0.724599608, 0.724599608, 0.724599608, 0.724599608, 0.727676653, 0.727676653, 0.727676653, 0.727676653, 0.713283369, 0.713283369, 0.713283369, 0.713283369, 0.702679432, 0.702679432, 0.702679432, 0.702679432, 0.707093519, 0.707093519, 0.707093519, 0.707093519, 0.690472405, 0.690472405, 0.690472405, 0.690472405, 0.685313329, 0.685313329, 0.685313329, 0.685313329, 0.679425038, 0.679425038, 0.679425038, 0.679425038, 0.667398045, 0.667398045, 0.667398045, 0.667398045, 0.654670493, 0.654670493, 0.654670493, 0.654670493, 0.639193325, 0.639193325, 0.639193325, 0.639193325, 0.615437052, 0.615437052, 0.615437052, 0.615437052, 0.587684866, 0.587684866, 0.587684866, 0.587684866, 0.563533574, 0.563533574, 0.563533574, 0.563533574, 0.232571203, 0.232571203, 0.232571203, 0.232571203, 0.232571203, 0.232571203];
    lambda_p(2, :) = [0.664120302, 0.664120302, 0.664120302, 0.664120302, 0.471075544, 0.471075544, 0.471075544, 0.471075544, 0.463235167, 0.463235167, 0.463235167, 0.463235167, 0.533073818, 0.533073818, 0.533073818, 0.533073818, 0.451374354, 0.451374354, 0.451374354, 0.451374354, 0.605783386, 0.605783386, 0.605783386, 0.605783386, 0.417647161, 0.417647161, 0.417647161, 0.417647161, 0.421474431, 0.421474431, 0.421474431, 0.421474431, 0.471715951, 0.471715951, 0.471715951, 0.471715951, 0.447568956, 0.447568956, 0.447568956, 0.447568956, 0.507029106, 0.507029106, 0.507029106, 0.507029106, 0.480841769, 0.480841769, 0.480841769, 0.480841769, 0.436940865, 0.436940865, 0.436940865, 0.436940865, 0.392864232, 0.392864232, 0.392864232, 0.392864232, 0.397521824, 0.397521824, 0.397521824, 0.397521824, 0.368621428, 0.368621428, 0.368621428, 0.368621428, 0.344828721, 0.344828721, 0.344828721, 0.344828721, 0.185175603, 0.185175603, 0.185175603, 0.185175603, 0.171600211, 0.171600211, 0.171600211, 0.171600211, 0.171600211, 0.171600211];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp0_agespec';
    figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp0_agespec';
    figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp0_agespec';
    figures_report
    
    % 4: hospital, benchmark
    clear lambda_p_h lambda_p_i lambda_p_r lambda_p_f
    clear S_fig S_table
    Z = Z_data;
    
    equilibrium
    results_table
    
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    % 5: universal
    lambda_p = [1.02163575, 1.02163575, 1.02163575, 1.02163575, 0.693723342, 0.693723342, 0.693723342, 0.693723342, 0.76213465, 0.76213465, 0.76213465, 0.76213465, 0.692492853, 0.692492853, 0.692492853, 0.692492853, 0.731795457, 0.731795457, 0.731795457, 0.731795457, 0.736977949, 0.736977949, 0.736977949, 0.736977949, 0.736342921, 0.736342921, 0.736342921, 0.736342921, 0.650508524, 0.650508524, 0.650508524, 0.650508524, 0.750587187, 0.750587187, 0.750587187, 0.750587187, 0.589146655, 0.589146655, 0.589146655, 0.589146655, 0.803407409, 0.803407409, 0.803407409, 0.803407409, 0.597397014, 0.597397014, 0.597397014, 0.597397014, 0.632720655, 0.632720655, 0.632720655, 0.632720655, 0.757373036, 0.757373036, 0.757373036, 0.757373036, 0.605434979, 0.605434979, 0.605434979, 0.605434979, 0.629341488, 0.629341488, 0.629341488, 0.629341488, 0.568801335, 0.568801335, 0.568801335, 0.568801335, 0.535595526, 0.535595526, 0.535595526, 0.535595526, 0.138088665, 0.138088665, 0.138088665, 0.138088665, 0.138088665, 0.138088665];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp1_univ';
    figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp1_univ';
    figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp1_univ';
    figures_report
    
    % 6: age specific
    lambda_p = [0.970972761, 0.970972761, 0.970972761, 0.970972761, 0.756281823, 0.756281823, 0.756281823, 0.756281823, 0.761410488, 0.761410488, 0.761410488, 0.761410488, 0.743674495, 0.743674495, 0.743674495, 0.743674495, 0.739953194, 0.739953194, 0.739953194, 0.739953194, 0.729003481, 0.729003481, 0.729003481, 0.729003481, 0.709925859, 0.709925859, 0.709925859, 0.709925859, 0.713763071, 0.713763071, 0.713763071, 0.713763071, 0.720067898, 0.720067898, 0.720067898, 0.720067898, 0.705672443, 0.705672443, 0.705672443, 0.705672443, 0.665122385, 0.665122385, 0.665122385, 0.665122385, 0.691760699, 0.691760699, 0.691760699, 0.691760699, 0.670526562, 0.670526562, 0.670526562, 0.670526562, 0.666851067, 0.666851067, 0.666851067, 0.666851067, 0.622426222, 0.622426222, 0.622426222, 0.622426222, 0.623560851, 0.623560851, 0.623560851, 0.623560851, 0.602885815, 0.602885815, 0.602885815, 0.602885815, 0.550397061, 0.550397061, 0.550397061, 0.550397061, 0.228296792, 0.228296792, 0.228296792, 0.228296792, 0.228296792, 0.228296792];
    lambda_p(2, :) = [0.618623087, 0.618623087, 0.618623087, 0.618623087, 0.420712624, 0.420712624, 0.420712624, 0.420712624, 0.41798676, 0.41798676, 0.41798676, 0.41798676, 0.385726635, 0.385726635, 0.385726635, 0.385726635, 0.471669121, 0.471669121, 0.471669121, 0.471669121, 0.439686608, 0.439686608, 0.439686608, 0.439686608, 0.351395835, 0.351395835, 0.351395835, 0.351395835, 0.373529995, 0.373529995, 0.373529995, 0.373529995, 0.449737056, 0.449737056, 0.449737056, 0.449737056, 0.476490295, 0.476490295, 0.476490295, 0.476490295, 0.263617725, 0.263617725, 0.263617725, 0.263617725, 0.375752439, 0.375752439, 0.375752439, 0.375752439, 0.401828061, 0.401828061, 0.401828061, 0.401828061, 0.36608758, 0.36608758, 0.36608758, 0.36608758, 0.328501249, 0.328501249, 0.328501249, 0.328501249, 0.448608994, 0.448608994, 0.448608994, 0.448608994, 0.261494831, 0.261494831, 0.261494831, 0.261494831, 0.198201067, 0.198201067, 0.198201067, 0.198201067, 0.057887798, 0.057887798, 0.057887798, 0.057887798, 0.057887798, 0.057887798];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp1_agespec';
    figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp1_agespec';
    figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp1_agespec';
    figures_report
    
    % clears lambdas and turns off hospital constraints
    clear lambda_p_h lambda_p_i lambda_p_r lambda_p_f
    clear S_fig S_table
    Z = 1;
    
    % parameters: tw0, bL
    theta = 0.03327394;
    gamma = 0.63589003;
    lambda_d = 1.5654178;
    lambda_i = 4.502474 - lambda_d;
    alpha1 = 0;
    alpha2 = 0;
    Pi_star = 0.106669183;
    Pi0 = 12.127;
    b = 5.6; % calibrated to fit Swedish data
    
    % benchmark
    equilibrium
    results_table
    
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    % universal
    lambda_p = [0.038716133, 0.038716133, 0.038716133, 0.038716133, 2.67E-09, 2.67E-09, 2.67E-09, 2.67E-09, 7.61E-10, 7.61E-10, 7.61E-10, 7.61E-10, 1.98E-09, 1.98E-09, 1.98E-09, 1.98E-09, 1.52E-08, 1.52E-08, 1.52E-08, 1.52E-08, 0.193301581, 0.193301581, 0.193301581, 0.193301581, 0.370906591, 0.370906591, 0.370906591, 0.370906591, 0.447329999, 0.447329999, 0.447329999, 0.447329999, 0.467615116, 0.467615116, 0.467615116, 0.467615116, 0.459311708, 0.459311708, 0.459311708, 0.459311708, 0.436863652, 0.436863652, 0.436863652, 0.436863652, 0.408179365, 0.408179365, 0.408179365, 0.408179365, 0.376539647, 0.376539647, 0.376539647, 0.376539647, 0.342875951, 0.342875951, 0.342875951, 0.342875951, 0.306523213, 0.306523213, 0.306523213, 0.306523213, 0.266661869, 0.266661869, 0.266661869, 0.266661869, 0.220362759, 0.220362759, 0.220362759, 0.220362759, 0.162871723, 0.162871723, 0.162871723, 0.162871723, 0.057404445, 0.057404445, 0.057404445, 0.057404445, 0.057404445, 0.057404445];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw0_bL_hosp0_univ';
    figures_paper
    filename_figures_smallplot_suffix = 'tw0_bL_hosp0_univ';
    figures_smallplot
    filename_figures_report_suffix = 'tw0_bL_hosp0_univ';
    figures_report
    
    % age specific
    lambda_p = [0.055232141, 0.055232141, 0.055232141, 0.055232141, 8.07E-05, 8.07E-05, 8.07E-05, 8.07E-05, 6.43E-05, 6.43E-05, 6.43E-05, 6.43E-05, 0.000110794, 0.000110794, 0.000110794, 0.000110794, 0.004431451, 0.004431451, 0.004431451, 0.004431451, 0.199556837, 0.199556837, 0.199556837, 0.199556837, 0.396278136, 0.396278136, 0.396278136, 0.396278136, 0.475121581, 0.475121581, 0.475121581, 0.475121581, 0.461982259, 0.461982259, 0.461982259, 0.461982259, 0.467705172, 0.467705172, 0.467705172, 0.467705172, 0.441989902, 0.441989902, 0.441989902, 0.441989902, 0.408709256, 0.408709256, 0.408709256, 0.408709256, 0.372691064, 0.372691064, 0.372691064, 0.372691064, 0.338250319, 0.338250319, 0.338250319, 0.338250319, 0.315506614, 0.315506614, 0.315506614, 0.315506614, 0.274033373, 0.274033373, 0.274033373, 0.274033373, 0.229751093, 0.229751093, 0.229751093, 0.229751093, 0.153083384, 0.153083384, 0.153083384, 0.153083384, 0.04713839, 0.04713839, 0.04713839, 0.04713839, 0.04713839, 0.04713839];
    lambda_p(2, :) = [0.212362146, 0.212362146, 0.212362146, 0.212362146, 0.049440507, 0.049440507, 0.049440507, 0.049440507, 0.327561012, 0.327561012, 0.327561012, 0.327561012, 0.478266865, 0.478266865, 0.478266865, 0.478266865, 0.788892041, 0.788892041, 0.788892041, 0.788892041, 0.52493742, 0.52493742, 0.52493742, 0.52493742, 0.281528385, 0.281528385, 0.281528385, 0.281528385, 0.422136588, 0.422136588, 0.422136588, 0.422136588, 0.540399397, 0.540399397, 0.540399397, 0.540399397, 0.494277207, 0.494277207, 0.494277207, 0.494277207, 0.47238721, 0.47238721, 0.47238721, 0.47238721, 0.564035957, 0.564035957, 0.564035957, 0.564035957, 0.414417054, 0.414417054, 0.414417054, 0.414417054, 0.52726148, 0.52726148, 0.52726148, 0.52726148, 0.4793209, 0.4793209, 0.4793209, 0.4793209, 0.310711425, 0.310711425, 0.310711425, 0.310711425, 0.290083054, 0.290083054, 0.290083054, 0.290083054, 0.027365673, 0.027365673, 0.027365673, 0.027365673, 0.16815866, 0.16815866, 0.16815866, 0.16815866, 0.16815866, 0.16815866];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw0_bL_hosp0_agespec';
    figures_paper
    filename_figures_smallplot_suffix = 'tw0_bL_hosp0_agespec';
    figures_smallplot
    filename_figures_report_suffix = 'tw0_bL_hosp0_agespec';
    figures_report
    
elseif do==6 % vaccine takes a long time to arrive
    
    filename_results_table = 'tw1_bH_vacc5_10_15';
    
    t_vaccine_vec = [5, 10, 15];
    
    for i_t_vaccine = 1:length(t_vaccine_vec)
        
        t_vaccine = t_vaccine_vec(i_t_vaccine) * 52;
        T = t_vaccine + 10;
        
        disp(i_t_vaccine)
        
        clear S_table S_Pi S_c_eq
        clear lambda_p_i xi_p
        
        % benchmark
        equilibrium
        results_table

        S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
        S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);
        S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);

        % epidemiological
        flag_epidemiological = 1;

        equilibrium
        results_table

        flag_epidemiological = 0;

        % age ext. partial
        flag_Pi = 1;
        flag_fake_young = 1;

        equilibrium
        results_table

        flag_Pi = 0;
        flag_fake_young = 0;

        % age ext. general
        flag_fake_young = 1;

        equilibrium
        results_table

        flag_fake_young = 0;
        
    end
    
%     t_vaccine = 15 * 52;
%     T = t_vaccine + 10;
%     
%     filename_results_table = 'tw1_bH_vacc15';
%     
%     % benchmark
%     equilibrium
%     results_table
%     
%     S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
%     S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);
%     S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
%     
%     % epidemiological
%     flag_epidemiological = 1;
%     
%     equilibrium
%     results_table
%     
%     flag_epidemiological = 0;
%     
%     % age ext. partial
%     flag_Pi = 1;
%     flag_fake_young = 1;
%     
%     equilibrium
%     results_table
%     
%     flag_Pi = 0;
%     flag_fake_young = 0;
%     
%     % age ext. general
%     flag_fake_young = 1;
%     
%     equilibrium
%     results_table
%     
%     flag_fake_young = 0;
    
elseif do==7 % non-targeted moments
    
    % runs benchmark
    equilibrium
    
    % leisure in the no-COVID world
    disp('l_h_term * 112 = ')
    disp(l_h_term * 112)
    
    % weekly growth of infections
    M = frac_young * M_i_all(i_young, :) + frac_old * M_i_all(i_old, :);
    weekly_growth = (M(:, 2:T) ./ M(:, 1:T-1) - 1) * 100;
    daily_growth = ((M(:, 2:T) ./ M(:, 1:T-1)).^(1/7) - 1) * 100;
    
    disp('weekly growth of infections (%)')
    disp(weekly_growth(1:10))
    
    disp('daily growth of infections (%)')
    disp(daily_growth(1:10))
    
    % old/all
    M = M_d(i_old, :) * frac_old ./ ...
        (M_d(i_old, :) * frac_old + M_d(i_young, :) * frac_young);
    
    disp('old dead/all dead')
    disp(M(1:10))
    
elseif do==8 % no altruism
    
    filename_results_table = 'tw1_bH_no_altruism';
    lambda_i = 0;
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    % epidemiological
    flag_epidemiological = 1;
    
    equilibrium
    results_table
    
elseif do==9 % Ferguson
    
    filename_results_table = 'tw1_bH_ferguson';
    
    phi(i_nosymptoms, i_young) = 0.991289141;
    phi(i_nosymptoms, i_old) = 0.893010253;
    phi(i_symptoms, i_young) = 0.284090909;
    phi(i_symptoms, i_old) = 0.284090909;
    delta(i_young) = 0.372680985;
    delta(i_old) = 0.371014263;
    Pi0 = 13.7;
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    % epidemiological
    flag_epidemiological = 1;
    
    equilibrium
    results_table
    
elseif do==10 % policies with hospital constraints
    
    filename_results_table = 'tw1_bH_hosp';
    Z = Z_data;
    policies
    
elseif do==11 % no natural deaths for the old
    
    filename_results_table = 'tw1_bH_no_natural_deaths';
    
    % with natural deaths
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    no natural deaths
    survival_young = 1;
    survival_old = 1;

    % survival vector
    Delta = zeros(n_age, 1);
    Delta(i_young) = survival_young;
    Delta(i_old) = survival_old;

    % effective discount (beta)
    beta = time_discount * Delta;
    
    % runs simulation
    equilibrium
    results_table
    
elseif do==12 % partial equilibrium analysis w.r.t. optimal policy
    
    filename_results_table = 'tw1_bH_optimal_partial';
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);

    % optimal: age specific
    lambda_p = [0.970508218, 0.970508218, 0.970508218, 0.970508218, 0.758040022, 0.758040022, 0.758040022, 0.758040022, 0.75806864, 0.75806864, 0.75806864, 0.75806864, 0.745269026, 0.745269026, 0.745269026, 0.745269026, 0.724599608, 0.724599608, 0.724599608, 0.724599608, 0.727676653, 0.727676653, 0.727676653, 0.727676653, 0.713283369, 0.713283369, 0.713283369, 0.713283369, 0.702679432, 0.702679432, 0.702679432, 0.702679432, 0.707093519, 0.707093519, 0.707093519, 0.707093519, 0.690472405, 0.690472405, 0.690472405, 0.690472405, 0.685313329, 0.685313329, 0.685313329, 0.685313329, 0.679425038, 0.679425038, 0.679425038, 0.679425038, 0.667398045, 0.667398045, 0.667398045, 0.667398045, 0.654670493, 0.654670493, 0.654670493, 0.654670493, 0.639193325, 0.639193325, 0.639193325, 0.639193325, 0.615437052, 0.615437052, 0.615437052, 0.615437052, 0.587684866, 0.587684866, 0.587684866, 0.587684866, 0.563533574, 0.563533574, 0.563533574, 0.563533574, 0.232571203, 0.232571203, 0.232571203, 0.232571203, 0.232571203, 0.232571203];
    lambda_p(2, :) = [0.664120302, 0.664120302, 0.664120302, 0.664120302, 0.471075544, 0.471075544, 0.471075544, 0.471075544, 0.463235167, 0.463235167, 0.463235167, 0.463235167, 0.533073818, 0.533073818, 0.533073818, 0.533073818, 0.451374354, 0.451374354, 0.451374354, 0.451374354, 0.605783386, 0.605783386, 0.605783386, 0.605783386, 0.417647161, 0.417647161, 0.417647161, 0.417647161, 0.421474431, 0.421474431, 0.421474431, 0.421474431, 0.471715951, 0.471715951, 0.471715951, 0.471715951, 0.447568956, 0.447568956, 0.447568956, 0.447568956, 0.507029106, 0.507029106, 0.507029106, 0.507029106, 0.480841769, 0.480841769, 0.480841769, 0.480841769, 0.436940865, 0.436940865, 0.436940865, 0.436940865, 0.392864232, 0.392864232, 0.392864232, 0.392864232, 0.397521824, 0.397521824, 0.397521824, 0.397521824, 0.368621428, 0.368621428, 0.368621428, 0.368621428, 0.344828721, 0.344828721, 0.344828721, 0.344828721, 0.185175603, 0.185175603, 0.185175603, 0.185175603, 0.171600211, 0.171600211, 0.171600211, 0.171600211, 0.171600211, 0.171600211];

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end

    equilibrium
    results_table
    
    S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);

    % partial eq.
    clear lambda_p_i lambda_p_h lambda_p_r lambda_p_f
    flag_Pi = 1;
    
    equilibrium
    results_table
    
    filename_figures_smallplot_suffix = 'tw1_bH_optimal_partial';
    figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_optimal_partial';
    figures_report
    
elseif do==13 % late lockdown
    
    filename_results_table = 'tw1_bH_lateSH';
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);

    % late lockdown
    new_lambda = 32.3598330005704; % 90% increase in time at home
    duration = 26;
    
    lambda_p_i(:,:) = 0;
    lambda_p_h(:,:) = 0;
    lambda_p_r(:,:) = 0;
    lambda_p_f(:,:) = 0;
    
    lambda_p_i(:, 9:9+duration) = new_lambda;
    lambda_p_h(:, 9:9+duration) = new_lambda;
    lambda_p_r(:, 9:9+duration) = new_lambda;
    lambda_p_f(:, 9:9+duration) = new_lambda;

    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_lateSH';
    figures_paper
    
elseif do==14 % simulates some optimal lockdowns (1 week)
    
    filename_results_table = 'optimal_1week';
    
    % 1: benchmark
    equilibrium
    results_table
    
    S_fig = struct('M_i_all', M_i_all, 'M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp);
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
    
    % 2: universal
    lambda_p = [1.58521406, 0.803243892, 0.791357832, 0.773862409, 0.762806053, 0.755595071, 0.750525764, 0.746236633, 0.742561398, 0.739235146, 0.73630333, 0.73323226, 0.73050359, 0.727774969, 0.72525351, 0.722710109, 0.720342685, 0.717981249, 0.715832832, 0.713606808, 0.711450296, 0.709363588, 0.707385647, 0.705378178, 0.703366036, 0.701538823, 0.699674097, 0.697715893, 0.69592753, 0.694018593, 0.692159155, 0.690413679, 0.688695398, 0.686735588, 0.6849916, 0.683150101, 0.681363086, 0.679506398, 0.677553648, 0.675663457, 0.673732699, 0.671684859, 0.669768145, 0.66751448, 0.665557762, 0.663289617, 0.661012704, 0.658670236, 0.65607989, 0.6536066, 0.65088091, 0.648064367, 0.6451501, 0.641840787, 0.638485476, 0.634900788, 0.631083748, 0.626811623, 0.622439426, 0.617521492, 0.612041499, 0.606174338, 0.599575143, 0.592285647, 0.583831432, 0.574316138, 0.563256679, 0.550509428, 0.535189931, 0.516838099, 0.49408114, 0.46564904, 0.428821987, 0.380030714, 0.31338532, 0.218081334, 0.074851168, 0.000162362];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp0_univ_1week';
    %figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp0_univ_1week';
    %figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp0_univ_1week';
    %figures_report
    
    % 3: age specific
%     lambda_p = [1.85189353, 0.710260839, 0.784422018, 0.770126339, 0.765972731, 0.757489474, 0.754329929, 0.751925572, 0.749947061, 0.744479374, 0.741131927, 0.739549884, 0.737751999, 0.735749025, 0.729545076, 0.73010245, 0.726835611, 0.724278411, 0.721390013, 0.718454504, 0.720203424, 0.717472274, 0.713657772, 0.712957255, 0.709620711, 0.708810279, 0.703385991, 0.716960661, 0.702069139, 0.698627192, 0.703185466, 0.691728894, 0.694937895, 0.694991816, 0.692492661, 0.689748432, 0.691312458, 0.687094477, 0.683797448, 0.685838356, 0.682347713, 0.676643556, 0.677348222, 0.674296704, 0.673143251, 0.669059536, 0.672648194, 0.666433021, 0.663281658, 0.659480844, 0.659004625, 0.654230641, 0.653373365, 0.648466147, 0.646577793, 0.643536504, 0.636324105, 0.633832768, 0.629052859, 0.623326729, 0.618749581, 0.614750626, 0.605788975, 0.601102853, 0.588717991, 0.57855993, 0.570192334, 0.555564925, 0.539392303, 0.527771931, 0.499313026, 0.469168225, 0.433385021, 0.386913678, 0.317710597, 0.219899414, 0.073950378, 0.005520763];
%     lambda_p(2, :) = [0.92474025, 0.722329404, 0.536398696, 0.556147817, 0.548602006, 0.544890789, 0.540581899, 0.534974595, 0.533525876, 0.532161585, 0.526310779, 0.528721275, 0.527866515, 0.530300708, 0.520625566, 0.517886676, 0.522194087, 0.514882463, 0.511902372, 0.513152898, 0.51445859, 0.509480049, 0.510945631, 0.512741282, 0.503663247, 0.496268049, 0.49941222, 0.506082774, 0.499331752, 0.494732209, 0.500425002, 0.493873097, 0.492008834, 0.493787191, 0.489885347, 0.493827639, 0.485540047, 0.486225658, 0.482968468, 0.48254652, 0.479737594, 0.480293124, 0.477779836, 0.478245634, 0.470951501, 0.472668831, 0.47398343, 0.468649132, 0.46443955, 0.467468631, 0.461757845, 0.461316201, 0.457269283, 0.459221205, 0.448889459, 0.452028895, 0.440146159, 0.443815262, 0.437516956, 0.432718372, 0.434395423, 0.424166716, 0.420393169, 0.410899612, 0.413910438, 0.40020856, 0.395329069, 0.375924693, 0.364063283, 0.350210783, 0.323490536, 0.307162348, 0.275574735, 0.244215925, 0.187819469, 0.113806584, 0.000406882, 0.029355847];
    lambda_p = [1.8185434, 0.721996808, 0.787572037, 0.772584232, 0.765966776, 0.760328664, 0.756039114, 0.752209708, 0.748944961, 0.745696235, 0.742649625, 0.739778431, 0.737201439, 0.734600472, 0.73208844, 0.729548129, 0.727391598, 0.725000243, 0.722839001, 0.720694865, 0.71858115, 0.71659909, 0.714529672, 0.712540563, 0.710659096, 0.70874489, 0.706861036, 0.705062538, 0.703228611, 0.701433716, 0.699477694, 0.697857026, 0.695873246, 0.69420826, 0.692325899, 0.690583692, 0.688708425, 0.686718328, 0.684999398, 0.68309028, 0.681053797, 0.679043621, 0.677043247, 0.674915661, 0.672787718, 0.670557056, 0.668348371, 0.66590882, 0.663387502, 0.660832753, 0.658085677, 0.655221465, 0.652138508, 0.648971507, 0.64547989, 0.641981221, 0.638101932, 0.633799949, 0.629290811, 0.624340101, 0.618853333, 0.612903519, 0.606240069, 0.598786802, 0.590388998, 0.580740118, 0.569594804, 0.556662515, 0.541214787, 0.522486452, 0.499558587, 0.470798267, 0.433595147, 0.384338265, 0.316832029, 0.220263343, 0.075109494, 0.000462144];
    lambda_p(2, :) = [0.909976433, 0.725530255, 0.543798508, 0.560969444, 0.552320156, 0.54805362, 0.544507129, 0.541581659, 0.538800718, 0.536255559, 0.53402191, 0.531687497, 0.52956849, 0.527419386, 0.525484733, 0.523460598, 0.521737908, 0.519949821, 0.518059929, 0.516538761, 0.514770142, 0.513124641, 0.51162695, 0.510132948, 0.508275225, 0.507037358, 0.505417864, 0.504021812, 0.502547181, 0.501008453, 0.499618856, 0.498320909, 0.496707433, 0.49530895, 0.49389689, 0.492606866, 0.490895703, 0.489574252, 0.48805453, 0.486497455, 0.484884904, 0.483353416, 0.481674705, 0.480068284, 0.478366163, 0.476686785, 0.474737876, 0.473029785, 0.470850596, 0.469038642, 0.466542684, 0.464473261, 0.461906946, 0.459358242, 0.45676878, 0.453950646, 0.450839167, 0.44746343, 0.443861746, 0.440012504, 0.435439529, 0.430864917, 0.425474734, 0.419772866, 0.413023502, 0.405340418, 0.396739277, 0.386482707, 0.374560802, 0.359620089, 0.341624153, 0.319378491, 0.290511519, 0.253165333, 0.203517662, 0.136420351, 0.05041662, 0.002132444];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    filename_figures_paper_suffix = 'tw1_bH_hosp0_agespec_1week';
    %figures_paper
    filename_figures_smallplot_suffix = 'tw1_bH_hosp0_agespec_1week';
    %figures_smallplot
    filename_figures_report_suffix = 'tw1_bH_hosp0_agespec_1week';
    %figures_report
    
elseif do==15 % optimal policy (4 week) without teleworking
    
    filename_results_table = 'tw0_bH_optimal_4week';
    
    alpha1 = 0;
    alpha2 = 0;
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
    
    % universal
    lambda_p = [2.49096675, 2.49096675, 2.49096675, 2.49096675, 2.17577774, 2.17577774, 2.17577774, 2.17577774, 2.15974947, 2.15974947, 2.15974947, 2.15974947, 2.14399225, 2.14399225, 2.14399225, 2.14399225, 2.12909973, 2.12909973, 2.12909973, 2.12909973, 2.11459644, 2.11459644, 2.11459644, 2.11459644, 2.10069159, 2.10069159, 2.10069159, 2.10069159, 2.08490396, 2.08490396, 2.08490396, 2.08490396, 2.06910058, 2.06910058, 2.06910058, 2.06910058, 2.05147756, 2.05147756, 2.05147756, 2.05147756, 2.03251737, 2.03251737, 2.03251737, 2.03251737, 2.00973293, 2.00973293, 2.00973293, 2.00973293, 1.9833065, 1.9833065, 1.9833065, 1.9833065, 1.95084288, 1.95084288, 1.95084288, 1.95084288, 1.90767423, 1.90767423, 1.90767423, 1.90767423, 1.84917198, 1.84917198, 1.84917198, 1.84917198, 1.76061859, 1.76061859, 1.76061859, 1.76061859, 1.62284328, 1.62284328, 1.62284328, 1.62284328, 0.68204571, 0.68204571, 0.68204571, 0.68204571, 0.68204571, 0.68204571];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    % age specific
    lambda_p = [2.55749457, 2.55749457, 2.55749457, 2.55749457, 2.22981048, 2.22981048, 2.22981048, 2.22981048, 2.21256824, 2.21256824, 2.21256824, 2.21256824, 2.19599454, 2.19599454, 2.19599454, 2.19599454, 2.1799321, 2.1799321, 2.1799321, 2.1799321, 2.1643316, 2.1643316, 2.1643316, 2.1643316, 2.14870382, 2.14870382, 2.14870382, 2.14870382, 2.13273008, 2.13273008, 2.13273008, 2.13273008, 2.11635604, 2.11635604, 2.11635604, 2.11635604, 2.09713241, 2.09713241, 2.09713241, 2.09713241, 2.07710352, 2.07710352, 2.07710352, 2.07710352, 2.05456443, 2.05456443, 2.05456443, 2.05456443, 2.02637353, 2.02637353, 2.02637353, 2.02637353, 1.99232022, 1.99232022, 1.99232022, 1.99232022, 1.94899768, 1.94899768, 1.94899768, 1.94899768, 1.88810585, 1.88810585, 1.88810585, 1.88810585, 1.79702443, 1.79702443, 1.79702443, 1.79702443, 1.65627909, 1.65627909, 1.65627909, 1.65627909, 0.692561497, 0.692561497, 0.692561497, 0.692561497, 0.692561497, 0.692561497];
    lambda_p(2, :) = [1.45378229, 1.45378229, 1.45378229, 1.45378229, 1.29595446, 1.29595446, 1.29595446, 1.29595446, 1.2828639, 1.2828639, 1.2828639, 1.2828639, 1.27377331, 1.27377331, 1.27377331, 1.27377331, 1.26459638, 1.26459638, 1.26459638, 1.26459638, 1.25655435, 1.25655435, 1.25655435, 1.25655435, 1.24925891, 1.24925891, 1.24925891, 1.24925891, 1.24064036, 1.24064036, 1.24064036, 1.24064036, 1.23091213, 1.23091213, 1.23091213, 1.23091213, 1.21883531, 1.21883531, 1.21883531, 1.21883531, 1.2072884, 1.2072884, 1.2072884, 1.2072884, 1.19437258, 1.19437258, 1.19437258, 1.19437258, 1.17898486, 1.17898486, 1.17898486, 1.17898486, 1.1584335, 1.1584335, 1.1584335, 1.1584335, 1.13315063, 1.13315063, 1.13315063, 1.13315063, 1.09620112, 1.09620112, 1.09620112, 1.09620112, 1.04060585, 1.04060585, 1.04060585, 1.04060585, 0.945027098, 0.945027098, 0.945027098, 0.945027098, 0.38259806, 0.38259806, 0.38259806, 0.38259806, 0.38259806, 0.38259806];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
elseif do==16 % figures_optimal (tw1_bH)
    
    filename_figures_optimal_suffix = 'tw1_bH';
    C_fig = cell(2, 4);

    % benchmark
    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 4} = time_outside_hours;

    % age-specific optimal
%     lambda_p = [1.85189353, 0.710260839, 0.784422018, 0.770126339, 0.765972731, 0.757489474, 0.754329929, 0.751925572, 0.749947061, 0.744479374, 0.741131927, 0.739549884, 0.737751999, 0.735749025, 0.729545076, 0.73010245, 0.726835611, 0.724278411, 0.721390013, 0.718454504, 0.720203424, 0.717472274, 0.713657772, 0.712957255, 0.709620711, 0.708810279, 0.703385991, 0.716960661, 0.702069139, 0.698627192, 0.703185466, 0.691728894, 0.694937895, 0.694991816, 0.692492661, 0.689748432, 0.691312458, 0.687094477, 0.683797448, 0.685838356, 0.682347713, 0.676643556, 0.677348222, 0.674296704, 0.673143251, 0.669059536, 0.672648194, 0.666433021, 0.663281658, 0.659480844, 0.659004625, 0.654230641, 0.653373365, 0.648466147, 0.646577793, 0.643536504, 0.636324105, 0.633832768, 0.629052859, 0.623326729, 0.618749581, 0.614750626, 0.605788975, 0.601102853, 0.588717991, 0.57855993, 0.570192334, 0.555564925, 0.539392303, 0.527771931, 0.499313026, 0.469168225, 0.433385021, 0.386913678, 0.317710597, 0.219899414, 0.073950378, 0.005520763];
%     lambda_p(2, :) = [0.92474025, 0.722329404, 0.536398696, 0.556147817, 0.548602006, 0.544890789, 0.540581899, 0.534974595, 0.533525876, 0.532161585, 0.526310779, 0.528721275, 0.527866515, 0.530300708, 0.520625566, 0.517886676, 0.522194087, 0.514882463, 0.511902372, 0.513152898, 0.51445859, 0.509480049, 0.510945631, 0.512741282, 0.503663247, 0.496268049, 0.49941222, 0.506082774, 0.499331752, 0.494732209, 0.500425002, 0.493873097, 0.492008834, 0.493787191, 0.489885347, 0.493827639, 0.485540047, 0.486225658, 0.482968468, 0.48254652, 0.479737594, 0.480293124, 0.477779836, 0.478245634, 0.470951501, 0.472668831, 0.47398343, 0.468649132, 0.46443955, 0.467468631, 0.461757845, 0.461316201, 0.457269283, 0.459221205, 0.448889459, 0.452028895, 0.440146159, 0.443815262, 0.437516956, 0.432718372, 0.434395423, 0.424166716, 0.420393169, 0.410899612, 0.413910438, 0.40020856, 0.395329069, 0.375924693, 0.364063283, 0.350210783, 0.323490536, 0.307162348, 0.275574735, 0.244215925, 0.187819469, 0.113806584, 0.000406882, 0.029355847];
    lambda_p = [1.8185434, 0.721996808, 0.787572037, 0.772584232, 0.765966776, 0.760328664, 0.756039114, 0.752209708, 0.748944961, 0.745696235, 0.742649625, 0.739778431, 0.737201439, 0.734600472, 0.73208844, 0.729548129, 0.727391598, 0.725000243, 0.722839001, 0.720694865, 0.71858115, 0.71659909, 0.714529672, 0.712540563, 0.710659096, 0.70874489, 0.706861036, 0.705062538, 0.703228611, 0.701433716, 0.699477694, 0.697857026, 0.695873246, 0.69420826, 0.692325899, 0.690583692, 0.688708425, 0.686718328, 0.684999398, 0.68309028, 0.681053797, 0.679043621, 0.677043247, 0.674915661, 0.672787718, 0.670557056, 0.668348371, 0.66590882, 0.663387502, 0.660832753, 0.658085677, 0.655221465, 0.652138508, 0.648971507, 0.64547989, 0.641981221, 0.638101932, 0.633799949, 0.629290811, 0.624340101, 0.618853333, 0.612903519, 0.606240069, 0.598786802, 0.590388998, 0.580740118, 0.569594804, 0.556662515, 0.541214787, 0.522486452, 0.499558587, 0.470798267, 0.433595147, 0.384338265, 0.316832029, 0.220263343, 0.075109494, 0.000462144];
    lambda_p(2, :) = [0.909976433, 0.725530255, 0.543798508, 0.560969444, 0.552320156, 0.54805362, 0.544507129, 0.541581659, 0.538800718, 0.536255559, 0.53402191, 0.531687497, 0.52956849, 0.527419386, 0.525484733, 0.523460598, 0.521737908, 0.519949821, 0.518059929, 0.516538761, 0.514770142, 0.513124641, 0.51162695, 0.510132948, 0.508275225, 0.507037358, 0.505417864, 0.504021812, 0.502547181, 0.501008453, 0.499618856, 0.498320909, 0.496707433, 0.49530895, 0.49389689, 0.492606866, 0.490895703, 0.489574252, 0.48805453, 0.486497455, 0.484884904, 0.483353416, 0.481674705, 0.480068284, 0.478366163, 0.476686785, 0.474737876, 0.473029785, 0.470850596, 0.469038642, 0.466542684, 0.464473261, 0.461906946, 0.459358242, 0.45676878, 0.453950646, 0.450839167, 0.44746343, 0.443861746, 0.440012504, 0.435439529, 0.430864917, 0.425474734, 0.419772866, 0.413023502, 0.405340418, 0.396739277, 0.386482707, 0.374560802, 0.359620089, 0.341624153, 0.319378491, 0.290511519, 0.253165333, 0.203517662, 0.136420351, 0.05041662, 0.002132444];

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 1} = time_outside_hours;
    C_fig{1, 1} = lambda_p;

    S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);

    % age specific (partial equilibrium)
    clear lambda_p_i lambda_p_h lambda_p_r lambda_p_f

    flag_Pi = 1;

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 3} = time_outside_hours;

    flag_Pi = 0;

    % universal optimal
    lambda_p = [1.58521406, 0.803243892, 0.791357832, 0.773862409, 0.762806053, 0.755595071, 0.750525764, 0.746236633, 0.742561398, 0.739235146, 0.73630333, 0.73323226, 0.73050359, 0.727774969, 0.72525351, 0.722710109, 0.720342685, 0.717981249, 0.715832832, 0.713606808, 0.711450296, 0.709363588, 0.707385647, 0.705378178, 0.703366036, 0.701538823, 0.699674097, 0.697715893, 0.69592753, 0.694018593, 0.692159155, 0.690413679, 0.688695398, 0.686735588, 0.6849916, 0.683150101, 0.681363086, 0.679506398, 0.677553648, 0.675663457, 0.673732699, 0.671684859, 0.669768145, 0.66751448, 0.665557762, 0.663289617, 0.661012704, 0.658670236, 0.65607989, 0.6536066, 0.65088091, 0.648064367, 0.6451501, 0.641840787, 0.638485476, 0.634900788, 0.631083748, 0.626811623, 0.622439426, 0.617521492, 0.612041499, 0.606174338, 0.599575143, 0.592285647, 0.583831432, 0.574316138, 0.563256679, 0.550509428, 0.535189931, 0.516838099, 0.49408114, 0.46564904, 0.428821987, 0.380030714, 0.31338532, 0.218081334, 0.074851168, 0.000162362];
    lambda_p = repmat(lambda_p, [2, 1]);

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 2} = time_outside_hours;
    C_fig{1, 2} = lambda_p(1, :);

    % makes figures
    figures_optimal
    
elseif do==17 % figures_optimal (tw0_bH)
    
    alpha1 = 0;
    alpha2 = 0;
    filename_figures_optimal_suffix = 'tw0_bH';
    C_fig = cell(2, 4);

    % benchmark
    equilibrium
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    
    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 4} = time_outside_hours;

    % age-specific optimal
    lambda_p = [2.55749457, 2.55749457, 2.55749457, 2.55749457, 2.22981048, 2.22981048, 2.22981048, 2.22981048, 2.21256824, 2.21256824, 2.21256824, 2.21256824, 2.19599454, 2.19599454, 2.19599454, 2.19599454, 2.1799321, 2.1799321, 2.1799321, 2.1799321, 2.1643316, 2.1643316, 2.1643316, 2.1643316, 2.14870382, 2.14870382, 2.14870382, 2.14870382, 2.13273008, 2.13273008, 2.13273008, 2.13273008, 2.11635604, 2.11635604, 2.11635604, 2.11635604, 2.09713241, 2.09713241, 2.09713241, 2.09713241, 2.07710352, 2.07710352, 2.07710352, 2.07710352, 2.05456443, 2.05456443, 2.05456443, 2.05456443, 2.02637353, 2.02637353, 2.02637353, 2.02637353, 1.99232022, 1.99232022, 1.99232022, 1.99232022, 1.94899768, 1.94899768, 1.94899768, 1.94899768, 1.88810585, 1.88810585, 1.88810585, 1.88810585, 1.79702443, 1.79702443, 1.79702443, 1.79702443, 1.65627909, 1.65627909, 1.65627909, 1.65627909, 0.692561497, 0.692561497, 0.692561497, 0.692561497, 0.692561497, 0.692561497];
    lambda_p(2, :) = [1.45378229, 1.45378229, 1.45378229, 1.45378229, 1.29595446, 1.29595446, 1.29595446, 1.29595446, 1.2828639, 1.2828639, 1.2828639, 1.2828639, 1.27377331, 1.27377331, 1.27377331, 1.27377331, 1.26459638, 1.26459638, 1.26459638, 1.26459638, 1.25655435, 1.25655435, 1.25655435, 1.25655435, 1.24925891, 1.24925891, 1.24925891, 1.24925891, 1.24064036, 1.24064036, 1.24064036, 1.24064036, 1.23091213, 1.23091213, 1.23091213, 1.23091213, 1.21883531, 1.21883531, 1.21883531, 1.21883531, 1.2072884, 1.2072884, 1.2072884, 1.2072884, 1.19437258, 1.19437258, 1.19437258, 1.19437258, 1.17898486, 1.17898486, 1.17898486, 1.17898486, 1.1584335, 1.1584335, 1.1584335, 1.1584335, 1.13315063, 1.13315063, 1.13315063, 1.13315063, 1.09620112, 1.09620112, 1.09620112, 1.09620112, 1.04060585, 1.04060585, 1.04060585, 1.04060585, 0.945027098, 0.945027098, 0.945027098, 0.945027098, 0.38259806, 0.38259806, 0.38259806, 0.38259806, 0.38259806, 0.38259806];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 1} = time_outside_hours;
    C_fig{1, 1} = lambda_p;

    S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);

    % age specific (partial equilibrium)
    clear lambda_p_i lambda_p_h lambda_p_r lambda_p_f

    flag_Pi = 1;

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 3} = time_outside_hours;

    flag_Pi = 0;

    % universal optimal
    lambda_p = [2.49096675, 2.49096675, 2.49096675, 2.49096675, 2.17577774, 2.17577774, 2.17577774, 2.17577774, 2.15974947, 2.15974947, 2.15974947, 2.15974947, 2.14399225, 2.14399225, 2.14399225, 2.14399225, 2.12909973, 2.12909973, 2.12909973, 2.12909973, 2.11459644, 2.11459644, 2.11459644, 2.11459644, 2.10069159, 2.10069159, 2.10069159, 2.10069159, 2.08490396, 2.08490396, 2.08490396, 2.08490396, 2.06910058, 2.06910058, 2.06910058, 2.06910058, 2.05147756, 2.05147756, 2.05147756, 2.05147756, 2.03251737, 2.03251737, 2.03251737, 2.03251737, 2.00973293, 2.00973293, 2.00973293, 2.00973293, 1.9833065, 1.9833065, 1.9833065, 1.9833065, 1.95084288, 1.95084288, 1.95084288, 1.95084288, 1.90767423, 1.90767423, 1.90767423, 1.90767423, 1.84917198, 1.84917198, 1.84917198, 1.84917198, 1.76061859, 1.76061859, 1.76061859, 1.76061859, 1.62284328, 1.62284328, 1.62284328, 1.62284328, 0.68204571, 0.68204571, 0.68204571, 0.68204571, 0.68204571, 0.68204571];
    lambda_p = repmat(lambda_p, [2, 1]);

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end

    equilibrium

    time_outside_hours = 112 * (n_h + l_h);

    C_fig{2, 2} = time_outside_hours;
    C_fig{1, 2} = lambda_p(1, :);

    % makes figures
    figures_optimal
    
elseif do==18 % table: optimal policy with testing and quarantine
    
    filename_results_table = 'tw1_bH_optimal_TQ';
    
    % benchmark
    equilibrium
    results_table
    
    S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
    S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
    
    % universal
    xi_p = 0.5 * ones(n_age, 1);
    new_lambda = 32.3598330005704; % 90% lockdown for those who test positive
    
    lambda_p = [0.677907102, 0.44524075, 0.376632654, 0.34439994, 0.325874651, 0.314249987, 0.306069238, 0.299835859, 0.294912341, 0.290623122, 0.286610596, 0.28293704, 0.279615088, 0.276265675, 0.273597159, 0.270577998, 0.267901042, 0.26528934, 0.262852237, 0.260232385, 0.257942731, 0.256014695, 0.253480605, 0.251146616, 0.249263458, 0.247273578, 0.24537583, 0.24333524, 0.241488605, 0.239560609, 0.237901922, 0.236069546, 0.234098092, 0.232334095, 0.230602233, 0.228828736, 0.22705933, 0.225406858, 0.223407564, 0.22184006, 0.219703843, 0.218361277, 0.216333996, 0.214422132, 0.212530716, 0.210674677, 0.208636687, 0.206793911, 0.204486082, 0.202261735, 0.200210412, 0.198056635, 0.19557535, 0.19320853, 0.190529329, 0.188276066, 0.18541285, 0.182217576, 0.179227473, 0.175792942, 0.172585817, 0.168732768, 0.164449348, 0.160527698, 0.155491613, 0.150474637, 0.144705804, 0.138469364, 0.131492262, 0.123639254, 0.114400902, 0.104569441, 0.092041664, 0.078098883, 0.0609482, 0.04033551, 0.014384263, 0.001477183];
    lambda_p = repmat(lambda_p, [2, 1]);
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = max(new_lambda - lambda_i, 0);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table
    
    % age specific
%     lambda_p = [0.708477017, 0.44032757, 0.382379122, 0.352092413, 0.324788106, 0.318200148, 0.310064414, 0.307293295, 0.289355831, 0.290330055, 0.290745759, 0.290726604, 0.278251508, 0.272709189, 0.280422779, 0.27498344, 0.271579059, 0.267530614, 0.270044485, 0.260634797, 0.259320148, 0.258005054, 0.259426388, 0.255820297, 0.252110267, 0.251927418, 0.247028014, 0.24684321, 0.243376499, 0.240148824, 0.24110557, 0.238306642, 0.234806807, 0.231031929, 0.235516258, 0.232637803, 0.229083976, 0.227312823, 0.223329876, 0.222648752, 0.223277132, 0.220540591, 0.217573245, 0.215232136, 0.212829064, 0.213692086, 0.20773464, 0.206791139, 0.206442885, 0.200492918, 0.204978338, 0.201776362, 0.198196579, 0.190350407, 0.193718127, 0.188877056, 0.188071874, 0.182144199, 0.181121232, 0.178093418, 0.174611217, 0.168028561, 0.165709733, 0.163773284, 0.15385048, 0.150308984, 0.147467932, 0.143940555, 0.130697508, 0.122754162, 0.11831891, 0.108647911, 0.086788073, 0.06971503, 0.057791319, 0.046047297, 0.028544487, 0.02770843];
%     lambda_p(2, :) = [0.40421967, 0.31861445, 0.286995104, 0.262191457, 0.264485793, 0.215553015, 0.234882838, 0.245308083, 0.238565674, 0.237505145, 0.246755042, 0.223325191, 0.177860095, 0.168591382, 0.145835253, 0.169592623, 0.190143148, 0.190453295, 0.189871369, 0.174885713, 0.14566414, 0.181257098, 0.172600855, 0.168699409, 0.1385097, 0.147241386, 0.139844016, 0.141418969, 0.203413059, 0.214967074, 0.190236514, 0.178147989, 0.18335087, 0.193523034, 0.153721772, 0.214266761, 0.13368608, 0.176448101, 0.135737802, 0.152007178, 0.167594511, 0.187998338, 0.160177196, 0.163929293, 0.13901804, 0.129076069, 0.134685772, 0.133027451, 0.155112471, 0.159986495, 0.187714488, 0.140485097, 0.170086961, 0.156686232, 0.160033473, 0.159928264, 0.160604616, 0.154658791, 0.176721774, 0.149610658, 0.145707582, 0.172644361, 0.127949875, 0.145344306, 0.067493834, 0.050120194, 0.030115857, 0.04507035, 0.048665612, 0.087546684, 0.07802103, 0.06320248, 0.052060159, 0.011762357, 0.014016165, 0.016104933, 0.038524988, 0.035965378];
    lambda_p = [0.69604005, 0.448676266, 0.38030634, 0.34891018, 0.328206029, 0.317460979, 0.309057913, 0.303520852, 0.29757279, 0.292066863, 0.28934271, 0.285192191, 0.282349481, 0.277714611, 0.276111644, 0.273135552, 0.268670349, 0.2680156, 0.265743532, 0.263267496, 0.260098909, 0.257208711, 0.255864417, 0.255203269, 0.251259748, 0.249164694, 0.246717719, 0.246155355, 0.243824344, 0.242046539, 0.239467057, 0.239361526, 0.235688467, 0.233443044, 0.232740616, 0.229961607, 0.229097288, 0.227708702, 0.225330006, 0.223723914, 0.221265911, 0.22023187, 0.218505054, 0.216552557, 0.213847078, 0.211400693, 0.210833514, 0.208354046, 0.206533188, 0.203338692, 0.202035731, 0.199189852, 0.196696914, 0.195550734, 0.192273385, 0.188989341, 0.18674213, 0.183508319, 0.180820509, 0.177062377, 0.173189488, 0.169865664, 0.165350325, 0.162224074, 0.156711322, 0.150870307, 0.146799158, 0.140336902, 0.131680395, 0.125863474, 0.115178135, 0.105913959, 0.092181454, 0.078482312, 0.05964561, 0.041201912, 0.015888371, 0.002067805];
    lambda_p(2, :) = [0.406750896, 0.32823991, 0.264151093, 0.247629178, 0.232589199, 0.227531673, 0.213978073, 0.211389914, 0.206595409, 0.210250308, 0.203501249, 0.200478472, 0.19943133, 0.20470435, 0.200605449, 0.18669803, 0.190421076, 0.194026275, 0.190981253, 0.194011956, 0.181941156, 0.185585544, 0.182578541, 0.179188368, 0.179504814, 0.177696392, 0.171144817, 0.17035643, 0.175341905, 0.178290143, 0.176711339, 0.168681847, 0.169634194, 0.164710751, 0.172958606, 0.171356823, 0.158874923, 0.162327455, 0.155774814, 0.156952872, 0.161245495, 0.162688078, 0.163432663, 0.154912944, 0.1487629, 0.153162344, 0.145943109, 0.156711854, 0.149220556, 0.150699918, 0.147270623, 0.142718813, 0.137459977, 0.143752817, 0.152507789, 0.129261285, 0.152246551, 0.12118299, 0.132998226, 0.141085503, 0.126485773, 0.13108487, 0.124257482, 0.125472288, 0.083096844, 0.106426372, 0.098260883, 0.121632947, 0.097429055, 0.095787836, 0.095650209, 0.079340598, 0.057743701, 4.37E-06, 0.04359081, 0.02096658, 0.013297797, 0.011349973];

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = max(new_lambda - lambda_i, 0);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    results_table

elseif do==19 % figure: optimal policy with testing and quarantine
    
    filename_figures_optimal_suffix = 'tw1_bH_optimal_TQ';
    
    C_fig = cell(2, 2);
    
    % optimal lockdown with TQ
    xi_p = 0.5 * ones(n_age, 1); % 50% testing
    new_lambda = 32.3598330005704; % 90% lockdown for those who test positive
    
%     lambda_p = [0.708477017, 0.44032757, 0.382379122, 0.352092413, 0.324788106, 0.318200148, 0.310064414, 0.307293295, 0.289355831, 0.290330055, 0.290745759, 0.290726604, 0.278251508, 0.272709189, 0.280422779, 0.27498344, 0.271579059, 0.267530614, 0.270044485, 0.260634797, 0.259320148, 0.258005054, 0.259426388, 0.255820297, 0.252110267, 0.251927418, 0.247028014, 0.24684321, 0.243376499, 0.240148824, 0.24110557, 0.238306642, 0.234806807, 0.231031929, 0.235516258, 0.232637803, 0.229083976, 0.227312823, 0.223329876, 0.222648752, 0.223277132, 0.220540591, 0.217573245, 0.215232136, 0.212829064, 0.213692086, 0.20773464, 0.206791139, 0.206442885, 0.200492918, 0.204978338, 0.201776362, 0.198196579, 0.190350407, 0.193718127, 0.188877056, 0.188071874, 0.182144199, 0.181121232, 0.178093418, 0.174611217, 0.168028561, 0.165709733, 0.163773284, 0.15385048, 0.150308984, 0.147467932, 0.143940555, 0.130697508, 0.122754162, 0.11831891, 0.108647911, 0.086788073, 0.06971503, 0.057791319, 0.046047297, 0.028544487, 0.02770843];
%     lambda_p(2, :) = [0.40421967, 0.31861445, 0.286995104, 0.262191457, 0.264485793, 0.215553015, 0.234882838, 0.245308083, 0.238565674, 0.237505145, 0.246755042, 0.223325191, 0.177860095, 0.168591382, 0.145835253, 0.169592623, 0.190143148, 0.190453295, 0.189871369, 0.174885713, 0.14566414, 0.181257098, 0.172600855, 0.168699409, 0.1385097, 0.147241386, 0.139844016, 0.141418969, 0.203413059, 0.214967074, 0.190236514, 0.178147989, 0.18335087, 0.193523034, 0.153721772, 0.214266761, 0.13368608, 0.176448101, 0.135737802, 0.152007178, 0.167594511, 0.187998338, 0.160177196, 0.163929293, 0.13901804, 0.129076069, 0.134685772, 0.133027451, 0.155112471, 0.159986495, 0.187714488, 0.140485097, 0.170086961, 0.156686232, 0.160033473, 0.159928264, 0.160604616, 0.154658791, 0.176721774, 0.149610658, 0.145707582, 0.172644361, 0.127949875, 0.145344306, 0.067493834, 0.050120194, 0.030115857, 0.04507035, 0.048665612, 0.087546684, 0.07802103, 0.06320248, 0.052060159, 0.011762357, 0.014016165, 0.016104933, 0.038524988, 0.035965378];
    lambda_p = [0.69604005, 0.448676266, 0.38030634, 0.34891018, 0.328206029, 0.317460979, 0.309057913, 0.303520852, 0.29757279, 0.292066863, 0.28934271, 0.285192191, 0.282349481, 0.277714611, 0.276111644, 0.273135552, 0.268670349, 0.2680156, 0.265743532, 0.263267496, 0.260098909, 0.257208711, 0.255864417, 0.255203269, 0.251259748, 0.249164694, 0.246717719, 0.246155355, 0.243824344, 0.242046539, 0.239467057, 0.239361526, 0.235688467, 0.233443044, 0.232740616, 0.229961607, 0.229097288, 0.227708702, 0.225330006, 0.223723914, 0.221265911, 0.22023187, 0.218505054, 0.216552557, 0.213847078, 0.211400693, 0.210833514, 0.208354046, 0.206533188, 0.203338692, 0.202035731, 0.199189852, 0.196696914, 0.195550734, 0.192273385, 0.188989341, 0.18674213, 0.183508319, 0.180820509, 0.177062377, 0.173189488, 0.169865664, 0.165350325, 0.162224074, 0.156711322, 0.150870307, 0.146799158, 0.140336902, 0.131680395, 0.125863474, 0.115178135, 0.105913959, 0.092181454, 0.078482312, 0.05964561, 0.041201912, 0.015888371, 0.002067805];
    lambda_p(2, :) = [0.406750896, 0.32823991, 0.264151093, 0.247629178, 0.232589199, 0.227531673, 0.213978073, 0.211389914, 0.206595409, 0.210250308, 0.203501249, 0.200478472, 0.19943133, 0.20470435, 0.200605449, 0.18669803, 0.190421076, 0.194026275, 0.190981253, 0.194011956, 0.181941156, 0.185585544, 0.182578541, 0.179188368, 0.179504814, 0.177696392, 0.171144817, 0.17035643, 0.175341905, 0.178290143, 0.176711339, 0.168681847, 0.169634194, 0.164710751, 0.172958606, 0.171356823, 0.158874923, 0.162327455, 0.155774814, 0.156952872, 0.161245495, 0.162688078, 0.163432663, 0.154912944, 0.1487629, 0.153162344, 0.145943109, 0.156711854, 0.149220556, 0.150699918, 0.147270623, 0.142718813, 0.137459977, 0.143752817, 0.152507789, 0.129261285, 0.152246551, 0.12118299, 0.132998226, 0.141085503, 0.126485773, 0.13108487, 0.124257482, 0.125472288, 0.083096844, 0.106426372, 0.098260883, 0.121632947, 0.097429055, 0.095787836, 0.095650209, 0.079340598, 0.057743701, 4.37E-06, 0.04359081, 0.02096658, 0.013297797, 0.011349973];
    
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);

    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = max(new_lambda - lambda_i, 0);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    
    time_outside_hours = 112 * (n_h + l_h);
    C_fig{2, 1} = time_outside_hours;
    C_fig{1, 1} = lambda_p;
    
    % optimal lockdown without TQ
    xi_p(:) = 0;
%     lambda_p = [1.85189353, 0.710260839, 0.784422018, 0.770126339, 0.765972731, 0.757489474, 0.754329929, 0.751925572, 0.749947061, 0.744479374, 0.741131927, 0.739549884, 0.737751999, 0.735749025, 0.729545076, 0.73010245, 0.726835611, 0.724278411, 0.721390013, 0.718454504, 0.720203424, 0.717472274, 0.713657772, 0.712957255, 0.709620711, 0.708810279, 0.703385991, 0.716960661, 0.702069139, 0.698627192, 0.703185466, 0.691728894, 0.694937895, 0.694991816, 0.692492661, 0.689748432, 0.691312458, 0.687094477, 0.683797448, 0.685838356, 0.682347713, 0.676643556, 0.677348222, 0.674296704, 0.673143251, 0.669059536, 0.672648194, 0.666433021, 0.663281658, 0.659480844, 0.659004625, 0.654230641, 0.653373365, 0.648466147, 0.646577793, 0.643536504, 0.636324105, 0.633832768, 0.629052859, 0.623326729, 0.618749581, 0.614750626, 0.605788975, 0.601102853, 0.588717991, 0.57855993, 0.570192334, 0.555564925, 0.539392303, 0.527771931, 0.499313026, 0.469168225, 0.433385021, 0.386913678, 0.317710597, 0.219899414, 0.073950378, 0.005520763];
%     lambda_p(2, :) = [0.92474025, 0.722329404, 0.536398696, 0.556147817, 0.548602006, 0.544890789, 0.540581899, 0.534974595, 0.533525876, 0.532161585, 0.526310779, 0.528721275, 0.527866515, 0.530300708, 0.520625566, 0.517886676, 0.522194087, 0.514882463, 0.511902372, 0.513152898, 0.51445859, 0.509480049, 0.510945631, 0.512741282, 0.503663247, 0.496268049, 0.49941222, 0.506082774, 0.499331752, 0.494732209, 0.500425002, 0.493873097, 0.492008834, 0.493787191, 0.489885347, 0.493827639, 0.485540047, 0.486225658, 0.482968468, 0.48254652, 0.479737594, 0.480293124, 0.477779836, 0.478245634, 0.470951501, 0.472668831, 0.47398343, 0.468649132, 0.46443955, 0.467468631, 0.461757845, 0.461316201, 0.457269283, 0.459221205, 0.448889459, 0.452028895, 0.440146159, 0.443815262, 0.437516956, 0.432718372, 0.434395423, 0.424166716, 0.420393169, 0.410899612, 0.413910438, 0.40020856, 0.395329069, 0.375924693, 0.364063283, 0.350210783, 0.323490536, 0.307162348, 0.275574735, 0.244215925, 0.187819469, 0.113806584, 0.000406882, 0.029355847];
    lambda_p = [1.8185434, 0.721996808, 0.787572037, 0.772584232, 0.765966776, 0.760328664, 0.756039114, 0.752209708, 0.748944961, 0.745696235, 0.742649625, 0.739778431, 0.737201439, 0.734600472, 0.73208844, 0.729548129, 0.727391598, 0.725000243, 0.722839001, 0.720694865, 0.71858115, 0.71659909, 0.714529672, 0.712540563, 0.710659096, 0.70874489, 0.706861036, 0.705062538, 0.703228611, 0.701433716, 0.699477694, 0.697857026, 0.695873246, 0.69420826, 0.692325899, 0.690583692, 0.688708425, 0.686718328, 0.684999398, 0.68309028, 0.681053797, 0.679043621, 0.677043247, 0.674915661, 0.672787718, 0.670557056, 0.668348371, 0.66590882, 0.663387502, 0.660832753, 0.658085677, 0.655221465, 0.652138508, 0.648971507, 0.64547989, 0.641981221, 0.638101932, 0.633799949, 0.629290811, 0.624340101, 0.618853333, 0.612903519, 0.606240069, 0.598786802, 0.590388998, 0.580740118, 0.569594804, 0.556662515, 0.541214787, 0.522486452, 0.499558587, 0.470798267, 0.433595147, 0.384338265, 0.316832029, 0.220263343, 0.075109494, 0.000462144];
    lambda_p(2, :) = [0.909976433, 0.725530255, 0.543798508, 0.560969444, 0.552320156, 0.54805362, 0.544507129, 0.541581659, 0.538800718, 0.536255559, 0.53402191, 0.531687497, 0.52956849, 0.527419386, 0.525484733, 0.523460598, 0.521737908, 0.519949821, 0.518059929, 0.516538761, 0.514770142, 0.513124641, 0.51162695, 0.510132948, 0.508275225, 0.507037358, 0.505417864, 0.504021812, 0.502547181, 0.501008453, 0.499618856, 0.498320909, 0.496707433, 0.49530895, 0.49389689, 0.492606866, 0.490895703, 0.489574252, 0.48805453, 0.486497455, 0.484884904, 0.483353416, 0.481674705, 0.480068284, 0.478366163, 0.476686785, 0.474737876, 0.473029785, 0.470850596, 0.469038642, 0.466542684, 0.464473261, 0.461906946, 0.459358242, 0.45676878, 0.453950646, 0.450839167, 0.44746343, 0.443861746, 0.440012504, 0.435439529, 0.430864917, 0.425474734, 0.419772866, 0.413023502, 0.405340418, 0.396739277, 0.386482707, 0.374560802, 0.359620089, 0.341624153, 0.319378491, 0.290511519, 0.253165333, 0.203517662, 0.136420351, 0.05041662, 0.002132444];

    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    equilibrium
    
    time_outside_hours = 112 * (n_h + l_h);
    C_fig{2, 2} = time_outside_hours;
    C_fig{1, 2} = lambda_p;
    
    figures_optimal
    
elseif do==20 % computes CEV of some lockdowns compared to epidemiological
    
    % benchmark
    equilibrium
    
    % epidemiological
    flag_epidemiological = 1;
    equilibrium
    S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
    flag_epidemiological = 0;
    
    % benchmark again
    equilibrium
    
    for i_age = 1:n_age
        factor_c_eq = c_eq(V_h_private, i_age);
        x = (factor_c_eq - 1) * 100;
        disp(x)
    end
    
    % lockdowns
    mild = 0.383675452740587;
    strict = 32.3598330005704;
    duration_vec = [4, 4, 26, 26, 26, 26];
    lambda_vec = [mild, strict, mild, strict, mild, strict];
    who_vec = [1, 1, 1, 1, 2, 2]; % 1: all, 2: old
    
    for i_sims = 1:length(duration_vec)
        disp(i_sims)
        duration = duration_vec(i_sims);
        new_lambda = lambda_vec(i_sims);
        who = who_vec(i_sims);
        
        lambda_p_i(:,:) = 0;
        lambda_p_h(:,:) = 0;
        lambda_p_r(:,:) = 0;
        lambda_p_f(:,:) = 0;
        
        if who==1
            lambda_p_i(:, 1:duration) = max(new_lambda - lambda_i, 0);
            lambda_p_h(:, 1:duration) = new_lambda;
            lambda_p_r(:, 1:duration) = new_lambda;
            lambda_p_f(:, 1:duration) = new_lambda;
        else
            lambda_p_i(i_old, 1:duration) = max(new_lambda - lambda_i, 0);
            lambda_p_h(i_old, 1:duration) = new_lambda;
            lambda_p_r(i_old, 1:duration) = new_lambda;
            lambda_p_f(i_old, 1:duration) = new_lambda;
        end
        
        equilibrium
        
        for i_age = 1:n_age
            factor_c_eq = c_eq(V_h_private, i_age);
            x = (factor_c_eq - 1) * 100;
            disp(x)
        end
        
    end
    
elseif do==21 % social interactions by old in ND and at peak
    
    % benchmark
    equilibrium
    
    % steady state hours
    total_hours = 0;
    
    for i_age = 1:n_age
        total_hours = total_hours + ...
            frac_age(i_age) * (n_h(i_age, end) + l_h(i_age, end));
    end
    
    hours_old = frac_old * (n_h(i_old, end) + l_h(i_old, end));
    
    disp(hours_old / total_hours)
    
    % at peak
    total_hours = 0;
    
    for i_age = 1:n_age
        total_hours = total_hours + ...
            frac_age(i_age) * (n_h(i_age, t_peak_I) + l_h(i_age, t_peak_I));
    end
    
    hours_old = frac_old * (n_h(i_old, t_peak_I) + l_h(i_old, t_peak_I));
    
    disp(hours_old / total_hours)
    
end










%toc

















