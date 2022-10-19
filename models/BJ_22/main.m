
% clc
% clear
% close all

% tic

global N T q q_m1 rho_r_s rho_r_a w delta_n delta_h delta_s delta_l A ...
    beta alpha1 alpha2 gamma rho_d initial_infected varphi Pi_pq Pi_pn ...
    h_over_n_plus_h n_nd m_nd flag_internalize_q Pi_q flag_planner ...
    delta_d planner_obj_n

%-------------------------------------------------------------------------%
%                                 Moments                                 %
%-------------------------------------------------------------------------%

R0_data = 2.5; % R0 for COVID
h_over_n_plus_h = 1 - 0.9823667407; % h/(n+h) in the no-disease case

%-------------------------------------------------------------------------%
%                          Technical parameters                           %
%-------------------------------------------------------------------------%

% T = 120; % time horizon
T = 52 * 3; % time horizon
% T = 52 * 1.5; % time horizon

% fmincon parameters
max_iter = intmax;
tol = 1e-10;

% opts_fmincon = optimoptions('fmincon', ...
%           'Algorithm', 'active-set', ...  
%           'MaxIter', max_iter, ...
%           'MaxFunEvals', max_iter, ...
%           'FunctionTolerance', tol, ...
%           'OptimalityTolerance', tol, ...
%           'StepTolerance', tol, ...
%           'ConstraintTolerance', tol, ...
%           'TolCon', tol, ...
%           'Display', 'None');

%-------------------------------------------------------------------------%
%                               Parameters                                %
%-------------------------------------------------------------------------%

alpha1 = 1; % production function concavity (inside)
alpha2 = 0.7; % production function concavity (outside)
gamma = 0.957; % teleworking productivity

N = 1; % nr. of workers
%initial_infected = 0.001; % initial number of infected workers
helper=load('inf_ini.mat');
initial_infected=helper.helper;
beta = 0.96^(1/52); % time preference
delta_n = ones(T, 1); % relative wage rate
delta_h = 1; % relative wage rate
delta_l = 1; % relative wage rate
delta_s = 1; % relative wage rate
rho_r_s = 0.5; % probability of recovery (symptomatic)
rho_r_a = 0.5; % probability of recovery (asymptomatic)
rho_d = 0.0025;
varphi = 1/2; % fraction of symptomatic
A = ones(T, 1); % tfp
% A = [linspace(1, 0.9, T/25)'; linspace(0.9, 0.9, 2*T/25)'; linspace(0.9, 1, 9*T/25)'; ones(T * 15/25, 1)];
% delta_n = [linspace(1, 5, T/25)'; linspace(5, 5, 2*T/25)'; linspace(5, 1, 9*T/25)'; ones(T * 13/25, 1)];

Pi_pq = 1; % infection rate parameter (related to intercept of p)
Pi_pn = 0.7;%0.75; %1/N * 2/3; % infection rate parameter (slope) (0.45/1.2) 2/3

T_a = 1/rho_r_a;
T_i = 1 + (1-varphi) * T_a;

% Pi_q = (R0_data/((1-initial_infected) * T_i) - Pi_pn)/Pi_pq; % Pi_q to fit R0_data
% fprintf('Pi_q: %.8f\n', Pi_q)

Pi_q = 0.55; %1/4;

Pi2 = Pi_q;

% R0 = (1-initial_infected) * T_i * (Pi_q * Pi_pq + Pi_pn);
R0 = T_i * (Pi_q * Pi_pq + Pi_pn);
fprintf('R0: %.8f\n', R0)

% qq = Pi_q * initial_infected;
% pp = Pi_pq * qq + Pi_pn * initial_infected;
% R0 = T_i * pp/initial_infected % transmissions per initially infected

%-------------------------------------------------------------------------%
%                Sets some parameters by matching moments                 %
%-------------------------------------------------------------------------%

if alpha1==1 % h is zero in the no-disease case
    
    h_over_n_plus_h = 0;
    
    w = A(T)*alpha2/(delta_n(T) * N^(1-alpha2)); % w s.t. optimal demand is N
    
    n_nd = N;
    m_nd = 0;
    
    fprintf('w: %.8f\n', w)
    
else % there is inada for h
    
    % finds w such that n + h = N in no-COVID world
    n_nd = (1-h_over_n_plus_h) * N;
    m_nd = h_over_n_plus_h * N;
    
    w = A(T) * (n_nd^alpha1 + gamma * m_nd^alpha1)^(alpha2-1) * alpha1 * alpha2 * n_nd^(alpha1-1);
    
    fprintf('w: %.8f\n', w)

end

%-------------------------------------------------------------------------%
%                               Some flags                                %
%-------------------------------------------------------------------------%

flag_fixed_choices = 0; % if 1, uses fixed choices for the firm

flag_fix_q = 0; % if 1, uses q_fixed and q_m1_fixed for q and q_m1

flag_endog_A_1 = 0; % if 1, A is a function of the number of sick workers
                    % uses variable A_min

A_min = 0.83;

flag_endog_A_2 = 0; % if 1, A is a function of the number of sick workers
                    % uses slope_A
                    
slope_A = 0;

flag_endog_delta_n = 0; % if 1, delta_n is a function of the n. of sick
                        % workers. Uses variable_delta_n_max

delta_n_max = 1.01;

flag_internalize_q = 0; % if 1, firm internalizes its impact on q

flag_planner = 0; % if 1, uses planner's objective function

delta_d = 0;

flag_use_guess = 0;

planner_obj_n = 1;

%-------------------------------------------------------------------------%
%                              Runs something                             %
%-------------------------------------------------------------------------%

% variable that says what the code will do
do = 1;

if do==1 % finds general equilibrium
    
    % calculates general equilibrium
    equilibrium
    
    % starts a cell
    C = cell(1, 1);
    
    % saves output variables in this cell
    [~, C{1}] = compute_profits(solution);
    
    % makes figures
 %   i_fig = 1;
%     filename_fig_suffix = 'benchmark'; % define this variable if you want to save figure
%  figures_July9
    
elseif do==2 % simulates several scenarios
    
    policies
    
elseif do==3 % some figures
    
    % starts a cell
    C = cell(5, 1);
    
    % equil 1: benchmark
    equilibrium
    [~, C{1}] = compute_profits(solution);
    
    % equil. 2: fixed choices
    flag_fixed_choices = 1;
    equilibrium
    [~, C{2}] = compute_profits(solution);
    flag_fixed_choices = 0;
    
    % equil. 3: delta_l = 0
    delta_l_bm = delta_l;
    delta_l = 0;
    equilibrium
    [~, C{3}] = compute_profits(solution);
    delta_l = delta_l_bm;
    
    % equil. 4: A shock
    A_min = A_min;
    flag_endog_A_1 = 1;
    equilibrium
    [~, C{4}] = compute_profits(solution);
    
    % equil. 5: large A shock
    A_min = 0.5;
    flag_endog_A_1 = 1;
    equilibrium
    [~, C{5}] = compute_profits(solution);
    
    % makes figures
    i_vec = [1, 3, 4, 5];
    i_overlay_vec = [2, 1, 1, 1];
    fig_desc_vec = {'benchmark', 'delta_l', 'tfp', 'tfplarge'};
    
    for ii = 1:length(i_vec)
        for i_overlay = 0:1
            if ii==1
                leg_out = {'Benchmark', 'Fixed choices'};
            elseif ii==4
                leg_out = {'Large changes in $A$', 'Benchmark'};
            else
                leg_out = {'Counterfactual', 'Benchmark'};
            end
            
            if i_overlay==0
                i_fig = i_vec(ii);
                leg_out = '';
            else
                i_fig = [i_vec(ii), i_overlay_vec(ii)];
            end
            
            filename_fig_suffix = [fig_desc_vec{ii}, '_overlay', num2str(i_overlay)];
            figures_July13
            close all
        end
    end
    
elseif do==4 % TFP shock with benchmark q
    
    % starts a cell
    C = cell(3, 1);
    
    % minimum A
    A_min = 0.833;
    flag_endog_A_1 = 0;
    
    % calculates general equilibrium
    equilibrium
    
    % saves benchmark q
    q_fixed = q;
    q_m1_fixed = q_m1;
    
    % saves output variables
    [~, C{1}] = compute_profits(solution);
    
    % TFP shock
    s = C{1}.s;
    s_max = max(s);
    A = 1 + s * (A_min - 1)/s_max;
    
    % fixes benchmark q
    flag_fix_q = 1;
    
    % finds choices of firm
    equilibrium
    
    % saves output variables
    [~, C{2}] = compute_profits(solution);
    
    % sets variables for equilibrium with TFP shock
    flag_fix_q = 0;
    flag_endog_A_1 = 1;
    
    % calculates general equilibrium
    equilibrium
    
    % saves output variables
    [~, C{3}] = compute_profits(solution);
    
    % adds some variables to structs
    for i = 1:3
        C{i}.incubation = C{i}.n_tilde + C{i}.m_tilde;
        C{i}.exposed = [C{i}.incubation(2:T); C{i}.incubation(T)];
        C{i}.susceptible = N - C{i}.s - C{i}.a - C{i}.r_s - C{i}.r_a - C{i}.d - C{i}.incubation - C{i}.exposed;
        C{i}.available = N - C{i}.d - C{i}.s;
    end
    
    % figures 1
    i_fig = [2, 1];
    leg_out = {'Counterfactual', 'Benchmark'};
    t0_fig = 1;
    t1_fig = 60;
    figures_July13
    
    % figures 2
    i_fig = [3, 2];
    leg_out = {'General equilibrium', 'Partial equilibrium'};
    t0_fig = 1;
    t1_fig = 60;
    figures_July13
    
    % figures 3
    i_fig = [3, 1];
    leg_out = {'TFP shock (g.e.)', 'Benchmark'};
    t0_fig = 1;
    t1_fig = 60;
    figures_July13
    
elseif do==5 % policies
    
    policies
    
elseif do==6 % policies with endogenous A
    
    flag_endog_A_1 = 1;
    A_min = 0.8315;
    
    policies
    
elseif do==7 % planner's problem (and the rest)
    
    planner_obj_n = 1;
    
    % vector with delta values
    delta_d_vec = [0, 0, 0, 1e2, 1e3, 1e10];
%     delta_d_vec = [0, 0, 0, 100];
%     delta_d_vec = [701, 702, 703, 704];
    
    % starts a cell
    C = cell(length(delta_d_vec), 1);
    
    % runs planner's problem for many deltas
    for i_delta_d_vec = 1:length(delta_d_vec)
        
        if i_delta_d_vec==1 % benchmark
            flag_internalize_q = 0;
            flag_planner = 0;
        elseif i_delta_d_vec==2 % benchmark, internalize q
            flag_internalize_q = 1;
            flag_planner = 0;
        else
            flag_internalize_q = 1;
            flag_planner = 1;
            
            delta_d = delta_d_vec(i_delta_d_vec);
        end
        
        % solves planner's problem
        equilibrium
        
        % saves output variables in a cell
        [~, C{i_delta_d_vec}] = compute_profits(solution);

        % adds delta_d to C
        C{i_delta_d_vec}.delta_d = delta_d;
        
        % use this solution as guess for next simulation problem
%         flag_use_guess = 1;
    end

    filename = ['mat/planner_obj_n_1_T_', num2str(T), '.mat'];

    save(filename)
    
elseif do==8 % different gamma and delta_n
    

    flag_internalize_q = 0;
    flag_planner = 0;
    
    % vector with delta values
    gamma = linspace(0.70, 0.97, 20);
    delta_n_vec = linspace(1, 1.30, 20);    
    param_vec = combvec(gamma, delta_n_vec);
    
    
    % starts a cell
    C = cell(length(param_vec), 1);
    
    % runs planner's problem for many deltas
    for i_param_vec = 1:length(param_vec)
        
        gamma = param_vec(1, i_param_vec);
        delta_n = param_vec(2, i_param_vec).* ones(T, 1);        
        % solves planner's problem
        equilibrium
        
        % saves output variables in a cell
        [~, C{i_param_vec}] = compute_profits(solution);

        % adds delta_d to C
        C{i_param_vec}.delta_h = delta_h;
        C{i_param_vec}.delta_n = delta_n;        
        
    end    

%     mesh(linspace(0.70, 0.97, 10), delta_n_vec, reshape(X1(14,1:end), [C,C]))

end





 
% toc

%-------------------------------------------------------------------------%
%                                Functions                                %
%-------------------------------------------------------------------------%

% this function evaluates the FOC for choices of the firm in the no-disease
% case
function out = n_h_nd(in)
    
    global A alpha1 alpha2 gamma delta_h delta_n w T

    % "in" is a guess for h

    % transforms in to positive
    h = in^2;
    
    % gets n
    n = (delta_h/(gamma * delta_n))^(1/(1-alpha1)) * h;
    
    % evaluates FOC
    out = A(T) * alpha2 * (n^alpha1 + gamma * h^alpha1)^(alpha2 - 1) * ...
        gamma * alpha1 * h^(alpha1 - 1) - delta_h * w;

end

% this functions finds parameters such that optimal demand of firm in the
% no-disease case fits h/(n+h) = n_over_n_plus_h and n + h = N
function out = calib_nd(in)
    
    global gamma w delta_n delta_h N h_over_n_plus_h alpha1

    % defines parameters
    gamma = 0.5 * (1 + sin(in(1)));
    w = in(2)^2;
    
    % finds optimal choices of n and h
    h = fzero(@n_h_nd, 1);
    h = h^2;
    n = (delta_h/(gamma * delta_n))^(1/(1-alpha1)) * h;
    
    % evaluates conditions to be fitted
    out(1) = h/(n+h) - h_over_n_plus_h;
    out(2) = n + h - N;
    
end