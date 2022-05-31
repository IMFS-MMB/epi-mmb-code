
%-------------------------------------------------------------------------%
%                          Income of old agents                           %
%-------------------------------------------------------------------------%

if alpha1 > 1 % there's teleworking if there's no COVID
    n = n_data;
    v = v_data;
else % there's no teleworking if there's no COVID
    n = total_work_time_data;
    v = 0;
end

income_ND = w * (n + alpha1 * v - alpha2 * v^2);
w_bar = income_old_over_income_young * income_ND; % pension for old folks

%-------------------------------------------------------------------------%
%                            Initialize arrays                            %
%-------------------------------------------------------------------------%

V_h_term = zeros(n_age, 1); % terminal value functions of healthy
V_i_term = zeros(n_age, 1); % terminal value functions of infected
V_s_term = zeros(n_age, 1); % terminal value functions of symptoms

V_h = zeros(n_age, T); % value function for healthy agents, per age
V_f = zeros(n_age, T); % value function for 'fever' agents, per age
V_i = zeros(n_age, T); % value function for infected agents, per age
V_r = zeros(n_age, T); % value function for recovered agents, per age
V_s = zeros(n_age, T); % value function for symptoms agents, per age

V_h_private = zeros(n_age, T); % 'private' values don't use lambda_p
V_f_private = zeros(n_age, T);
V_i_private = zeros(n_age, T);
V_r_private = zeros(n_age, T);
V_s_private = zeros(n_age, T);

V_h_tomorrow = zeros(n_age, 1); % value function in the next period
V_f_tomorrow = zeros(n_age, 1);
V_i_tomorrow = zeros(n_age, 1);
V_r_tomorrow = zeros(n_age, 1);
V_s_tomorrow = zeros(n_age, 1);

V_h_private_tomorrow = zeros(n_age, 1);
V_f_private_tomorrow = zeros(n_age, 1);
V_i_private_tomorrow = zeros(n_age, 1);
V_r_private_tomorrow = zeros(n_age, 1);
V_s_private_tomorrow = zeros(n_age, 1);

n_h = zeros(n_age, T); % policy functions
n_f = zeros(n_age, T);
n_i = zeros(n_age, T);
n_r = zeros(n_age, T);
v_h = zeros(n_age, T);
v_f = zeros(n_age, T);
v_i = zeros(n_age, T);
v_r = zeros(n_age, T);
l_h = zeros(n_age, T);
l_f = zeros(n_age, T);
l_i = zeros(n_age, T);
l_r = zeros(n_age, T);
x_h = zeros(n_age, T);
x_f = zeros(n_age, T);
x_i = zeros(n_age, T);
x_r = zeros(n_age, T);
d_h = zeros(n_age, T);
d_f = zeros(n_age, T);
d_i = zeros(n_age, T);
d_r = zeros(n_age, T);
c_h = zeros(n_age, T);
c_f = zeros(n_age, T);
c_i = zeros(n_age, T);
c_r = zeros(n_age, T);

n_h_case3 = zeros(1, T); % case 3
n_f_case3 = zeros(1, T);
n_i_case3 = zeros(1, T);
n_r_case3 = zeros(1, T);
v_h_case3 = zeros(1, T);
v_f_case3 = zeros(1, T);
v_i_case3 = zeros(1, T);
v_r_case3 = zeros(1, T);
l_h_case3 = zeros(1, T);
l_f_case3 = zeros(1, T);
l_i_case3 = zeros(1, T);
l_r_case3 = zeros(1, T);
x_h_case3 = zeros(1, T);
x_f_case3 = zeros(1, T);
x_i_case3 = zeros(1, T);
x_r_case3 = zeros(1, T);
d_h_case3 = zeros(1, T);
d_f_case3 = zeros(1, T);
d_i_case3 = zeros(1, T);
d_r_case3 = zeros(1, T);
c_h_case3 = zeros(1, T);
c_f_case3 = zeros(1, T);
c_i_case3 = zeros(1, T);
c_r_case3 = zeros(1, T);

c_h_term = zeros(n_age, 1); % choices in the terminal situation
x_h_term = zeros(n_age, 1);
n_h_term = zeros(n_age, 1);
v_h_term = zeros(n_age, 1);
l_h_term = zeros(n_age, 1);
d_h_term = zeros(n_age, 1);

c_i_term = zeros(n_age, 1); % choices in the terminal situation
x_i_term = zeros(n_age, 1);
n_i_term = zeros(n_age, 1);
v_i_term = zeros(n_age, 1);
l_i_term = zeros(n_age, 1);
d_i_term = zeros(n_age, 1);

M_h = zeros(n_age, T); % measure of healthy agents
M_fh = zeros(n_age, T); % measure of healthy 'fever' agents
M_fi = zeros(n_age, T); % measure of infected 'fever' agents
M_i = zeros(n_age, T); % measure of infected agents
M_r = zeros(n_age, T); % measure of recovered agents
M_rn = zeros(n_age, T); % doesn't take into account natural deaths, non-decreasing over time
M_d = zeros(n_age, T); % measure of deceased agents (only covid deaths)
M_dn = zeros(n_age, T); % measure of deceased agents (only natural deaths)
M_t = zeros(n_age, T); % measure of tested agents
N_c_s = zeros(n_age, T); % number of cases (stock)
N_c_f = zeros(n_age, T); % number of cases (flow)

pi_h_vec = zeros(n_age, T-1);
pif_h_vec = zeros(n_age, T-1);
pi_f_vec = zeros(n_age, T-1);
pif_f_vec = zeros(n_age, T-1);

belief_vec = zeros(n_age, T);

delta_vec_new = zeros(n_age, T);

delta_vec = zeros(n_age, T);
delta_vec(i_young, :) = delta_1(i_young);
delta_vec(i_old, :) = delta_1(i_old);

Pi_hat = zeros(n_age, 2); % related to covid infection
I_new = zeros(T, 1); % aggregate time outside of infected agents

% policy variables
if exist('lambda_p_i', 'var')==0
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
end

if exist('xi_p', 'var')==0
    xi_p = zeros(n_age, 1);
end

if exist('factor_epidemiological', 'var')==0
    factor_epidemiological = 1;
end

%-------------------------------------------------------------------------%
%                   Variables in the terminal situation                   %
%-------------------------------------------------------------------------%

% healthy
for i_age = 1:n_age
    
    % solves problem
    L = 0;
    lambda_other = 0;
    flag_fast = 0;
    
    if i_age==i_young
        n_guess = 0;
        v_guess = 0;
        x_guess = 0;
        l_guess = 0;
        [c, x, n, v, l, d] = ...
            solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
            alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    else
        l_guess = 0;
        x_guess = 0;
        [c, x, n, v, l, d] = ...
            solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, ...
            l_guess, x_guess);
    end
    
    % stores choice variables
    c_h_term(i_age) = c;
    x_h_term(i_age) = x;
    n_h_term(i_age) = n;
    v_h_term(i_age) = v;
    l_h_term(i_age) = l;
    d_h_term(i_age) = d;
    
    % computes value
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    
    if flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    V_h_term(i_age) = u/(1-beta(i_age2));
end

% symptoms
for i_age = 1:n_age
    if flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    % there are free hospital beds in terminal situation
    delta_t = delta_1(i_age2);
    
%     V_s_term(i_age) = beta(i_age2) * phi(i_symptoms, i_age2) * v_h_term(i_age) / ...
%         (1 - beta(i_age2) * (1 - phi(i_symptoms, i_age2)) * (1 - delta_t));
    V_s_term(i_age) = beta(i_age2) * phi(i_symptoms, i_age2) * V_h_term(i_age) / ...
        (1 - beta(i_age2) * (1 - phi(i_symptoms, i_age2)) * (1 - delta_t));
end

% infected
for i_age = 1:n_age
    % solves problem
    L = 0;
    lambda_other = lambda_i;
    flag_fast = 0;
    
    if i_age==i_young
        n_guess = 0;
        v_guess = 0;
        x_guess = 0;
        l_guess = 0;
        [c, x, n, v, l, d] = ...
            solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
            alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    else
        l_guess = 0;
        x_guess = 0;
        [c, x, n, v, l, d] = ...
            solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, ...
            l_guess, x_guess);
    end
    
    % stores choice variables
    c_i_term(i_age) = c;
    x_i_term(i_age) = x;
    n_i_term(i_age) = n;
    v_i_term(i_age) = v;
    l_i_term(i_age) = l;
    d_i_term(i_age) = d;
    
    % computes value
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    
    if flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    V_i_term(i_age) = u + beta(i_age2) * phi(i_nosymptoms, i_age2) * V_h_term(i_age) + ...
                       beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_term(i_age);
                            
    V_i_term(i_age) = V_i_term(i_age) / ...
        (1 - beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * (1 - alpha(i_age2)));
end

%-------------------------------------------------------------------------%
%                    Distributions in the first period                    %
%-------------------------------------------------------------------------%
helper=load('inf_ini.mat');
initial_sick=helper.helper;

ratio = (frac_young * (n_h_term(i_young) + l_h_term(i_young))) / ...
        (frac_old *   (n_h_term(i_old)   + l_h_term(i_old)));
eps_old = initial_sick/(1 + ratio);
eps_young = initial_sick - eps_old;

M_fi(i_young, 1) = eps_young/frac_young;
M_fi(i_old, 1) = eps_old/frac_old;

for i_age = 1:n_age
    pi = (n_h_term(i_age) + l_h_term(i_age)) * Pi_star;
    
    M_h(i_age, 1) = (1 - pi) * (1 - M_fi(i_age, 1));
    M_fh(i_age, 1) = pi * (1 - M_fi(i_age, 1));
end

%-------------------------------------------------------------------------%
%                         Preparing for main loop                         %
%-------------------------------------------------------------------------%

% initial guess for general equilibrium variable
if flag_Pi_guess==1
    % does nothing (uses M_s and Pi stored in the memory)
else
    M_s = zeros(n_age, T);
    
    if flag_Pi==1
        Pi = S_Pi.Pi;
        I = S_Pi.I;
    else
        I = zeros(T, 1); % aggregate infectiousness at each point in time
        I(1:40) = linspace(0, 0.03, 40)';
        I(41:80) = linspace(0.03, 0, 40)';
        
        Pi = zeros(n_age, T);
        Pi(i_young, :) = 1 - exp(-Pi0 * I);
        Pi(i_old, :) = 1 - exp(-Pi0 * I);
    end
end

% if flag_fast = 1, the code doesn't evaluate cases 1 and 3 if case 2 gives
% variables inside bounds (in the young agent's problem)
flag_fast = 1;

%-------------------------------------------------------------------------%
%                    Distributions in the first period                    %
%-------------------------------------------------------------------------%

norm = 1;
iter = 0;

while norm > tol  || flag_fast==1
    iter = iter + 1;
    
    if flag_Pi==1
        M_s = S_Pi.M_s;
    end
    
    % if converged, makes one more iteration with slow code
    if norm <= tol
        flag_fast = 0;
    end
    
    %---------------------------------------%
    % Backward induction on value functions %
    %---------------------------------------%
    
    for t=T:-1:1
        
        % value in the next period
        if t==T
            for i_age = 1:n_age
                V_h_tomorrow(i_age) = V_h_term(i_age);
                V_i_tomorrow(i_age) = V_i_term(i_age);
                V_r_tomorrow(i_age) = V_h_term(i_age);
                V_f_tomorrow(i_age) = V_h_term(i_age);
                V_s_tomorrow(i_age) = V_s_term(i_age);
                
                V_h_private_tomorrow(i_age) = V_h_term(i_age);
                V_i_private_tomorrow(i_age) = V_i_term(i_age);
                V_r_private_tomorrow(i_age) = V_h_term(i_age);
                V_f_private_tomorrow(i_age) = V_h_term(i_age);
                V_s_private_tomorrow(i_age) = V_s_term(i_age);
            end
        else
            for i_age = 1:n_age
                V_h_tomorrow(i_age) = V_h(i_age, t+1);
                V_i_tomorrow(i_age) = V_i(i_age, t+1);
                V_r_tomorrow(i_age) = V_r(i_age, t+1);
                V_f_tomorrow(i_age) = V_f(i_age, t+1);
                V_s_tomorrow(i_age) = V_s(i_age, t+1);
                
                V_h_private_tomorrow(i_age) = V_h_private(i_age, t+1);
                V_i_private_tomorrow(i_age) = V_i_private(i_age, t+1);
                V_r_private_tomorrow(i_age) = V_r_private(i_age, t+1);
                V_f_private_tomorrow(i_age) = V_f_private(i_age, t+1);
                V_s_private_tomorrow(i_age) = V_s_private(i_age, t+1);
            end
        end
        
        % infectiousness in the previous period
        if t==1
            Pi_yesterday = [0, 0];
        else
            Pi_yesterday = [Pi(i_young, t-1), Pi(i_old, t-1)];
        end
        
        % value functions
        for i_age = 1:n_age
            
            % if we want young agents to use some variables of the old
            if flag_fake_young==1 && i_age==i_young
                i_age2 = i_old;
            else
                i_age2 = i_age;
            end
            
            % loop to solve problems of resistant, infected, healthy, and
            % fever
            for i_health_status = 1:4
                
                % sets some variables that depend on health status
                if i_health_status==1 % resistant
                    lambda_other = lambda_p_r(i_age, t);
                    lambda_other_private = 0;
                    L = 0;
                    
                    if i_age==i_young
                        n_guess = n_r_case3(i_age, t);
                        v_guess = v_r_case3(i_age, t);
                        x_guess = x_r_case3(i_age, t);
                        l_guess = l_r_case3(i_age, t);
                    else
                        x_guess = x_r(i_age, t);
                        l_guess = l_r(i_age, t);
                    end
                elseif i_health_status==2 % infected
                    lambda_other = lambda_i + lambda_p_i(i_age, t);
                    lambda_other_private = lambda_i;
                    L = 0;
                    
                    if i_age==i_young
                        n_guess = n_i_case3(i_age, t);
                        v_guess = v_i_case3(i_age, t);
                        x_guess = x_i_case3(i_age, t);
                        l_guess = l_i_case3(i_age, t);
                    else
                        x_guess = x_i(i_age, t);
                        l_guess = l_i(i_age, t);
                    end
                elseif i_health_status==3 % healthy
                    lambda_other = lambda_p_h(i_age, t);
                    lambda_other_private = 0;
                    L = beta(i_age2) * ...
                        (xi_p(i_age) * Pi(i_age, t) * (V_i_tomorrow(i_age) - V_h_tomorrow(i_age)) + ...
                        ((1 - xi_p(i_age)) * (Pi(i_age, t) + Pi_star) * (V_f_tomorrow(i_age) - V_h_tomorrow(i_age))));
                    
                    if i_age==i_young
                        n_guess = n_h_case3(i_age, t);
                        v_guess = v_h_case3(i_age, t);
                        x_guess = x_h_case3(i_age, t);
                        l_guess = l_h_case3(i_age, t);
                    else
                        x_guess = x_h(i_age, t);
                        l_guess = l_h(i_age, t);
                    end
                else % fever
                    belief = Pi_star / (Pi_yesterday(i_age) + Pi_star);
                    lambda_other = belief * 0 + ...
                        (1 - belief) * lambda_i + lambda_p_f(i_age, t);
                    lambda_other_private = belief * 0 + ...
                        (1 - belief) * lambda_i;
                    L = belief * beta(i_age2) * ...
                        (xi_p(i_age) * Pi(i_age, t) * (V_i_tomorrow(i_age) - V_h_tomorrow(i_age)) + ...
                        ((1 - xi_p(i_age)) * (Pi(i_age, t) + Pi_star) * (V_f_tomorrow(i_age) - V_h_tomorrow(i_age))));
                    
                    if i_age==i_young
                        n_guess = n_f_case3(i_age, t);
                        v_guess = v_f_case3(i_age, t);
                        x_guess = x_f_case3(i_age, t);
                        l_guess = l_f_case3(i_age, t);
                    else
                        x_guess = x_f(i_age, t);
                        l_guess = l_f(i_age, t);
                    end
                end
                
                % solves problem of agent (if not epidemiological model)
                if flag_epidemiological==1
                    if i_age==i_young
                        factor = factor_epidemiological;
                    else
                        factor = 1;
                    end
                    
                    c = c_h_term(i_age);
                    x = x_h_term(i_age);
                    n = n_h_term(i_age) * factor;
                    v = v_h_term(i_age);
                    l = l_h_term(i_age) * factor;
                    d = d_h_term(i_age);
                    c3 = 0; x3 = 0; n3 = 0; v3 = 0; l3 = 0; d3 = 0;
                else
                    if i_age==i_young
                        [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
                            solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
                            alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
                    else
                        [c, x, n, v, l, d] = ...
                            solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, ...
                            l_guess, x_guess);
                    end
                end
                
                % stores some variables and computes values
                if i_health_status==1 % resistant
                    
                    c_r(i_age, t) = c;
                    x_r(i_age, t) = x;
                    n_r(i_age, t) = n;
                    v_r(i_age, t) = v;
                    l_r(i_age, t) = l;
                    d_r(i_age, t) = d;
                    
                    c_r_case3(i_age, t) = c3;
                    x_r_case3(i_age, t) = x3;
                    n_r_case3(i_age, t) = n3;
                    v_r_case3(i_age, t) = v3;
                    l_r_case3(i_age, t) = l3;
                    d_r_case3(i_age, t) = d3;
                    
                    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
                    u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                    
                    V_r(i_age, t) = u + beta(i_age2) * V_r_tomorrow(i_age);
                    V_r_private(i_age, t) = u_private + beta(i_age2) * V_r_private_tomorrow(i_age);
                elseif i_health_status==2 % infected
                    c_i(i_age, t) = c;
                    x_i(i_age, t) = x;
                    n_i(i_age, t) = n;
                    v_i(i_age, t) = v;
                    l_i(i_age, t) = l;
                    d_i(i_age, t) = d;
                    
                    c_i_case3(i_age, t) = c3;
                    x_i_case3(i_age, t) = x3;
                    n_i_case3(i_age, t) = n3;
                    v_i_case3(i_age, t) = v3;
                    l_i_case3(i_age, t) = l3;
                    d_i_case3(i_age, t) = d3;
                    
                    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
                    u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                    
                    V_i(i_age, t) = u + ...
                        beta(i_age2) * phi(i_nosymptoms, i_age2) * V_r_tomorrow(i_age) + ...
                        beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_tomorrow(i_age) + ...
                        beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * (1-alpha(i_age2)) * V_i_tomorrow(i_age);

                    V_i_private(i_age, t) = u_private + ...
                        beta(i_age2) * phi(i_nosymptoms, i_age2) * V_r_private_tomorrow(i_age) + ...
                        beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_private_tomorrow(i_age) + ...
                        beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * (1-alpha(i_age2)) * V_i_private_tomorrow(i_age);
                elseif i_health_status==3 % healthy
                    c_h(i_age, t) = c;
                    x_h(i_age, t) = x;
                    n_h(i_age, t) = n;
                    v_h(i_age, t) = v;
                    l_h(i_age, t) = l;
                    d_h(i_age, t) = d;
                    
                    c_h_case3(i_age, t) = c3;
                    x_h_case3(i_age, t) = x3;
                    n_h_case3(i_age, t) = n3;
                    v_h_case3(i_age, t) = v3;
                    l_h_case3(i_age, t) = l3;
                    d_h_case3(i_age, t) = d3;
                    
                    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
                    u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                    
                    pi = (n + l) * Pi(i_age, t);
                    pif = pi + (n + l) * Pi_star;
                    
                    % V_h
                    CV1 = pi * V_i_tomorrow(i_age) + (1 - pi) * V_h_tomorrow(i_age);
                    CV2 = pif * V_f_tomorrow(i_age) + (1 - pif) * V_h_tomorrow(i_age);
                    
                    V_h(i_age, t) = u + ...
                        beta(i_age2) * xi_p(i_age) * CV1 + ...
                        beta(i_age2) * (1 - xi_p(i_age)) * CV2;
                    
                    % V_h_private
                    CV1 = pi * V_i_private_tomorrow(i_age) + (1 - pi) * V_h_private_tomorrow(i_age);
                    CV2 = pif * V_f_private_tomorrow(i_age) + (1 - pif) * V_h_private_tomorrow(i_age);
                    
                    V_h_private(i_age, t) = u_private + ...
                        beta(i_age2) * xi_p(i_age) * CV1 + ...
                        beta(i_age2) * (1 - xi_p(i_age)) * CV2;
                else % fever
                    c_f(i_age, t) = c;
                    x_f(i_age, t) = x;
                    n_f(i_age, t) = n;
                    v_f(i_age, t) = v;
                    l_f(i_age, t) = l;
                    d_f(i_age, t) = d;
                    
                    c_f_case3(i_age, t) = c3;
                    x_f_case3(i_age, t) = x3;
                    n_f_case3(i_age, t) = n3;
                    v_f_case3(i_age, t) = v3;
                    l_f_case3(i_age, t) = l3;
                    d_f_case3(i_age, t) = d3;
                    
                    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
                    u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                    
                    pi = (n + l) * Pi(i_age, t);
                    pif = pi + (n + l) * Pi_star;
                    
                    % V_f
                    CV1 = phi(i_nosymptoms, i_age2) * V_r_tomorrow(i_age) ...
                        + (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_tomorrow(i_age) ...
                        + (1 - phi(i_nosymptoms, i_age2)) * (1 - alpha(i_age2)) * V_i_tomorrow(i_age);

                    CV2 = xi_p(i_age) * (pi * V_i_tomorrow(i_age) + (1 - pi) * V_h_tomorrow(i_age)) ...
                        + (1 - xi_p(i_age)) * (pif * V_f_tomorrow(i_age) + (1 - pif) * V_h_tomorrow(i_age));
    
                    V_f(i_age, t) = u + ...
                        beta(i_age2) * (1 - belief) * CV1 + ...
                        beta(i_age2) * belief * CV2;

                    % V_f_private
                    CV1 = phi(i_nosymptoms, i_age2) * V_r_private_tomorrow(i_age) ...
                        + (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_private_tomorrow(i_age) ...
                        + (1 - phi(i_nosymptoms, i_age2)) * (1 - alpha(i_age2)) * V_i_private_tomorrow(i_age);

                    CV2 = xi_p(i_age) * (pi * V_i_private_tomorrow(i_age) + (1 - pi) * V_h_private_tomorrow(i_age)) ...
                        + (1 - xi_p(i_age)) * (pif * V_f_private_tomorrow(i_age) + (1 - pif) * V_h_private_tomorrow(i_age));

                    V_f_private(i_age, t) = u_private + ...
                        beta(i_age2) * (1 - belief) * CV1 + ...
                        beta(i_age2) * belief * CV2;

                    belief_vec(i_age, t) = belief;
                end
            end
            
            % symptoms
%             U = frac_young * M_s(i_young, t) * (1 - phi(i_symptoms, i_young)) + ... 
%                 frac_old * M_s(i_old, t) * (1 - phi(i_symptoms, i_old));
            
            U = frac_young * M_s(i_young, t)  + frac_old * M_s(i_old, t);
            
            delta_t = delta_1(i_age2) * min(Z/U, 1) + delta_2(i_age2) * max((U-Z)/U, 0);
            
            delta_vec_new(i_age, t) = delta_t;
            
            % smooth delta?
            if Z==1
                delta_vec(i_age, t) = delta_t;
            else
                weight_new = 0.5;
                delta_vec(i_age, t) = weight_new * delta_t + ...
                                   (1-weight_new) * delta_vec(i_age, t);
                delta_t = delta_vec(i_age, t);
            end
            
            V_s(i_age, t) = beta(i_age2) * (phi(i_symptoms, i_age2) * V_r_tomorrow(i_age) + ...
                (1 - phi(i_symptoms, i_age2)) * (1 - delta_t) * V_s_tomorrow(i_age));

            V_s_private(i_age, t) = beta(i_age2) * (phi(i_symptoms, i_age2) * V_r_private_tomorrow(i_age) + ...
                (1 - phi(i_symptoms, i_age2)) * (1-delta_t) * V_s_private_tomorrow(i_age));
        end
    end
    
    %------------------------------------%
    % Forward induction on distributions %
    %------------------------------------%
    
    M_s_other = M_s;
    
    for t = 1:T-1
        for i_age = 1:n_age
            % pis
            pi_h = (n_h(i_age, t) + l_h(i_age, t)) * Pi(i_age, t);
            pif_h = pi_h + (n_h(i_age, t) + l_h(i_age, t)) * Pi_star;

            pi_f = (n_f(i_age, t) + l_f(i_age, t)) * Pi(i_age, t);
            pif_f = pi_f + (n_f(i_age, t) + l_f(i_age, t)) * Pi_star;
            
            % stores probabilities in a vector to plot later
            pi_h_vec(i_age, t) = pi_h;
            pif_h_vec(i_age, t) = pif_h;
            pi_f_vec(i_age, t) = pi_f;
            pif_f_vec(i_age, t) = pif_f;
            
            % probability to die if develops symptoms
            U = frac_young * M_s_other(i_young, t) * (1 - phi(i_symptoms, i_young)) + ...
                frac_old * M_s_other(i_old, t) * (1 - phi(i_symptoms, i_old));
            
            delta_t = delta_1(i_age) * min(Z/U, 1) + delta_2(i_age) * max((U-Z)/U, 0);
            
            % masses
            M_h(i_age, t+1) = M_h(i_age, t) * (1 - xi_p(i_age)) * (1 - pif_h) ...
                + M_h(i_age, t) * xi_p(i_age) * (1 - pi_h) ...
                + M_fh(i_age, t) * xi_p(i_age) * (1 - pi_f) ... 
                + M_fh(i_age, t) * (1 - xi_p(i_age)) * (1 - pif_f);
            M_h(i_age, t+1) = M_h(i_age, t+1) * Delta(i_age);

            M_fh(i_age, t+1) = M_h(i_age, t) * (1 - xi_p(i_age)) * pif_h * (Pi_star/(Pi(i_age, t) + Pi_star)) ...
                + M_fh(i_age, t) * (1 - xi_p(i_age)) * pif_f * (Pi_star/(Pi(i_age, t) + Pi_star));
            M_fh(i_age, t+1) = M_fh(i_age, t+1) * Delta(i_age);

            M_fi(i_age, t+1) = M_h(i_age, t) * (1 - xi_p(i_age)) * pif_h * (Pi(i_age, t)/(Pi(i_age, t) + Pi_star)) ...
                + M_fh(i_age, t) * (1 - xi_p(i_age)) * pif_f * (Pi(i_age, t)/(Pi(i_age, t) + Pi_star));
            M_fi(i_age, t+1) = M_fi(i_age, t+1) * Delta(i_age);

            M_i(i_age, t+1) = M_h(i_age, t) * xi_p(i_age) * pi_h ...
                + M_fh(i_age, t) * xi_p(i_age) * pi_f ...
                + (M_fi(i_age, t) + M_i(i_age, t)) * (1 - phi(i_nosymptoms, i_age)) * (1-alpha(i_age));
            M_i(i_age, t+1) = M_i(i_age, t+1) * Delta(i_age);
            
            M_s(i_age, t+1) = (M_fi(i_age, t) + M_i(i_age, t)) * (1 - phi(i_nosymptoms, i_age)) * alpha(i_age) ...
                + M_s(i_age, t)*(1 - delta_t) * (1 - phi(i_symptoms, i_age));
            M_s(i_age, t+1) = M_s(i_age, t+1) * Delta(i_age);

            M_r(i_age, t+1) = (M_fi(i_age, t) + M_i(i_age, t)) * phi(i_nosymptoms, i_age) ...
                + M_s(i_age, t) * phi(i_symptoms, i_age) ...
                + M_r(i_age, t);
            M_r(i_age, t+1) = M_r(i_age, t+1) * Delta(i_age);
        end
    end
    
    % computes new Pi_hat
    if flag_Pi==1
        Pi_new = S_Pi.Pi;
        I_new = S_Pi.I;
    else
        Pi_hat(:, :) = 0;

        for t = 1:T
            if t >= t_vaccine

                Pi_hat(:, t) = 0;

            else
                % 1 - zeta term
                sum1 = 0;

                for i_age = 1:n_age
                    n_plus_l_f = n_f(i_age, t) + l_f(i_age, t);
                    n_plus_l_i = n_i(i_age, t) + l_i(i_age, t);

                    if i_age==i_young
                        frac = frac_young;
                    else
                        frac = frac_old;
                    end

                    sum1 = sum1 + frac * (n_plus_l_f * M_fi(i_age, t) ...
                                        + n_plus_l_i * M_i(i_age, t) ...
                                        + n_plus_l_s * M_s(i_age, t));
                end

                I_new(t) = sum1;

                for i_age = 1:n_age
                    % zeta term
                    n_plus_l_f = n_f(i_age, t) + l_f(i_age, t);
                    n_plus_l_i = n_i(i_age, t) + l_i(i_age, t);

                    if i_age==i_young
                        frac = frac_young;
                    else
                        frac = frac_old;
                    end

                    sum2 = frac * (n_plus_l_f * M_fi(i_age, t) ...
                                 + n_plus_l_i * M_i(i_age, t) ...
                                 + n_plus_l_s * M_s(i_age, t));

                    % computes new Pi_hat
                    Pi_hat(i_age, t) = (1-zeta) * Pi0 * sum1 + ...
                        (zeta * Pi0 / vartheta(i_age)) * sum2;
                end

            end
        end

        % new Pi
        Pi_new = 1 - exp(-Pi_hat);
    end
    
    % checks for convergence, updates norm
    X_new = [Pi_new(:); M_s(:); delta_vec_new(:)];
    X_old = [Pi(:); M_s_other(:); delta_vec(:)];
    
    norm = max(abs(X_new - X_old));
    
    if flag_Pi==1
        norm = 0;
    end
    
    % prints some convergence info
%     fprintf('norm: %.8e\n', norm)

    % updates Pi
    if Z < 1
        weight_new = 0.5;
    else
        weight_new = 1;
    end
    
    Pi = weight_new * Pi_new + (1 - weight_new) * Pi;
    I = I_new;
    
    % smooth M_s
    if Z < 1
        weight_new = 0.5;
        M_s = weight_new * M_s + (1 - weight_new) * M_s_other;
    end
end

fprintf('max(Pi(:, T)) = %.8e\n', max(Pi(:, T)))

% computes some distribution variables that were not needed in the main
% loop

% mass of agents with fever (cold or covid)
M_f = M_fh + M_fi;

% mass of infected agents (active cases)
M_i_all = M_i + M_fi + M_s;

for t = 1:T-1
    for i_age = 1:n_age
        % delta_t
        U = frac_young * M_s_other(i_young, t) * (1 - phi(i_symptoms, i_young)) + ...
            frac_old * M_s_other(i_old, t) * (1 - phi(i_symptoms, i_old));
        
        delta_t = delta_1(i_age) * min(Z/U, 1) + delta_2(i_age) * max((U-Z)/U, 0);
        
        % number of tests per period
        M_t(i_age, t+1) = M_h(i_age, t) * Delta(i_age) * pif_h_vec(i_age, t) * xi_p(i_age) + ...
                         M_fh(i_age, t) * Delta(i_age) * pif_f_vec(i_age, t) * xi_p(i_age);
        
        % cumulative mass of recovered agents (doesn't take into account
        % natural deaths)
        M_rn(i_age, t+1) = (M_fi(i_age, t) + M_i(i_age, t)) * phi(i_nosymptoms, i_age) ...
            + M_s(i_age, t) * phi(i_symptoms, i_age);
        M_rn(i_age, t+1) = M_rn(i_age, t+1) * Delta(i_age);
        M_rn(i_age, t+1) = M_rn(i_age, t+1) + M_rn(i_age, t);
        
        % mass of dead agents (only covid deaths)
        M_d(i_age, t+1) = M_d(i_age, t) ...
            + M_s(i_age, t) * Delta(i_age) * delta_t * (1 - phi(i_symptoms, i_age));
        
        % mass of deceased agents (only natural deaths)     
        M_dn(i_age, t+1) = M_dn(i_age, t) ...
            + M_h(i_age, t) * (1 - Delta(i_age)) ...
            + M_f(i_age, t) * (1 - Delta(i_age)) ...
            + M_i(i_age, t) * (1 - Delta(i_age)) ...
            + M_r(i_age, t) * (1 - Delta(i_age)) ...
            + M_s(i_age, t) * (1 - Delta(i_age));
        
        % mass of covid cases (flux)
        N_c_f(i_age, t+1) = M_h(i_age, t) * pi_h_vec(i_age, t) + M_fh(i_age, t) * pi_f_vec(i_age, t);
        N_c_f(i_age, t+1) = N_c_f(i_age, t+1) * Delta(i_age); % some die before catching COVID

        % mass of covid cases (stock)
        N_c_s(i_age, t+1) = N_c_f(i_age, t+1) + N_c_s(i_age, t);
    end
end

% checks if distributions sum up to one (ok)
% M_h + M_fh + M_fi + M_i + M_s + M_r + M_d + M_dn

% gdp
tau_h = alpha1 * v_h - alpha2 * v_h.^2;
tau_f = alpha1 * v_f - alpha2 * v_f.^2;
tau_i = alpha1 * v_i - alpha2 * v_i.^2;
tau_r = alpha1 * v_r - alpha2 * v_r.^2;

gdp = w * ((n_h(i_young, :) + tau_h(i_young, :)) .* M_h(i_young,:) ...
         + (n_f(i_young, :) + tau_f(i_young, :)) .* M_f(i_young,:) ...
         + (n_i(i_young, :) + tau_i(i_young, :)) .* M_i(i_young,:) ...
         + (n_r(i_young, :) + tau_r(i_young, :)) .* M_r(i_young,:));

% gdp per capita
gdp_pc = gdp ./ (1 - M_d(i_young, :) - M_dn(i_young, :));

% time to peak
[~, t_peak_I] = max(I);
  
% calculates average welfare in first period with COVID
E_W = 0;

for i_age = 1:n_age
    E_W = E_W + frac_age(i_age) * ...
        (M_h(i_age, 1) * V_h_private(i_age, 1) + ...
        (M_fh(i_age, 1) + M_fi(i_age, 1)) * V_f_private(i_age, 1));
end



















