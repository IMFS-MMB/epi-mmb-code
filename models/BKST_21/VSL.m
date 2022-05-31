
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

% choose age group
i_age = i_young;
% i_age = i_old;

if i_age==i_young
    flag_young = 1;
    delta_bm = 1 - survival_young;
else
    flag_young = 0;
    delta_bm = 1 - survival_old;
end

%-------------------------------------------------------------------------%
%                    1: calibrated b to fit a given VSL                   %
%-------------------------------------------------------------------------%

% counterfactual situation: young dies
N = 13; % if quarterly (n. of weeks in a quarter)
% N = 52; % if yearly (n. of weeks in a year)

% increase in prob. to die per quarter
p = 1/10000;

% new weekly prob. of death
delta_new = 1 - ((1 - delta_bm)^N - p)^(1/N);

% some numbers in USD
income_quarter_USD_data = 15000; % if quarterly
% income_quarter_USD_data = 60000; % if yearly
% willingness_to_accept_USD_data = 1000;
willingness_to_accept_USD_data = 930;

%---%
% 1 %
%---%

disp('---------------------')
disp('          1          ')
disp('---------------------')

% finds w and b such that consumer is indifferent between counterfactual
% and no-COVID, and such that new quarterly income is 15,000 + 930
myfun = @(in) fun(in, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_bm, delta_new, income_quarter_USD_data, ...
    willingness_to_accept_USD_data, N, time_discount);

guess = [1, 3];

% solution = fsolve(myfun, guess);
solution = fminsearch(myfun, guess);

[~, w_new, b_new, V_ND, income_quarter_USD, V] = ...
    fun(solution, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_bm, delta_new, income_quarter_USD_data, ...
    willingness_to_accept_USD_data, N, time_discount)

%-------------------------------------------------------------------------%
%       2: for a given b, how much the agent is willing to accept?        %
%-------------------------------------------------------------------------%

disp('---------------------')
disp('          2          ')
disp('---------------------')

% b = 11.0049;
% b = 4.63;
% b = 4.45;
% b = 11.5469;
b = 5.6;
% b = 11.0121;
w = 1;

% no-COVID world
lambda_other = 0;
L = 0;
n_guess = 0;
v_guess = 0;
x_guess = 0;
l_guess = 0;
flag_fast = 0;

if flag_young==1
    [c, x, n, v, l, d] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
else
    [c, x, n, v, l, d] = ...
        solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, l_guess, x_guess);
end

% income
if flag_young==1
    income_ND = w * (n + alpha1 * v - alpha2 * v^2);
else
    income_ND = w_bar;
end

income_quarter_ND = income_ND * N

% value
u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
V_ND = u/(1 - time_discount * (1 - delta_bm));

% finds new w that makes agent indifferent btw. no-COVID and positive delta
myfun = @(in) fun2(in, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_new, b, V_ND, N, w_bar, flag_young, time_discount);

guess = 1;

solution = fzero(myfun, guess);

[out, w_new, income_quarter] = ...
    fun2(solution, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_new, b, V_ND, N, w_bar, flag_young, time_discount)

income_quarter_USD = income_quarter * income_quarter_USD_data / income_quarter_ND

disp('income_quarter_USD - income_quarter_USD_data =')
disp(income_quarter_USD - income_quarter_USD_data)

%-------------------------------------------------------------------------%
%                             Local functions                             %
%-------------------------------------------------------------------------%

function [out, w_new, b, V_ND, income_quarter_USD, V] = ...
    fun(in, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_bm, delta_new, income_quarter_USD_data, ...
    willingness_to_accept_USD_data, N, time_discount)
    
    % defines input variables
    w_new = in(1)^2;
    b = in(2)^2;
    
    % no-COVID world
    w = 1;
    lambda_other = 0;
    L = 0;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;

    [c, x, n, v, l, d] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);

    % income
    income_ND = w * (n + alpha1 * v - alpha2 * v^2);
    income_quarter_ND = income_ND * N;

    % value
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    V_ND = u/(1 - time_discount * (1 - delta_bm));
    
    % counterfactual (with death probability)
    w = w_new;
    lambda_other = 0;
    L = 0;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;

    [c, x, n, v, l, d] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    
    % income
    income = w * (n + alpha1 * v - alpha2 * v^2);
    income_quarter = income * N;
    income_quarter_USD = income_quarter * income_quarter_USD_data / income_quarter_ND;
    
    % value
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    V = u/(1 - time_discount * (1 - delta_new));
    
    % output variable
    out = zeros(2, 1);
    out(1) = income_quarter_USD/(income_quarter_USD_data + willingness_to_accept_USD_data) - 1;
    out(2) = V/V_ND - 1;
    
    % fminsearch
    out = sum(out.^2);
end


function [out, w_new, income_quarter] = ...
    fun2(in, lambda_d, gamma, theta, rho, alpha1, alpha2, ...
    delta_new, b, V_ND, N, w_bar, flag_young, time_discount)
    
    % defines input variables
    w_new = in(1)^2;
    
    % counterfactual (with death probability)
    if flag_young==1
        w = w_new;
    else
        w_bar = w_new;
    end
    
    lambda_other = 0;
    L = 0;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;

    if flag_young==1
        [c, x, n, v, l, d] = ...
            solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
            alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    else
        [c, x, n, v, l, d] = ...
            solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, l_guess, x_guess);
    end
    
    % value
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    V = u/(1 - time_discount * (1 - delta_new));
    
    % income
    if flag_young==1
        income = w * (n + alpha1 * v - alpha2 * v^2);
    else
        income = w_bar;
    end
    
    income_quarter = income * N;
    
    % output variable
    out = V/V_ND - 1;
end













