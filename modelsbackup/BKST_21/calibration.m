

%-------------------------------------------------------------------------%
%                         Steady state calibration                        %
%-------------------------------------------------------------------------%

% some targets
n_relative_data = 1 - 0.36; % n falls by 36% in counterfactual
income_relative_data = 1 - 0.1; % income falls by 10% in counterfactual

% prepares to calibrate
guess = ones(6, 1);
fun = @(in) calib_ND(in, n_data, v_data, l_data, x_over_income_data, n_relative_data, income_relative_data, w, rho);

% use fminsearch
solution = fminsearch(fun, guess);

% gets calibrated parameters
[~, theta, gamma, lambda_d, alpha1, alpha2, L2, x_over_income, n, v, l, n_relative, income_relative] = ...
    calib_ND(solution, n_data, v_data, l_data, x_over_income_data, n_relative_data, income_relative_data, w, rho);

% shows a table with calibrated parameter variables
K = zeros(6, 1);
K(:, 1) = [theta, gamma, lambda_d, alpha1, alpha2, L2];
K = array2table(K);
K.Properties.VariableNames = {'Parameter'};
K.Properties.RowNames = {'theta', 'gamma', 'lambda_d', 'alpha1', 'alpha2', 'L'};
disp(K)

% shows a table with calibration results
K = zeros(6, 2);
K(:, 1) = [x_over_income_data; n_data; v_data; l_data; n_relative_data; income_relative_data];
K(:, 2) = [x_over_income; n; v; l; n_relative; income_relative];
K = array2table(K);
K.Properties.VariableNames = {'Data', ' Model'};
K.Properties.RowNames = {'x/income', ' n', 'v', 'l', 'n_relative', 'income_relative'};
disp(K)

% variables in the s.s.
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

K = zeros(6, 2);
K(:, 1) = [c, x, n, v, l, d];

% variables in the counterfactual
lambda_other = 0;
L = L2;
n_guess = 0;
v_guess = 0;
x_guess = 0;
l_guess = 0;
flag_fast = 0;

[c, x, n, v, l, d] = ...
    solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);

K(:, 2) = [c, x, n, v, l, d];

% plot table
K = array2table(K);
K.Properties.VariableNames = {'no_COVID', ' counterfactual'};
K.Properties.RowNames = {'c', ' x', 'n', 'v', 'l', 'd'};
disp(K)

%-------------------------------------------------------------------------%
%                                lambda_i                                 %
%-------------------------------------------------------------------------%

% finds lambda_other such that time spent at home increases by x%
% time_home_relative_data = 1.5; % 50%
% time_home_relative_data = 1.25; % 25%
% time_home_relative_data = 1.26; % 26%
% time_home_relative_data = 1.75; % 75%
time_home_relative_data = 1.90; % 90%

% first, gets d in steady state
lambda_other = 0;
L = 0;
n_guess = 0;
v_guess = 0;
x_guess = 0;
l_guess = 0;
flag_fast = 0;

[~, ~, ~, v_ND, ~, d_ND] = ...
    solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);

time_home_ND = v_ND + d_ND

% gets variables
guess = 1;
fun = @(in) calib_lambda_other(in, w, lambda_d, gamma, theta, rho, ...
    alpha1, alpha2, time_home_relative_data, time_home_ND);

solution = fsolve(fun, guess);

[out, lambda_other, time_home] = calib_lambda_other(solution, w, ...
    lambda_d, gamma, theta, rho, alpha1, alpha2, time_home_relative_data, time_home_ND)

disp('time_home/time_home_ND')
disp(time_home/time_home_ND)

%-------------------------------------------------------------------------%
%                             Local functions                             %
%-------------------------------------------------------------------------%

% calibrates parameters in the no-disease (ND) world
function [out, theta, gamma, lambda_d, alpha1, alpha2, L, x_over_income, n, v, l, n_relative, income_relative] = ...
    calib_ND(in, n_data, v_data, l_data, x_over_income_data, n_relative_data, income_relative_data, w, rho)
    
    % opens input vector
    theta = in(1);
    gamma = in(2);
    lambda_d = in(3);
    alpha1 = in(4);
    alpha2 = in(5);
    L2 = in(6);

    % converts to desired range (theta)
    lb = 0;
    ub = 1;
    theta = (lb + ub)/2 + sin(theta) * (ub - lb)/2;
    
    % converts to desired range (others)
    gamma = gamma^2;
    lambda_d = lambda_d^2;
    alpha1 = alpha1^2;
    alpha2 = alpha2^2;
    L2 = -L2^2;
    
    % solves problem in no-COVID world (simulation 1)
    lambda_other = 0;
    L = 0;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;
    
    [~, x1, n1, v1, l1, ~] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    
    % some variables in simulation 1
    tau1 = alpha1 * v1 - alpha2 * v1^2;
    income1 = w * (n1 + tau1);
    x_over_income = x1/income1;
    
    % solves problem with a given L (simulation 2)
    lambda_other = 0;
    L = L2;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;
    
    [~, ~, n2, v2, ~, ~] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    
    % some variables in simulation 2
    tau2 = alpha1 * v2 - alpha2 * v2^2;
    income2 = w * (n2 + tau2);
    n_relative = n2/n1;
    income_relative = income2/income1;

    % defines output variable
    n = n1;
    v = v1;
    l = l1;
    
    out = zeros(6, 1);
    out(1) = x_over_income/x_over_income_data - 1;
    out(2) = n/n_data - 1;
    out(3) = v/v_data - 1;
    out(4) = l/l_data - 1;
    out(5) = n_relative/n_relative_data - 1;
    out(6) = income_relative/income_relative_data - 1;
    
    weights = ones(6, 1);
    weights(5) = 3;
    weights(6) = 7;
    
    out = sum(weights .* out.^2);
end

% calibrates lambda_other
function [out, lambda_other, time_home] = ...
    calib_lambda_other(in, w, lambda_d, gamma, theta, rho, ...
    alpha1, alpha2, time_home_relative_data, time_home_ND)

    % defines and transforms inputs
    lambda_other = in(1)^2;
    
    % solves problem of agent
    L = 0;
    n_guess = 0;
    v_guess = 0;
    x_guess = 0;
    l_guess = 0;
    flag_fast = 0;
    
    [~, ~, ~, v, ~, d] = ...
        solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    
    % defines output variables
    time_home = d + v;
    time_home_relative = time_home/time_home_ND;
    
    out = time_home_relative/time_home_relative_data - 1;
end

















