
% variable to count number of iterations
planner_iter = 0;

% starts file to store info for each iteration
filename = 'planner.csv';

dlmwrite(filename, 0, 'delimiter', ',');

% policy only changes in each 'n_weeks_planner' weeks
n_weeks_planner = 4;
% n_weeks_planner = 1;
length_guess = floor(min(t_vaccine, T) / n_weeks_planner);

% loads initial guess
% load('planner_guess')

% is the lockdown age-specific?
n_age_lambda_p = n_age; % yes
% n_age_lambda_p = 1; % no (universal)

% initial guess
% guess = ones(n_age_lambda_p, length_guess);
% guess = sqrt(0.5 * ones(n_age_lambda_p, length_guess));
% guess = sqrt(lambda_p);
% guess = zeros(n_age_lambda_p, length_guess);
lambda_p = [1.8185434, 0.721996808, 0.787572037, 0.772584232, 0.765966776, 0.760328664, 0.756039114, 0.752209708, 0.748944961, 0.745696235, 0.742649625, 0.739778431, 0.737201439, 0.734600472, 0.73208844, 0.729548129, 0.727391598, 0.725000243, 0.722839001, 0.720694865, 0.71858115, 0.71659909, 0.714529672, 0.712540563, 0.710659096, 0.70874489, 0.706861036, 0.705062538, 0.703228611, 0.701433716, 0.699477694, 0.697857026, 0.695873246, 0.69420826, 0.692325899, 0.690583692, 0.688708425, 0.686718328, 0.684999398, 0.68309028, 0.681053797, 0.679043621, 0.677043247, 0.674915661, 0.672787718, 0.670557056, 0.668348371, 0.66590882, 0.663387502, 0.660832753, 0.658085677, 0.655221465, 0.652138508, 0.648971507, 0.64547989, 0.641981221, 0.638101932, 0.633799949, 0.629290811, 0.624340101, 0.618853333, 0.612903519, 0.606240069, 0.598786802, 0.590388998, 0.580740118, 0.569594804, 0.556662515, 0.541214787, 0.522486452, 0.499558587, 0.470798267, 0.433595147, 0.384338265, 0.316832029, 0.220263343, 0.075109494, 0.000462144];
lambda_p(2, :) = [0.909976433, 0.725530255, 0.543798508, 0.560969444, 0.552320156, 0.54805362, 0.544507129, 0.541581659, 0.538800718, 0.536255559, 0.53402191, 0.531687497, 0.52956849, 0.527419386, 0.525484733, 0.523460598, 0.521737908, 0.519949821, 0.518059929, 0.516538761, 0.514770142, 0.513124641, 0.51162695, 0.510132948, 0.508275225, 0.507037358, 0.505417864, 0.504021812, 0.502547181, 0.501008453, 0.499618856, 0.498320909, 0.496707433, 0.49530895, 0.49389689, 0.492606866, 0.490895703, 0.489574252, 0.48805453, 0.486497455, 0.484884904, 0.483353416, 0.481674705, 0.480068284, 0.478366163, 0.476686785, 0.474737876, 0.473029785, 0.470850596, 0.469038642, 0.466542684, 0.464473261, 0.461906946, 0.459358242, 0.45676878, 0.453950646, 0.450839167, 0.44746343, 0.443861746, 0.440012504, 0.435439529, 0.430864917, 0.425474734, 0.419772866, 0.413023502, 0.405340418, 0.396739277, 0.386482707, 0.374560802, 0.359620089, 0.341624153, 0.319378491, 0.290511519, 0.253165333, 0.203517662, 0.136420351, 0.05041662, 0.002132444];
guess = sqrt(lambda_p);

% if we want testing and quarantine
% xi_p = 0.5 * ones(n_age, 1); % 50% testing
% new_lambda = 2.45096990297207; % 75% lockdown for those who test positive
% new_lambda = 32.3598330005704; % 90% lockdown for those who test positive

% with selective mixing
zeta = zeta_data;

% if no testing and quarantine
xi_p = zeros(n_age, 1);

% solves planner's problem
opts = optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6);
solution = fminsearch(@fun, guess, opts);

%-------------------------------------------------------------------------%
%                              Local function                             %
%-------------------------------------------------------------------------%

% function that takes some lambda_p, calculates general equilibrium and 
% gets welfare
function out = fun(in)
    
    % gets some variables from base workspace
    t_vaccine = evalin('base', 't_vaccine');
    T = evalin('base', 'T');
    n_age = evalin('base', 'n_age');
    n_weeks_planner = evalin('base', 'n_weeks_planner');
    filename = evalin('base', 'filename');
    xi_p = evalin('base', 'xi_p');
    
    if xi_p(1) > 0
        lambda_i = evalin('base', 'lambda_i');
        new_lambda = evalin('base', 'new_lambda');
    end

    % updates planner_iter
    evalin('base', 'planner_iter = planner_iter + 1;')
    
    planner_iter = evalin('base', 'planner_iter');
    fprintf('planner_iter: %i\n', planner_iter)

    % defines lambda
    in = in.^2; % transforms input variable to positive
    
    if size(in, 1)==1
        in = repmat(in, [n_age, 1]);
    end
    
    lambda_p = zeros(n_age, min(t_vaccine, T));
    
    for i_age = 1:n_age
        i0 = 1;
        i1 = n_weeks_planner;
        
        for t = 1:length(in)
            if t < length(in)
                lambda_p(i_age, i0:i1) = in(i_age, t);
            else
                lambda_p(i_age, i0:end) = in(i_age, t);
            end

            i0 = i1 + 1;
            i1 = i0 + n_weeks_planner - 1;
        end
    end
    
    % defines lambda_p vectors
    lambda_p_i = zeros(n_age, T);
    lambda_p_h = zeros(n_age, T);
    lambda_p_r = zeros(n_age, T);
    lambda_p_f = zeros(n_age, T);
    
    for i_age = 1:n_age
        lambda_p_h(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_r(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        lambda_p_f(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
    end
    
    if xi_p(1) > 0
        for i_age = 1:n_age
            lambda_p_i(i_age, 1:length(lambda_p)) = max(new_lambda - lambda_i, 0);
        end
    else
        for i_age = 1:n_age
            lambda_p_i(i_age, 1:length(lambda_p)) = lambda_p(i_age, :);
        end
    end
    
    % sends lambda_p variables to base workspace
    assignin('base', 'lambda_p_i', lambda_p_i)
    assignin('base', 'lambda_p_h', lambda_p_h)
    assignin('base', 'lambda_p_r', lambda_p_r)
    assignin('base', 'lambda_p_f', lambda_p_f)
    
    % calculates equilibrium
    evalin('base', 'equilibrium')
    
    % gets some variables from base
    E_W = evalin('base', 'E_W');
    V_h_private = evalin('base', 'V_h_private');
    
    % stores information in csv file
    for i_age = 1:n_age
        row = [planner_iter, E_W, V_h_private(i_age, 1), lambda_p(i_age, :)];
        dlmwrite(filename, row, 'delimiter', ',', 'precision', 9, '-append');
    end
    
    % output variable
    out = -E_W;
    
    % uses Pi and M_s calculated in this iteration as initial guesses for
    % next iteration
%     evalin('base', 'flag_Pi_guess = 1;')
end





