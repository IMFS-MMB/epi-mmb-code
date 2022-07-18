

% distribution of agents across age groups
frac_old = 0.214;
% frac_old = 0.16;
frac_young = 1 - frac_old;
frac_age = [frac_young, frac_old];

% % Ferguson
% phi(i_nosymptoms, i_young) = 0.991289141;
% phi(i_nosymptoms, i_old) = 0.893010253;
% phi(i_symptoms, i_young) = 0.284090909;
% phi(i_symptoms, i_old) = 0.284090909;
% delta(i_young) = 0.372680985;
% delta(i_old) = 0.371014263;

% case 1: with teleworking
Pi0 = 13.425; % baseline parameters
% Pi0 = 13.7; % Ferguson
n = [0.3284, 0];
l = [0.154464286, 0.205531943];
l_bar = [0.158, 0.158];

% % case 2: without teleworking
% Pi0 = 12.127;
% n = [0.357, 0];
% l = [0.154464286, 0.205531943];
% l_bar = [0.158, 0.158];

%-------------------------------------------------------------------------%
%                              Calculates R0                              %
%-------------------------------------------------------------------------%

avg_time_out = 0;

for i_age = 1:n_age
    avg_time_out = avg_time_out + ...
        frac_age(i_age) * (n(i_age) + l(i_age));
end

R0_vec = zeros(n_age, 1);

for i_age = 1:n_age
    T_i = 1/(phi(i_nosymptoms, i_age) + (1 - phi(i_nosymptoms, i_age)) * alpha(i_age));
    T_s = 1/(phi(i_symptoms, i_age) + (1 - phi(i_symptoms, i_age)) * delta_1(i_age));
    prob_s = (1 - phi(i_nosymptoms, i_age)) * alpha(i_age) / (1 - (1 - phi(i_nosymptoms, i_age)) * (1 - alpha(i_age)));
    R0_vec(i_age) = Pi0 * avg_time_out * (T_i * (n(i_age) + l(i_age)) + prob_s * T_s * l_bar(i_age));
end

R0 = 0;

for i_age = 1:n_age
    R0 = R0 + frac_age(i_age) * R0_vec(i_age);
end

fprintf('R0 = %.9f\n', R0)

%-------------------------------------------------------------------------%
%                           Calculates Pi_star                            %
%-------------------------------------------------------------------------%

n_colds_per_year = 3;
n_weeks_in_year = 52;
weekly_rate = 1 / (n_weeks_in_year / n_colds_per_year);

Pi_star_hat = weekly_rate / (n(i_young) + l(i_young));
Pi_star = 1 - exp(-Pi_star_hat);

fprintf('Pi_star = %.9f\n', Pi_star)











