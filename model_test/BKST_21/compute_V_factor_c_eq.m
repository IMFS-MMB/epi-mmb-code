
% this code receives a factor_c_eq greater than zero (that multiplies
% consumption in all periods) and returns the value of agents for this new
% consumption in all periods

% pre-allocation
V_h_term_new = zeros(n_age, 1);
V_i_term_new = zeros(n_age, 1);
V_s_term_new = zeros(n_age, 1);
V_h_private_new = zeros(n_age, T);
V_f_private_new = zeros(n_age, T);
V_i_private_new = zeros(n_age, T);
V_r_private_new = zeros(n_age, T);
V_s_private_new = zeros(n_age, T);

%-------------------------------------------------------------------------%
%                   Variables in the terminal situation                   %
%-------------------------------------------------------------------------%

% healthy
for i_age = 1:n_age
    
    % new choices
    c = S_c_eq.c_h_term(i_age) * factor_c_eq;
    x = S_c_eq.x_h_term(i_age);
    n = S_c_eq.n_h_term(i_age);
    v = S_c_eq.v_h_term(i_age);
    l = S_c_eq.l_h_term(i_age);
    d = S_c_eq.d_h_term(i_age);
    
    % computes value
    lambda_other = 0;
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    
    if S_c_eq.flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    V_h_term_new(i_age) = u/(1-beta(i_age2));
    
end

% symptoms
for i_age = 1:n_age
    if S_c_eq.flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    % there are free hospital beds in terminal situation
    delta_t = delta_1(i_age2);
    
    V_s_term_new(i_age) = beta(i_age2) * phi(i_symptoms, i_age2) * V_h_term_new(i_age) / ...
        (1 - beta(i_age2) * (1 - phi(i_symptoms, i_age2)) * (1 - delta_t));
end

% infected
for i_age = 1:n_age
    
    % new choices
    c = S_c_eq.c_i_term(i_age) * factor_c_eq;
    x = S_c_eq.x_i_term(i_age);
    n = S_c_eq.n_i_term(i_age);
    v = S_c_eq.v_i_term(i_age);
    l = S_c_eq.l_i_term(i_age);
    d = S_c_eq.d_i_term(i_age);
    
    % computes value
    lambda_other = lambda_i;
    u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta);
    
    if S_c_eq.flag_fake_young==1 && i_age==i_young
        i_age2 = i_old;
    else
        i_age2 = i_age;
    end
    
    V_i_term_new(i_age) = u + beta(i_age2) * phi(i_nosymptoms, i_age2) * V_h_term_new(i_age) + ...
                       beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_term_new(i_age);
                            
    V_i_term_new(i_age) = V_i_term_new(i_age) / ...
        (1 - beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * (1 - alpha(i_age2)));
end

%-------------------------------------------------------------------------%
%                     Variables during the transition                     %
%-------------------------------------------------------------------------%

for t=T:-1:1
        
    % value in the next period
    if t==T
        for i_age = 1:n_age
            V_h_private_tomorrow(i_age) = V_h_term_new(i_age);
            V_i_private_tomorrow(i_age) = V_i_term_new(i_age);
            V_r_private_tomorrow(i_age) = V_h_term_new(i_age);
            V_f_private_tomorrow(i_age) = V_h_term_new(i_age);
            V_s_private_tomorrow(i_age) = V_s_term_new(i_age);
        end
    else
        for i_age = 1:n_age
            V_h_private_tomorrow(i_age) = V_h_private_new(i_age, t+1);
            V_i_private_tomorrow(i_age) = V_i_private_new(i_age, t+1);
            V_r_private_tomorrow(i_age) = V_r_private_new(i_age, t+1);
            V_f_private_tomorrow(i_age) = V_f_private_new(i_age, t+1);
            V_s_private_tomorrow(i_age) = V_s_private_new(i_age, t+1);
        end
    end

    % infectiousness in the previous period
    if t==1
        Pi_yesterday = [0, 0];
    else
        Pi_yesterday = [S_c_eq.Pi(i_young, t-1), S_c_eq.Pi(i_old, t-1)];
    end

    % value functions
    for i_age = 1:n_age

        % if we want young agents to use some variables of the old
        if S_c_eq.flag_fake_young==1 && i_age==i_young
            i_age2 = i_old;
        else
            i_age2 = i_age;
        end

        % loop to compute new values of resistant, infected, healthy, and
        % fever
        for i_health_status = 1:4

            % sets some variables that depend on health status
            if i_health_status==1 % resistant
                lambda_other_private = 0;
                c = S_c_eq.c_r(i_age, t) * factor_c_eq;
                x = S_c_eq.x_r(i_age, t);
                n = S_c_eq.n_r(i_age, t);
                v = S_c_eq.v_r(i_age, t);
                l = S_c_eq.l_r(i_age, t);
                d = S_c_eq.d_r(i_age, t);
            elseif i_health_status==2 % infected
                lambda_other_private = lambda_i;
                c = S_c_eq.c_i(i_age, t) * factor_c_eq;
                x = S_c_eq.x_i(i_age, t);
                n = S_c_eq.n_i(i_age, t);
                v = S_c_eq.v_i(i_age, t);
                l = S_c_eq.l_i(i_age, t);
                d = S_c_eq.d_i(i_age, t);
            elseif i_health_status==3 % healthy
                lambda_other_private = 0;
                c = S_c_eq.c_h(i_age, t) * factor_c_eq;
                x = S_c_eq.x_h(i_age, t);
                n = S_c_eq.n_h(i_age, t);
                v = S_c_eq.v_h(i_age, t);
                l = S_c_eq.l_h(i_age, t);
                d = S_c_eq.d_h(i_age, t);
            else % fever
                belief = Pi_star / (Pi_yesterday(i_age) + Pi_star);
                lambda_other_private = belief * 0 + ...
                    (1 - belief) * lambda_i;
                c = S_c_eq.c_f(i_age, t) * factor_c_eq;
                x = S_c_eq.x_f(i_age, t);
                n = S_c_eq.n_f(i_age, t);
                v = S_c_eq.v_f(i_age, t);
                l = S_c_eq.l_f(i_age, t);
                d = S_c_eq.d_f(i_age, t);
            end

            % stores some variables and computes values
            if i_health_status==1 % resistant
                u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                V_r_private_new(i_age, t) = u_private + beta(i_age2) * V_r_private_tomorrow(i_age);
            elseif i_health_status==2 % infected
                u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);
                V_i_private_new(i_age, t) = u_private + ...
                    beta(i_age2) * phi(i_nosymptoms, i_age2) * V_r_private_tomorrow(i_age) + ...
                    beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_private_tomorrow(i_age) + ...
                    beta(i_age2) * (1 - phi(i_nosymptoms, i_age2)) * (1-alpha(i_age2)) * V_i_private_tomorrow(i_age);
            elseif i_health_status==3 % healthy
                u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);

                pi = (n + l) * S_c_eq.Pi(i_age, t);
                pif = pi + (n + l) * Pi_star;

                % V_h_private
                CV1 = pi * V_i_private_tomorrow(i_age) + (1 - pi) * V_h_private_tomorrow(i_age);
                CV2 = pif * V_f_private_tomorrow(i_age) + (1 - pif) * V_h_private_tomorrow(i_age);

                V_h_private_new(i_age, t) = u_private + ...
                    beta(i_age2) * S_c_eq.xi_p(i_age) * CV1 + ...
                    beta(i_age2) * (1 - S_c_eq.xi_p(i_age)) * CV2;
            else % fever
                u_private = utility(c, x, l, d, v, lambda_d, lambda_other_private, b, gamma, rho, theta);

                pi = (n + l) * S_c_eq.Pi(i_age, t);
                pif = pi + (n + l) * Pi_star;

                % V_f_private
                CV1 = phi(i_nosymptoms, i_age2) * V_r_private_tomorrow(i_age) ...
                    + (1 - phi(i_nosymptoms, i_age2)) * alpha(i_age2) * V_s_private_tomorrow(i_age) ...
                    + (1 - phi(i_nosymptoms, i_age2)) * (1 - alpha(i_age2)) * V_i_private_tomorrow(i_age);

                CV2 = S_c_eq.xi_p(i_age) * (pi * V_i_private_tomorrow(i_age) + (1 - pi) * V_h_private_tomorrow(i_age)) ...
                    + (1 - S_c_eq.xi_p(i_age)) * (pif * V_f_private_tomorrow(i_age) + (1 - pif) * V_h_private_tomorrow(i_age));

                V_f_private_new(i_age, t) = u_private + ...
                    beta(i_age2) * (1 - belief) * CV1 + ...
                    beta(i_age2) * belief * CV2;
            end
        end

        % symptoms
        delta_t = S_c_eq.delta_vec(i_age, t);
        V_s_private_new(i_age, t) = beta(i_age2) * (phi(i_symptoms, i_age2) * V_r_private_tomorrow(i_age) + ...
            (1 - phi(i_symptoms, i_age2)) * (1-delta_t) * V_s_private_tomorrow(i_age));
    end
end








