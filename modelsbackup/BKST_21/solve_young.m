function [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
    solve_young(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast)

%     if L > 1e-6
%         error('solve_young: L>0')
%     end
    
    % goes over several possibilities
    if alpha1==0
        % only case 1
        [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
            cases_1(w, lambda_d, lambda_other, L, gamma, theta, rho);
    elseif alpha1 > 0 && alpha1 < 1
        if L==0
            % only case 1
            [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
                cases_1(w, lambda_d, lambda_other, L, gamma, theta, rho);
        else
            % all three cases are possible
            [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
                cases_123(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
                alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
        end
    elseif alpha1==1
        if L==0
            % only case 1
            [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
                cases_1(w, lambda_d, lambda_other, L, gamma, theta, rho);
        else
            % cases 2 and 3 are possible
            [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
                cases_23(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
                alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
        end
    else % alpha1 > 1
        % cases 2 and 3 are possible
        [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
            cases_23(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
            alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast);
    end
    
end

%---------------------------------------------------------------------%
%                        Some local functions                         %
%---------------------------------------------------------------------%

function [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
    cases_1(w, lambda_d, lambda_other, L, gamma, theta, rho)
    
    % case 1: n > 0, v = 0
    [c, x, n, v, l, d, ~] = ...
        case1(w, lambda_d, lambda_other, L, gamma, theta, rho);
    
    % case 3: zero variables
    c3 = 0; x3 = 0; n3 = 0; v3 = 0; l3 = 0; d3 = 0;
end

function [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
    cases_23(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast)

    % case 2: n > 0, v > 0
    [c, x, n, v, l, d, H] = ...
        case2(w, lambda_d, lambda_other, L, gamma, theta, rho, alpha1, alpha2);

    if flag_fast==1 && H~=-999
        c3 = 0; x3 = 0; n3 = 0; v3 = 0; l3 = 0; d3 = 0;
        return % exits without going to case 3
    end

    % case 3: n = 0, v > 0
    [c3, x3, n3, v3, l3, d3, H3] = ...
        case3(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess);

    % checks if case 3 is better than case 2
    if H3 > H
        c = c3; x = x3; n = n3; v = v3; l = l3; d = d3;
    end
end

function [c, x, n, v, l, d, c3, x3, n3, v3, l3, d3] = ...
    cases_123(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess, flag_fast)
    
    % case 2: n > 0, v > 0
    [c, x, n, v, l, d, H] = ...
        case2(w, lambda_d, lambda_other, L, gamma, theta, rho, alpha1, alpha2);

    if flag_fast==1 && H~=-999
        c3 = 0; x3 = 0; n3 = 0; v3 = 0; l3 = 0; d3 = 0;
        return % exits without going to case 3
    end
    
    % case 3: n = 0, v > 0
    [c3, x3, n3, v3, l3, d3, H3] = ...
        case3(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
        alpha1, alpha2, n_guess, v_guess, x_guess, l_guess);

    % checks if case 3 is better than case 2
    if H3 > H
        c = c3; x = x3; n = n3; v = v3; l = l3; d = d3;
        H = H3;
    end
    
    % case 1: n > 0, v = 0
    [c1, x1, n1, v1, l1, d1, H1] = ...
        case1(w, lambda_d, lambda_other, L, gamma, theta, rho);
    
    % checks if case 3 is better than last best case
    if H1 > H
        c = c1; x = x1; n = n1; v = v1; l = l1; d = d1;
    end
end

%---------------------------------------------------------------------%
%                        Case 1: n > 0, v = 0                         %
%---------------------------------------------------------------------%

function [c, x, n, v, l, d, H] = ...
    case1(w, lambda_d, lambda_other, L, gamma, theta, rho)
    
    % in this case, we just need the sum of the two lambdas
    lambda = lambda_d + lambda_other;

    % no teleworking in this case
    v = 0;

    % psis
    psi1 = gamma^-1 * (1 + (theta/(1-theta))^(1/(rho-1)) * w^( rho/(rho-1)));
    psi2 = gamma^-1 * (1 + ((1-theta)/theta)^(1/(rho-1)) * w^(-rho/(rho-1)));
    
    % coefficients for quadratic equation
    r2 = L * psi2 * (psi1 + psi2 + psi1 * psi2); % quadratic term
    r1 = psi1 + psi2 + psi1 * psi2 * (1 - L + lambda); % linear term
    r0 = - psi1; % constant term
    
    if L==0 % if quadratic term is zero
        
        l_vec = -r0/r1;
        
    else % if quadratic term is non-zero
        
        l_vec = roots([r2 r1 r0]); % solves for the root of a quadratic eq.

        if ~isreal(l_vec)
            error('solve_young, case 1: not real roots')
        end
        
    end
    
    % computes other choice variables
    d_vec = lambda  * l_vec * psi2 ./ (1 + L * l_vec * psi2);
    n_vec = 1 - l_vec - d_vec;
    x_vec = w * n_vec/(1 + psi1);
    c_vec = w * n_vec - x_vec;
    
    % checks if variables are inside bounds
    if_vec = c_vec >= 0 & ...
             x_vec >= 0 & ...
             n_vec >= 0 & n_vec <= 1 & ...
             l_vec >= 0 & l_vec <= 1 & ...
             d_vec >= 0 & d_vec <= 1;
    
    % number of possible solutions
    sum_if_vec = sum(if_vec);
    
    % checks some cases
    if sum_if_vec==1
        l = l_vec(if_vec);
        d = d_vec(if_vec);
        n = n_vec(if_vec);
        x = x_vec(if_vec);
        c = c_vec(if_vec);
    elseif sum_if_vec==2 % only one feasible solution
        error('solve_young, case 1: two possible feasible roots')
    elseif sum_if_vec==0 % some variable out of bounds
        error('solve_young, case 1: some variable out of bounds')
    else
        error('solve_young: sum_if_vec is not 2, 1, or 0.')
    end
    
    % value
    b = 0;
    H = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta) + (n + l) * L;
    
    % checks FOCs
    FOC = zeros(3, 1);
    FOC(1) = -1/(w*n - x) + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho);
    FOC(2) = w/(w*n - x) - lambda/(1-n-l) + L;
    FOC(3) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda/(1-n-l) + L;

    if max(abs(FOC)) > 1e-6
        error('FOC not zero')
    end

end

%---------------------------------------------------------------------%
%                        Case 2: n > 0, v > 0                         %
%---------------------------------------------------------------------%

function [c, x, n, v, l, d, H] = ...
    case2(w, lambda_d, lambda_other, L, gamma, theta, rho, alpha1, alpha2)

    % psis
    psi1 = gamma^-1*(1+(theta/(1-theta))^(1/(rho-1))*w^(rho/(rho-1)));
    psi3 = ((1 - theta)/(theta * w))^(1/(rho-1));
    
    one_plus_psi1 = 1 + psi1;
    psi3m1 = psi3^-1;
    
    % coefficients of initial polynomials
    A0 = one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w);
    A1 = (-one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w));

    B0 = alpha1*one_plus_psi1*psi1;
    B1 = - psi1*(alpha1*one_plus_psi1 + 2*alpha2*one_plus_psi1 - alpha1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) + alpha1^2*psi3m1*w + alpha1*lambda_d*psi1);
    B2 = psi1*(2*alpha2*one_plus_psi1 - alpha2*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) - alpha1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) + alpha2*lambda_d*psi1 + 3*alpha1*alpha2*psi3m1*w);
    B3 = psi1*(alpha2*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) - 2*alpha2^2*psi3m1*w);

    C0 = - lambda_d*one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w);
    C1 = - (-lambda_d*one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w));

    D0 = one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) - alpha1*one_plus_psi1^2 - alpha1*one_plus_psi1*psi3m1*w;
    D1 = (one_plus_psi1*(alpha1*one_plus_psi1 + 2*alpha2*one_plus_psi1 + alpha1^2*psi3m1*w + alpha1*lambda_d*psi1) - one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) - one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) + psi3m1*w*(alpha1*one_plus_psi1 + 2*alpha2*one_plus_psi1 - alpha1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) + alpha1^2*psi3m1*w + alpha1*lambda_d*psi1));
    D2 = (one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) - one_plus_psi1*(2*alpha2*one_plus_psi1 + alpha2*lambda_d*psi1 + 3*alpha1*alpha2*psi3m1*w) - psi3m1*w*(2*alpha2*one_plus_psi1 - alpha2*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) - alpha1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) + alpha2*lambda_d*psi1 + 3*alpha1*alpha2*psi3m1*w));
    D3 = (2*alpha2^2*one_plus_psi1*psi3m1*w - psi3m1*w*(alpha2*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) - 2*alpha2^2*psi3m1*w));

    E0 = - lambda_other*one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w);
    E1 = - (-lambda_other*one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w));

    F0 = one_plus_psi1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) - alpha1*one_plus_psi1^2 - alpha1*one_plus_psi1*psi3m1*w;
    F1 = (one_plus_psi1*(alpha1*one_plus_psi1 + 2*alpha2*one_plus_psi1 + alpha1^2*psi3m1*w + alpha1*lambda_d*psi1) - one_plus_psi1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) + psi3m1*w*(alpha1*one_plus_psi1 + 2*alpha2*one_plus_psi1 - alpha1*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) + alpha1^2*psi3m1*w + alpha1*lambda_d*psi1));
    F2 = (- one_plus_psi1*(2*alpha2*one_plus_psi1 + alpha2*lambda_d*psi1 + 3*alpha1*alpha2*psi3m1*w) - psi3m1*w*(2*alpha2*one_plus_psi1 - alpha2*(alpha1*one_plus_psi1 + lambda_d*psi1 + alpha1*psi3m1*w) - alpha1*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) + alpha2*lambda_d*psi1 + 3*alpha1*alpha2*psi3m1*w));
    F3 = (2*alpha2^2*one_plus_psi1*psi3m1*w - psi3m1*w*(alpha2*(2*alpha2*one_plus_psi1 + 2*alpha2*psi3m1*w) - 2*alpha2^2*psi3m1*w));

    % coefficients of final polynomial
    r0 = A0*D0*F0 + C0*B0*F0 + E0*B0*D0 + L*B0*D0*F0;
    r1 = A0*D0*F1 + A0*D1*F0 + A1*D0*F0 + C0*B0*F1 + C0*B1*F0 + C1*B0*F0 + E0*B0*D1 + E0*B1*D0 + E1*B0*D0 + L*B0*D0*F1 + L*B0*D1*F0 + L*B1*D0*F0;
    r2 = A0*D0*F2 + A0*D1*F1 + A0*D2*F0 + A1*D0*F1 + A1*D1*F0 + C0*B0*F2 + C0*B1*F1 + C0*B2*F0 + C1*B0*F1 + C1*B1*F0 + E0*B0*D2 + E0*B1*D1 + E0*B2*D0 + E1*B0*D1 + E1*B1*D0 + L*B0*D0*F2 + L*B0*D1*F1 + L*B0*D2*F0 + L*B1*D0*F1 + L*B1*D1*F0 + L*B2*D0*F0;
    r3 = A0*D0*F3 + A0*D1*F2 + A0*D2*F1 + A0*D3*F0 + A1*D0*F2 + A1*D1*F1 + A1*D2*F0 + C0*B0*F3 + C0*B1*F2 + C0*B2*F1 + C0*B3*F0 + C1*B0*F2 + C1*B1*F1 + C1*B2*F0 + E0*B0*D3 + E0*B1*D2 + E0*B2*D1 + E0*B3*D0 + E1*B0*D2 + E1*B1*D1 + E1*B2*D0 + L*B0*D0*F3 + L*B0*D1*F2 + L*B0*D2*F1 + L*B0*D3*F0 + L*B1*D0*F2 + L*B1*D1*F1 + L*B1*D2*F0 + L*B2*D0*F1 + L*B2*D1*F0 + L*B3*D0*F0;
    r4 = A0*D1*F3 + A0*D2*F2 + A0*D3*F1 + A1*D0*F3 + A1*D1*F2 + A1*D2*F1 + A1*D3*F0 + C0*B1*F3 + C0*B2*F2 + C0*B3*F1 + C1*B0*F3 + C1*B1*F2 + C1*B2*F1 + C1*B3*F0 + E0*B1*D3 + E0*B2*D2 + E0*B3*D1 + E1*B0*D3 + E1*B1*D2 + E1*B2*D1 + E1*B3*D0 + L*B0*D1*F3 + L*B0*D2*F2 + L*B0*D3*F1 + L*B1*D0*F3 + L*B1*D1*F2 + L*B1*D2*F1 + L*B1*D3*F0 + L*B2*D0*F2 + L*B2*D1*F1 + L*B2*D2*F0 + L*B3*D0*F1 + L*B3*D1*F0;
    r5 = A0*D2*F3 + A0*D3*F2 + A1*D1*F3 + A1*D2*F2 + A1*D3*F1 + C0*B2*F3 + C0*B3*F2 + C1*B1*F3 + C1*B2*F2 + C1*B3*F1 + E0*B2*D3 + E0*B3*D2 + E1*B1*D3 + E1*B2*D2 + E1*B3*D1 + L*B0*D2*F3 + L*B0*D3*F2 + L*B1*D1*F3 + L*B1*D2*F2 + L*B1*D3*F1 + L*B2*D0*F3 + L*B2*D1*F2 + L*B2*D2*F1 + L*B2*D3*F0 + L*B3*D0*F2 + L*B3*D1*F1 + L*B3*D2*F0;
    r6 = A0*D3*F3 + A1*D2*F3 + A1*D3*F2 + C0*B3*F3 + C1*B2*F3 + C1*B3*F2 + E0*B3*D3 + E1*B2*D3 + E1*B3*D2 + L*B0*D3*F3 + L*B1*D2*F3 + L*B1*D3*F2 + L*B2*D1*F3 + L*B2*D2*F2 + L*B2*D3*F1 + L*B3*D0*F3 + L*B3*D1*F2 + L*B3*D2*F1 + L*B3*D3*F0;
    r7 = A1*D3*F3 + C1*B3*F3 + E1*B3*D3 + L*B1*D3*F3 + L*B2*D2*F3 + L*B2*D3*F2 + L*B3*D1*F3 + L*B3*D2*F2 + L*B3*D3*F1;
    r8 = L*B2*D3*F3 + L*B3*D2*F3 + L*B3*D3*F2;
    r9 = L*B3*D3*F3;
    
    % calculates roots
    v_vec = roots([r9, r8, r7, r6, r5, r4, r3, r2, r1, r0]);
    
    % computes other variables
    tau_vec = alpha1 * v_vec - alpha2 * v_vec.^2;
    tau_prime_vec = alpha1 - 2 * alpha2 * v_vec;

    num = (1-v_vec) * one_plus_psi1 .* tau_prime_vec - ...
        tau_vec .* tau_prime_vec * psi3m1 * w - tau_vec * lambda_d * psi1;
    den = lambda_d * psi1 + (one_plus_psi1 + psi3m1*w)*tau_prime_vec;

    n_vec = num ./ den;

    x_vec = w * (n_vec + tau_vec)/(1+psi1);
    l_vec = psi3^-1 * x_vec;
    c_vec = w * (n_vec + tau_vec) - x_vec;
    d_vec = 1 - n_vec - l_vec - v_vec;
    
%     % checks if imaginary
%     if any(~real([v_vec; n_vec; x_vec; l_vec; c_vec; d_vec]))
%         error('solve young, case 2: imaginary')
%     end
    
    % checks if variables are inside bounds
    if_vec = c_vec >= 1e-10 & ...
             x_vec >= 1e-10 & ...
             n_vec >= 1e-10 & n_vec <= 1-1e-10 & ...
             v_vec >= 1e-10 & v_vec <= 1-1e-10 & ...
             l_vec >= 1e-10 & l_vec <= 1-1e-10 & ...
             d_vec >= 1e-10 & d_vec <= 1-1e-10;
    
    sum_if_vec = sum(if_vec);
    
    % checks possible cases of roots
    if sum_if_vec==0
        
        % not a single root that implies in variables inside bounds
        % therefore, there is no interior solution for n and v
        H = -999;
        
        c = 0;
        x = 0;
        n = 0;
        v = 0;
        l = 0;
        d = 0;
        
    elseif sum_if_vec==1
        
        % this is a candidate for solution
        c = c_vec(if_vec);
        x = x_vec(if_vec);
        n = n_vec(if_vec);
        v = v_vec(if_vec);
        l = l_vec(if_vec);
        d = d_vec(if_vec);
        tau_prime = tau_prime_vec(if_vec);
        
        % value
        b = 0;
        H = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta) + (n + l) * L;
        
        % checks if imaginary
        if ~isreal([c, x, n, v, l, d])
            error('solve young, case 2: imaginary')
        end
        
        % checks FOCs
%         FOC = zeros(4, 1);
%         FOC(1) = w/c - lambda_d/d - lambda_other/(d+v) + L; % [n]
%         FOC(2) = w*tau_prime/c - lambda_d/d; % [v]
%         FOC(3) = -1/c + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho); % [x]
%         FOC(4) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda_d/d - lambda_other/(d+v) + L; % [l]

%         FOC(1) = (w/c - L) / (lambda_d/d + lambda_other/(d+v)) - 1; % [n]
%         FOC(2) = (w*tau_prime/c) / (lambda_d/d) - 1; % [v]
%         FOC(3) = (1/c) / (gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho)) - 1; % [x]
%         FOC(4) = (gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - L) / (lambda_d/d + lambda_other/(d+v)); % [l]
        
%         if max(abs(FOC)) > 1e-8
%         if max(abs(FOC)) > 1e-3
%             error('FOC not zero')
%         end
        
    else
        
        % sum_if_vec>1
        error('solve young, case 2: sum_if_vec>1')
        
    end
    
end

%---------------------------------------------------------------------%
%                        Case 3: n = 0, v > 0                         %
%---------------------------------------------------------------------%

function [c, x, n, v, l, d, H] = ...
    case3(w, lambda_d, lambda_other, L, gamma, theta, rho, ...
    alpha1, alpha2, n_guess, v_guess, x_guess, l_guess)
    
    % no on-site working in this case
    n = 0;
    
    % deals with initial guesses
    if n_guess==0 && v_guess>0
        % starts guess vector
        guess = zeros(3, 1);
        
        % v
        ub = alpha1/(2 * alpha2);
        ub = min(ub, 1);
        lb = 0;
        guess(1) = asin((v_guess - (lb + ub)/2) * ((ub - lb)/2)^-1);
        
        % l
        ub = 1 - v_guess;
        lb = 0;
        guess(2) = asin((l_guess - (lb + ub)/2) * ((ub - lb)/2)^-1);
        
        % computes some variables using v
        tau = alpha1 * v_guess - alpha2 * v_guess^2;

        % transforms x to desired range
        ub = w * tau;
        lb = 0;
        guess(3) = asin((x_guess - (lb + ub)/2) * ((ub - lb)/2)^-1);
    else
        guess = ones(3, 1); % guesses something
    end
    
    % finds solution
    opts = optimset('Display', 'off');
    solution = fsolve(@fun, guess, opts);
    
    % gets all variables
    [~, c, d, l, x, v] = ...
        FOCs_case3(solution, alpha1, alpha2, w, lambda_d, ...
        lambda_other, L, theta, rho, gamma);
    
    tau_prime = alpha1 - 2 * alpha2 * v;
    
    % checks if variables are inside bounds
    flag_bounds_ok = c >= 0 && ...
                     x >= 0 && ...
                     n >= 0 && n <= 1 && ...
                     v >= 0 && v <= 1 && ...
                     l >= 0 && l <= 1 && ...
                     d >= 0 && d <= 1;
    
    if flag_bounds_ok==0
        error('case3 error out of bounds')
    end
    
    % value
    b = 0;
    H = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta) + (n + l) * L;
    
    % checks FOCs
    FOC = zeros(3, 1);
    FOC(1) = w*tau_prime/c - lambda_d/d; % [v]
    FOC(2) = -1/c + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho); % [x]
    FOC(3) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda_d/d - lambda_other/(d+v) + L; % [l]
    
    if max(abs(FOC)) > 1e-6
        error('FOC not zero')
    end
    
    % internal function
    function out = fun(in)
        out = FOCs_case3(in, alpha1, alpha2, w, lambda_d, ...
            lambda_other, L, theta, rho, gamma);
    end
end

function [out, c, d, l, x, v] = ...
    FOCs_case3(in, alpha1, alpha2, w, lambda_d, lambda_other, L, theta, rho, gamma)

    % open input variables
    v = in(1);
    l = in(2);
    x = in(3);
    
    % transforms v to desired range
    ub = alpha1/(2 * alpha2);
    ub = min(ub, 1);
    lb = 0;
    v = (lb + ub)/2 + sin(v) * (ub - lb)/2;
    
    % transforms l to desired range
    ub = 1 - v;
    lb = 0;
    l = (lb + ub)/2 + sin(l) * (ub - lb)/2;
    
    % computes some variables using v
    tau = alpha1 * v - alpha2 * v^2;
    tau_prime = alpha1 - 2 * alpha2 * v;
    
    % transforms x to desired range
    ub = w * tau;
    lb = 0;
    x = (lb + ub)/2 + sin(x) * (ub - lb)/2;
    
    % computes other variables
    c = w * tau - x;
    d = 1 - v - l;
    
    % evaluates FOCs
    out = zeros(3, 1);
    out(1) = w*tau_prime/c - lambda_d/d; % [v]
    out(2) = -1/c + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho); % [x]
    out(3) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda_d/d - lambda_other/(d+v) + L; % [l]
end

% function [out, c, d, l, x, v] = ...
%     FOCs_case3(in, alpha1, alpha2, w, lambda_d, lambda_other, L, theta, rho, gamma)
% 
%     % open input variables
%     v = in(1);
%     x = in(2);
%     
%     % transforms v to desired range
%     ub = alpha1/(2 * alpha2);
%     lb = 0;
%     v = ub * 0.5 * (lb + 1 + sin(v));
%     
%     % computes some variables using v
%     tau = alpha1 * v - alpha2 * v^2;
%     tau_prime = alpha1 - 2 * alpha2 * v;
%     
%     % transforms x to desired range
%     lb = 0;
%     ub = w * tau;
%     x = ub * 0.5 * (lb + 1 + sin(x));
%     
%     % computes other variables
%     c = w * tau - x;
%     l = (w * tau_prime * (1-v) - lambda_d * c)/(w * tau_prime);
%     
%     if l < 0 || l > 1-v
%         lb = 0;
%         ub = 1-v;
%         l = ub * 0.5 * (lb + 1 + sin(l));
%     end
%     
%     d = 1 - v - l;
%     
%     % evaluates FOCs
%     out = zeros(2, 1);
%     out(1) = -1/c + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho); % [x]
%     out(2) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda_d/d - lambda_other/(d+v) + L; % [l]
% end

% function [out, c, d, l, x, v] = ...
%     FOCs_case3(in, alpha1, alpha2, w, lambda_d, lambda_other, L, theta, rho, gamma)
% 
%     % open input variables
%     x_over_l = in(1);
%     l = in(2);
%     
%     % transforms to desired range (x_over_l)
%     x_over_l = x_over_l^2;
%     
%     % some variables
%     xi1 = gamma^-1 * (1 + (1-theta)/theta * (x_over_l)^-rho);
%     xi2 = gamma^-1 * (1 + theta/(1-theta) * (x_over_l)^ rho);
%     
%     % v as a function of l
%     v_fun = @(l) 1 - l - ...
%         lambda_d * (1/(l*xi2) - lambda_other/(1-l) + L)^-1;
%     
%     % finds lower bounds of l
%     ub_v = min(alpha1/(2*alpha2), 1-l); % PAREI AQUI. cansado demais... d tá dando negativo, mas todo o resto tá ok...
%     fun = @(l) v_fun(l)/ub_v - 1;
%     lb = fzero(fun, [1e-6, 1 - 1e-6]);
%     
%     % transforms to desired range (l)
%     ub = 1;
%     l = (lb + ub)/2 + sin(l) * (ub - lb)/2;
%     
%     % computes some variables
%     v = v_fun(l);
%     tau = alpha1 * v - alpha2 * v^2;
%     tau_prime = alpha1 - 2 * alpha2 * v;
%     x = w * tau / (1+xi1);
%     x_over_l_new = x/l;
%     c = w * tau - x;
%     d = 1 - v - l;
%     
%     % output (must be zero)
%     out = zeros(2, 1);
%     out(1) = x_over_l - x_over_l_new;
%     out(2) = tau_prime * (1+xi1) / (tau * xi1) - lambda_d/d;
% end






















