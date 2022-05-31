function [c, x, n, v, l, d] = ...
    solve_old(lambda_d, lambda_other, gamma, theta, rho, w_bar, L, l_guess, x_guess)
    
    if L > 1e-6
        error('solve_old: L>0')
    end
    
    % relevant lambda in this case is the sum of lambda_d and lambda_other
    lambda = lambda_d + lambda_other;
    
    % finds equilibrium l_over_x
    if l_guess > 0 && x_guess > 0
        guess = sqrt(l_guess/x_guess);
    else
        guess = 1;
    end
    
    sqrt_l_over_x = fzero(@foc, guess);
    
    % computes variables
    l_over_x = sqrt_l_over_x^2;
    x = x_fun(l_over_x);
    l = l_fun(l_over_x);
    d = 1 - l;
    c = w_bar - x;
    n = 0;
    v = 0;
    
    % checks if inside bounds
    if c < 0 || x < 0 || l < 0 || l > 1 || d < 0 || d > 1
        error('solve_old: out of bounds')
    end
    
    % checks FOCs
    FOC = zeros(2, 1);
    FOC(1) = -1/(w_bar - x) + gamma*theta*x^(rho-1)/(theta*x^rho+(1-theta)*l^rho);
    FOC(2) = gamma*(1-theta)*l^(rho-1)/(theta*x^rho+(1-theta)*l^rho) - lambda/(1-l) + L;

    if max(abs(FOC)) > 1e-6
        error('FOC not zero')
    end
    
    % internal functions
    function out = x_fun(l_over_x)
        out = w_bar * (gamma^-1 * ...
            (1 + (1-theta)/theta * l_over_x^rho) + 1)^-1;
    end

    function out = l_fun(l_over_x)
        out = l_over_x * x_fun(l_over_x);
    end
    
    function out = foc(l_over_x)
        % transforms to positive
        l_over_x = l_over_x^2;
        
        % computes FOC
        out = gamma * (1 - theta) * l_fun(l_over_x)^(rho-1) / ...
              (theta * x_fun(l_over_x)^rho + ...
              (1 - theta) * l_fun(l_over_x)^rho) - ...
              lambda / (1 - l_fun(l_over_x)) + L;
    end
end




















