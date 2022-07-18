function u = utility(c, x, l, d, v, lambda_d, lambda_other, b, gamma, rho, theta)

    % leisure good
    g = (theta*x^rho + (1-theta)*l^rho)^(1/rho);

    % utility level
    u = log(c) + gamma * log(g) + lambda_d * log(d) + lambda_other * log(d+v) + b;

end