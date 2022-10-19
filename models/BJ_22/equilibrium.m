

% initial guess for q
if flag_fix_q==1
    q = q_fixed;
    q_m1_fix = q_m1_fixed;
else
    q = zeros(T, 1);
    q(1:20) = linspace(0, 0.04, 20);
    q(21:40) = linspace(0.04, 0, 20);

    q_m1 = Pi_q * initial_infected;
end

% guess for choices
if flag_use_guess==0
    if alpha1==1
        guess = ones(T, 2) * N;
    else
        guess = ones(T, 3) * (1 - h_over_n_plus_h) * N;
    end
end

% bounds for choices
if alpha1==1
    lb = zeros(T, 2);
    ub = ones(T, 2) * N;
else
    lb = zeros(T, 3);
    ub = ones(T, 3) * N;
end

% general equilibrium tolerance
tol_ge = 1e-6; %1e-7;
norm = tol_ge + 1;

% general equilibrium loop
while norm > tol_ge
    
    % gets firm's choices
    if flag_fixed_choices==1
        % fixed choices
        solution = guess;
    else
        % solves problem of firm
%         solution = fmincon(@compute_profits, guess, [], [], [], [], lb, ub, [], opts_fmincon);
        
%         opts = optimset('MaxFunEvals', 1e10, 'MaxIter', 1e10);
%         solution = fminsearch(@compute_profits, guess, opts);

max_iter = intmax;
tol = 1e-10;

opts_fmincon = optimoptions(@fmincon, ...
          'MaxIter', 1e10, ...
          'MaxFunEvals', 1e10, ...
          'FunctionTolerance', tol, ...
          'OptimalityTolerance', tol, ...
          'StepTolerance', tol, ...
          'ConstraintTolerance', tol, ...
          'TolCon', tol, ...
          'Display', 'off');

%         opts = optimset('MaxFunEvals', 1e10, 'MaxIter', 1e10);
%         solution = fmincon(@compute_profits, guess, [],[],[],[],lb,ub,[], opts);
        problem = createOptimProblem('fmincon',...
                'objective',@(x)compute_profits(x),...
                'x0', guess, 'lb', lb, 'ub', ub, ...
                'options', opts_fmincon);
        solution = fmincon(problem);
%                 problem = createOptimProblem('fmincon',...
%                 'objective',@(x)compute_profits(x),...
%                 'x0', solution, 'lb', lb, 'ub', ub, ...
%                 'options', opts_fmincon);
%         gs = GlobalSearch('Display','off');
%         ms = MultiStart('Display','off');
%         solution = run(gs, problem);
%         solution = run(ms,problem, 200);
%         solution = patternsearch(@compute_profits, guess, [],[],[],[],lb,ub,[], opts);
    end
    
    % use solution of this iteration as guess for next iteration
    guess = solution;
    
    % gets n. of infectious agents to compute new q
    [~, outstruct] = compute_profits(solution);
    
    % new q
    if flag_fix_q==1
        q_m1_new = q_m1_fixed;
        q_new = q_fixed;
    else
        infectious = outstruct.a + outstruct.n_tilde + outstruct.m_tilde;
        q_m1_new = Pi_q * initial_infected;
        q_new = Pi_q * infectious;
    end
    
    % defines A if it is endogenous
    if flag_endog_A_1==1
        s = outstruct.s;
        s_max = max(s);
        A_new = 1 + s * (A_min - 1)/s_max;
    elseif flag_endog_A_2==1
        s = outstruct.s;
        s_max = max(s);
        A_new = 1 + slope_A * s;
    else
        A_new = A;
    end
    
    % defines delta_n if it is endogenous
    if flag_endog_delta_n==1
        s = outstruct.s;
        s_max = max(s);
        delta_n_new =  s/s_max * delta_n_max + (s_max - s)/s_max * 1;
    else
        delta_n_new = delta_n;
    end
    
    % checks convergence of equilibrium variables
    V = [q_m1; q; A; delta_n];
    V_new = [q_m1_new; q_new; A_new; delta_n_new];
    
    norm = max(abs(V - V_new));
    
    % updates equilibrium variables
    q_m1 = q_m1_new;
    q = q_new;
    A = A_new;
    delta_n = delta_n_new;
    
    % display
    fprintf('norm: %.6e\n', norm)
    
end




