% V_in is the vector V_h_private that we want to use as target for the CEV

function factor_c_eq = c_eq(V_in, i_age)
    
    % defines target V
    assignin('base', 'V_h_private_target', V_in)
    assignin('base', 'i_age_c_eq', i_age)
    
    % finds consumption equivalent
    fzero(@fun, 1);
    factor_c_eq = evalin('base', 'factor_c_eq');
    
end

function distance = fun(in)

    % converts factor to positive
    factor_c_eq = in^2;

    % computes new V with higher consumption in all periods
    assignin('base', 'factor_c_eq', factor_c_eq)
    evalin('base', 'compute_V_factor_c_eq')

    % measures distance
    i_age = evalin('base', 'i_age_c_eq');
    V_h_private_new = evalin('base', 'V_h_private_new');
    V_h_private_target = evalin('base', 'V_h_private_target');
    distance = V_h_private_new(i_age, 1)/V_h_private_target(i_age, 1) - 1;
    
end