function out = imperfect_DynPath(N, MUs, par)

muc = MUs(:,1);
muw_s = MUs(:,2);
muw_i_tilde = MUs(:,3);
muw_r_tilde = MUs(:,4);
if par.ntypes > 3
    muw_i_star = MUs(:,5);
    if par.ntypes > 4 
        muw_r_star = MUs(:,6);
    end
end

N = reshape(N, par.H, par.ntypes);
ns          = N(:,1);
ni_tilde    = N(:,2);
nr_tilde    = N(:,3);
if par.ntypes > 3
    ni_star = N(:,4);
    if par.ntypes > 4
        nr_star = N(:,5);
    end
end


% Recovered and identified
lambda_br_tilde = par.theta/(par.A*par.phi_r_tilde)*nr_tilde./(1 + muw_r_tilde);
cr_tilde = (1 + muc).^(-1).*lambda_br_tilde.^(-1);
ur_tilde = log(cr_tilde) - par.theta/2.*(nr_tilde).^(2);

Gamma = (1 + muc).*cr_tilde - (1 + muw_r_tilde).*par.A*par.phi_r_tilde.*nr_tilde;

% Recovered, but not-identified
if par.ntypes == 5
    lambda_br_star = par.theta/(par.A*par.phi_r_star)*nr_star./(1 + muw_r_star);
    cr_star = (1 + muc).^(-1).*lambda_br_star.^(-1);
    ur_star = log(cr_star) - par.theta/2.*(nr_star).^(2);
end

% Infected and identified
lambda_bi_tilde = par.theta/(par.phi_tilde*par.A)*ni_tilde./(1 + muw_i_tilde);
ci_tilde = ((1 + muc).*lambda_bi_tilde).^(-1);
ui_tilde = log(ci_tilde) - par.theta/2*(ni_tilde).^(2);

% Infected, but not-identified
if par.chi < 1
    lambda_bi_star = par.theta/(par.phi_star*par.A)*ni_star./(1 + muw_i_star);
    ci_star = ((1 + muc).*lambda_bi_star).^(-1);
    ui_star = log(ci_star) - par.theta/2*(ni_star).^(2);
end

cs = ((1 + muw_s).*par.A.*ns + Gamma)./(1 + muc);
us = log(cs) - par.theta/2*(ns).^(2);

I = zeros(par.H,1); I(1) = par.I0;
S = zeros(par.H,1); S(1) = par.S0;
T = zeros(par.H,1); 
T1 = zeros(par.H,1); T2 = zeros(par.H,1); T3 = zeros(par.H,1); T4 = zeros(par.H,1); T5 = zeros(par.H,1); 
R = zeros(par.H,1);
D = zeros(par.H,1);
ICi = zeros(par.H,1);
INi = zeros(par.H,1);

pi_dt_tilde = zeros(par.H,1);
pi_dt_star  = zeros(par.H,1);

% Going forward
for t = 1:par.H
    % Different from ERT
    if par.chi < 1
        ICi(t) = I(t).*( par.chi*ci_tilde(t) + (1 - par.chi)*ci_star(t) );
        INi(t) = I(t).*( par.chi*ni_tilde(t) + (1 - par.chi)*ni_star(t) );
    elseif par.chi == 1
        ICi(t) = I(t).*ci_tilde(t);
        INi(t) = I(t).*ni_tilde(t);
    end
        
    T1(t) = par.pi_s1*(S(t)*cs(t))*(ICi(t));
    T2(t) = par.pi_s2*(S(t)*ns(t))*(INi(t));
    T3(t) = par.pi_s3*S(t)*I(t);
    T4(t) = par.pi_s4*(S(t)*cs(t))*(INi(t));
    T5(t) = par.pi_s5*(S(t)*ns(t))*(ICi(t));
    
    T(t) = T1(t) + T2(t) + T3(t) + T4(t) + T5(t);
        
    pi_dt_star(t) = par.pi_d + par.kappa*I(t)^2;
    pi_dt_tilde(t)  = par.pi_d + ( (par.chi - par.chi_bar)*I(t) )^2;
    
    S(t+1) = S(t) - T(t);
    I(t+1) = I(t) + T(t) - (par.pi_r + ( par.chi*pi_dt_tilde(t) + (1 - par.chi)*pi_dt_star(t) ))*I(t);
    R(t+1) = R(t) + par.pi_r*I(t);
    D(t+1) = D(t) + ( par.chi*pi_dt_tilde(t) + (1 - par.chi)*pi_dt_star(t) )*I(t);
end

% Terminal values
urss = log(par.c_ss) - par.theta/2*(par.n_ss)^2;
Ur_tilde = zeros(par.H+1,1); 
Ur_tilde_ss = urss/(1 - par.beta);
Ur_tilde(end) = Ur_tilde_ss;

if par.chi_AB < 1
    Ur_star = zeros(par.H+1,1); 
    Ur_star_ss = (urss + par.beta*par.chi_AB*Ur_tilde_ss)/(1 - par.beta*(1 - par.chi_AB));
    Ur_star(end) = Ur_star_ss;
end

ci_tilde_ss = par.phi_tilde*par.A*par.n_ss;
ui_tilde_ss = log(ci_tilde_ss) - par.theta/2*par.n_ss^2;
Ui_tilde_ss = ( ui_tilde_ss + par.beta*par.pi_r*Ur_tilde_ss )/( 1 - par.beta*(1 - par.pi_r - par.pi_d) );
Ui_tilde = zeros(par.H+1,1);
Ui_tilde(end) = Ui_tilde_ss;

if par.chi < 1
    ci_star_ss = par.phi_star*par.A*par.n_ss;
    ui_star_ss = log(ci_star_ss) - par.theta/2*par.n_ss^2;
    if par.chi_AB < 1
        Ui_star_ss = ( ui_star_ss + par.beta*par.pi_r*( (1 - par.chi_AB)*Ur_star_ss + par.chi_AB*Ur_tilde_ss ) )/( 1 - par.beta*(1 - par.pi_r - par.pi_d) );
    else
        Ui_star_ss = ( ui_star_ss + par.beta*par.pi_r*( Ur_tilde_ss ) )/( 1 - par.beta*(1 - par.pi_r - par.pi_d) );
    end
    Ui_star = zeros(par.H+1,1);
    Ui_star(end) = Ui_star_ss;
end

Us = zeros(par.H+1,1);
Us(end) = Ur_tilde_ss;

Us_true = zeros(par.H+1,1);
Us_true(end) = Ur_tilde_ss;

tau = T./S(1:par.H);

tau_hat = par.psi*tau;

% Going backward
for t = par.H:-1:1
    Ur_tilde(t) =  ur_tilde(t) + par.beta*Ur_tilde(t+1);
    
    Ui_tilde(t) = ui_tilde(t) + par.beta*( (1 - par.pi_r - pi_dt_tilde(t))*Ui_tilde(t+1) + par.pi_r*Ur_tilde(t+1) );

    if par.chi < 1
        if par.chi_AB < 1
            Ur_star(t) =  ur_star(t) + par.beta*( (1 - par.chi_AB)*Ur_star(t+1) + par.chi_AB*Ur_tilde(t+1) );
            Ui_star(t)  = ui_star(t) + par.beta*( (1 - par.pi_r - pi_dt_star(t))*( Ui_star(t+1) ) ...
                                       + par.pi_r*( (1 - par.chi_AB)*Ur_star(t+1) + par.chi_AB*Ur_tilde(t+1) ) );
        else 
            Ui_star(t)  = ui_star(t) + par.beta*( (1 - par.pi_r - pi_dt_star(t))*Ui_star(t+1) + par.pi_r*Ur_tilde(t+1) );
        end
        
        Us(t) = us(t) + par.beta*( (1 - tau_hat(t))*Us(t+1) + tau_hat(t)*( (1 - par.chi)*Ui_star(t+1) + par.chi*Ui_tilde(t+1) ) );
        Us_true(t) = us(t) + par.beta*( (1 - tau(t))*Us_true(t+1) + tau(t)*( (1 - par.chi)*Ui_star(t+1) + par.chi*Ui_tilde(t+1) ) );

    else
        Us(t)       = us(t) + par.beta*( (1 - tau_hat(t))*Us(t+1) + tau_hat(t)*( Ui_tilde(t+1) ) );
        Us_true(t)  = us(t) + par.beta*( (1 - tau(t))*Us_true(t+1) + tau(t)*( Ui_tilde(t+1) ) );
    end

    
end

if par.chi < 1
    lambda_tau_hat = par.beta*((1 - par.chi)*Ui_star(2:end,1) + par.chi*Ui_tilde(2:end,1) - Us(2:end,1));
else
    lambda_tau_hat = par.beta*(Ui_tilde(2:end,1) - Us(2:end,1));
end

lambda_bs = ( cs.^(-1) + lambda_tau_hat.*( par.pi_s1*par.psi*I(1:par.H).*ci_tilde + par.pi_s4*par.psi*I(1:par.H).*ni_tilde ) )./(1 + muc);

% Putting everything together
out.lambda_br_tilde = lambda_br_tilde;
out.cr_tilde = cr_tilde;
out.ur_tilde = ur_tilde;
out.Ur_tilde = Ur_tilde(1:par.H);

out.Gamma = Gamma;
out.lambda_bi_tilde = lambda_bi_tilde;
out.ci_tilde = ci_tilde;
out.ui_tilde = ui_tilde;
out.cs = cs;
out.us = us;
out.I = I(1:par.H);
out.S = S(1:par.H);
out.T = T;
out.T1 = T1;
out.T2 = T2;
out.T3 = T3;
out.T4 = T4;
out.T5 = T5;
out.R = R(1:par.H);
out.D = D(1:par.H);

out.Ui_tilde = Ui_tilde(1:par.H);
out.Us = Us(1:par.H);
out.Us_true = Us_true(1:par.H);

out.tau = tau;
out.tau_hat = tau_hat;
out.lambda_tau_hat = lambda_tau_hat;
out.lambda_bs = lambda_bs;

out.ICi = ICi;
out.INi = INi;
out.pi_dt_star = pi_dt_star;
out.pi_dt_tilde = pi_dt_tilde;

if par.ntypes == 4
    out.lambda_bi_star = lambda_bi_star;
    out.ci_star = ci_star;
    out.ui_star = ui_star;
    out.Ui_star = Ui_star(1:par.H);
end

if par.ntypes == 5
    out.lambda_br_star = lambda_br_star;
    out.cr_star = cr_star;
    out.ur_star = ur_star;
    out.Ur_star = Ur_star(1:par.H);
end


