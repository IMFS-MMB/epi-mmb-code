function F = imperfect_ObjFun(N, MUs, par)

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

X = imperfect_DynPath(N, MUs, par);
cr_tilde 	= X.cr_tilde;
Gamma 		= X.Gamma; 	

ICi = X.ICi;
ci_tilde 	= X.ci_tilde;		

cs 			= X.cs;		
I 			= X.I; 		
S 			= X.S; 		
R 			= X.R; 		

lambda_tau_hat  = X.lambda_tau_hat;
lambda_bs   = X.lambda_bs;

F = zeros(par.H, par.ntypes);

F(:,1) = (1 + muc).*ci_tilde - (1 + muw_i_tilde).*par.phi_tilde*par.A.*ni_tilde - Gamma;
F(:,2) = muc.*(S.*cs + ICi + R.*cr_tilde) ...
            - ( muw_s.*S.*ns + muw_i_tilde.*I.*par.phi_tilde.*ni_tilde + muw_r_tilde.*R.*nr_tilde )*par.A ...
            - Gamma.*(S + I + R);
F(:,3) = -par.theta*ns + (1 + muw_s).*par.A.*lambda_bs + lambda_tau_hat.*(par.pi_s2.*par.psi*I.*ni_tilde + par.pi_s5.*par.psi.*I.*ci_tilde);

if par.ntypes > 3
    ci_star = X.ci_star;
    F(:,2) = muc.*(S.*cs + ICi + R.*cr_tilde) ...
                - ( muw_s.*S.*ns + muw_i_tilde.*par.chi.*I.*par.phi_tilde.*ni_tilde + muw_i_star.*(1 - par.chi).*I.*par.phi_star.*ni_star ...
                + muw_r_tilde.*R.*nr_tilde )*par.A ...
                - Gamma.*(S + I + R);
    F(:,4) = (1 + muc).*ci_star - (1 + muw_i_star).*par.phi_star*par.A.*ni_star - Gamma;
    if par.ntypes > 4
        cr_star = X.cr_star;
        F(:,2) = muc.*(S.*cs + ICi + par.chi_AB*R.*cr_tilde + (1 - par.chi_AB)*R.*cr_star) ...
                - ( muw_s.*S.*ns + muw_i_tilde.*par.chi.*I.*par.phi_tilde.*ni_tilde + muw_i_star.*(1 - par.chi).*I.*par.phi_star.*ni_star ...
                + muw_r_tilde.*par.chi_AB.*R.*nr_tilde + muw_r_star.*(1 - par.chi_AB).*R.*nr_star )*par.A ...
                - Gamma.*(S + I + R);
        F(:,5) = (1 + muc).*cr_star - (1 + muw_r_star).*par.phi_r_star*par.A.*nr_star - Gamma;
    end
end

F = F(:);
