% -------------------------------------------------------------------------
% This version distinguishes psi and chi:
% psi measures the degree of misperception 
%   underestimation: 0 < psi < 1
%   overestimation: 1 < psi  
% chi measures a fration of identified/quarantined
% 
% When chi = 1, solve for 3 types of agents: ns, ni_tilde, and nr_tilde.
% When chi < 1, solve for 4 types of agents: ns, ni_tilde, ni_star, and nr_tilde.
% -------------------------------------------------------------------------

opt = optimoptions(@fsolve, 'Display', 'iter', 'TolFun', 1e-9);

% Need to have par in the same Workspace
if exist('par', 'var')~=1
    error('Need the parameter structure "par".')
end
% Need to have MUs in the same Workspace
if exist('MUs', 'var')~=1
    error('Need the parameter structure "MUs".')
end

if par.chi == 1
    par.ntypes = 3;
elseif par.chi < 1 && par.chi_AB == 1
    par.ntypes = 4;
elseif par.chi < 1 && par.chi_AB < 1
    par.ntypes = 5;
end

% Solve for the optimal path
% - Initial guess
N0 = par.n_ss*ones(par.ntypes*par.H,1);
% - Objective function
ObjFun = @(x) imperfect_ObjFun(x, MUs, par);
[N, ~, exitflag] = fsolve(ObjFun, N0, opt);
if exitflag ~= 1
    error('Incomplete Solution.')
%     warning('Incomplete Solution.')
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

cr_tilde = X.cr_tilde;
ci_tilde = X.ci_tilde;
cs = X.cs;		
I  = X.I; 		
S  = X.S; 		
R  = X.R; 		
D  = X.D;

ICi = X.ICi;

if par.chi == 1
    Hi = I.*par.phi_tilde.*ni_tilde;
else
    ci_star = X.ci_star;

    Hi = par.chi*I.*par.phi_tilde.*ni_tilde + (1 - par.chi)*I.*par.phi_star.*ni_star;
end

if par.chi_AB == 1
    RCr = R.*cr_tilde;
    Hr = R.*nr_tilde*par.phi_r_tilde;
else
    cr_star = X.cr_star;
    
    RCr = par.chi_AB*R.*cr_tilde + (1 - par.chi_AB)*R.*cr_star;
    Hr = par.chi_AB*R.*nr_tilde*par.phi_r_tilde + (1 - par.chi_AB)*R.*nr_star*par.phi_r_star;
end

% Aggregate consumption and hours
C = S.*cs + ICi + RCr;
H = S.*ns + Hi + Hr;

if par.chi == 1
    U0 = X.S(1)*X.Us_true(1) + X.I(1)*( X.Ui_tilde(1) );
else
    U0 = X.S(1)*X.Us_true(1) + X.I(1)*( par.chi*X.Ui_tilde(1) + (1 - par.chi)*X.Ui_star(1) );
end

pi_dt_tilde = X.pi_dt_tilde;
pi_dt_star = X.pi_dt_star;
 
if saveflag == 1
    savelist = {'par'; 'N'; 'ns'; 'cs'; 'ci_tilde'; 'cr_tilde'; ...
        'ni_tilde';  'nr_tilde';  'I'; 'S'; 'R'; 'D'; 'C'; 'H'; 'U0'; 'X'; ...
        'pi_dt_tilde'; 'pi_dt_star'};
    
    if par.ntypes > 3
        savelist = [savelist; {'ci_star'; 'ni_star'}];
        if par.ntypes > 4
            savelist = [savelist; {'cr_star'; 'nr_star'}];
        end
    end
    
    savelist = cellfun(@(x) [x, '_v' num2str(spec) ' = ' x ';'], savelist, 'UniformOutput', false);
    savelist = [savelist; ['save Results_Imperfect_v' num2str(spec) ' *_v' num2str(spec)]];
    
    for i = 1:length(savelist)
        eval(savelist{i})
    end
end


