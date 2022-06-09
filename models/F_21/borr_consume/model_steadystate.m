function [ys,params,check] = model_steadystate(ys,exo,M_,options_)
% function [ys,params,check] = full_model_uspackage_steadystate(ys,exo,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here

%% Enter model equations here
load('detss.mat','detss')
take_logs = 1;
varnames = fieldnames(detss);
var_nolog = {'nua'; 'nun'; 'nuu'; 'T'; 'rot'; 'govt_wage'; 'transfer'; 'bailout'; 'mu'; 'lbd'; 'spend'; 
    'spend_G'; 'spend_tau_l'; 'spend_varsigma'; 'spend_govt_wage'; 'spend_transfer'; 'income'; 
    'ci'; 'Gi'; 'm'};
for jj = 1:numel(varnames) 
    if sum(strcmp(varnames{jj}, var_nolog)) >= 1 || ~take_logs
        eval([varnames{jj},' = detss.',varnames{jj},';']);
    else
        eval([varnames{jj},' = log(detss.',varnames{jj},');']);
    end  
end

%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end

end