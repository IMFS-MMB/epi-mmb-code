function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 14);

T = two_sector_vcu_binding.static_resid_tt(T, y, x, params);

T(9) = (-1);
T(10) = getPowerDeriv(params(9)*y(1)/y(4),T(1),1);
T(11) = getPowerDeriv(y(1)*(1-params(9))/y(3),T(1),1);
T(12) = getPowerDeriv(T(6),1+params(3),1);
T(13) = getPowerDeriv(y(20)*y(9),params(1),1);
T(14) = getPowerDeriv(y(18)*y(19),1-params(1),1);

end
