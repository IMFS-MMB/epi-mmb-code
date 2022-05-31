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

assert(length(T) >= 39);

T = multisector.static_resid_tt(T, y, x, params);

T(28) = (-(params(6)/params(7)));
T(29) = (-1)/(y(2)*y(2));
T(30) = getPowerDeriv(y(2)/y(16),T(1),1);
T(31) = getPowerDeriv(y(2)/y(18),T(1),1);
T(32) = getPowerDeriv(y(2)/y(20),T(1),1);
T(33) = getPowerDeriv(y(2)/y(22),T(1),1);
T(34) = getPowerDeriv(y(2)/y(24),T(1),1);
T(35) = getPowerDeriv(y(2)/y(26),T(1),1);
T(36) = getPowerDeriv(y(2)/y(28),T(1),1);
T(37) = getPowerDeriv(y(2)/y(30),T(1),1);
T(38) = getPowerDeriv(y(2)/y(32),T(1),1);
T(39) = getPowerDeriv(T(23),params(5)/(params(5)-1),1);

end
