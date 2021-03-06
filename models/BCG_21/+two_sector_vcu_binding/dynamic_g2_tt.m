function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 29);

T = two_sector_vcu_binding.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(24) = getPowerDeriv(params(9)*y(115)/y(118),T(2),2);
T(25) = T(18)*T(18)*T(24)+T(12)*(-((-(params(9)*y(115)))*(y(118)+y(118))))/(y(118)*y(118)*y(118)*y(118));
T(26) = getPowerDeriv(y(115)*(1-params(9))/y(117),T(2),2);
T(27) = getPowerDeriv(T(7),1+params(3),2);
T(28) = getPowerDeriv(y(2)*y(134),params(1),2);
T(29) = getPowerDeriv(y(132)*y(133),1-params(1),2);

end
