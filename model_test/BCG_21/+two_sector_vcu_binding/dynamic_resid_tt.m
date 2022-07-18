function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 9);

T(1) = params(2)*params(6)/(y(247)-y(120)*params(2));
T(2) = params(3)/(1+params(3));
T(3) = (params(9)*y(115)/y(118))^T(2);
T(4) = (1-params(9))^T(2);
T(5) = 1/(1+params(3));
T(6) = params(9)^T(2);
T(7) = T(4)*y(117)^T(5)+T(6)*y(118)^T(5);
T(8) = (y(2)*y(134))^params(1);
T(9) = (y(132)*y(133))^(1-params(1));

end