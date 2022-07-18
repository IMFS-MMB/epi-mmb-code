function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 23);

T = two_sector_vcu_binding.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(10) = (-1);
T(11) = params(9)/y(118);
T(12) = getPowerDeriv(params(9)*y(115)/y(118),T(2),1);
T(13) = (1-params(9))/y(117);
T(14) = getPowerDeriv(y(115)*(1-params(9))/y(117),T(2),1);
T(15) = (-(y(115)*(1-params(9))))/(y(117)*y(117));
T(16) = T(4)*getPowerDeriv(y(117),T(5),1);
T(17) = getPowerDeriv(T(7),1+params(3),1);
T(18) = (-(params(9)*y(115)))/(y(118)*y(118));
T(19) = params(1)*T(12)*T(18);
T(20) = T(6)*getPowerDeriv(y(118),T(5),1);
T(21) = y(252)*(-(params(2)*params(6)*(-params(2))))/((y(247)-y(120)*params(2))*(y(247)-y(120)*params(2)));
T(22) = getPowerDeriv(y(2)*y(134),params(1),1);
T(23) = getPowerDeriv(y(132)*y(133),1-params(1),1);

end
