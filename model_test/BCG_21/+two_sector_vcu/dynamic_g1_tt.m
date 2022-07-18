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

assert(length(T) >= 19);

T = two_sector_vcu.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(14) = (-1);
T(15) = getPowerDeriv(params(9)*y(116)/y(119),T(5),1);
T(16) = getPowerDeriv(y(116)*(1-params(9))/y(118),T(5),1);
T(17) = getPowerDeriv(T(10),1+params(3),1);
T(18) = getPowerDeriv(y(3)*y(135),params(1),1);
T(19) = getPowerDeriv(y(133)*y(134),1-params(1),1);

end