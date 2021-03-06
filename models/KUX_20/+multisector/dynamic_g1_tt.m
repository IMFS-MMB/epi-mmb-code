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

assert(length(T) >= 39);

T = multisector.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(28) = (-(params(6)/params(7)));
T(29) = (-1)/(y(7)*y(7));
T(30) = getPowerDeriv(y(7)/y(21),T(1),1);
T(31) = getPowerDeriv(y(7)/y(23),T(1),1);
T(32) = getPowerDeriv(y(7)/y(25),T(1),1);
T(33) = getPowerDeriv(y(7)/y(27),T(1),1);
T(34) = getPowerDeriv(y(7)/y(29),T(1),1);
T(35) = getPowerDeriv(y(7)/y(31),T(1),1);
T(36) = getPowerDeriv(y(7)/y(33),T(1),1);
T(37) = getPowerDeriv(y(7)/y(35),T(1),1);
T(38) = getPowerDeriv(y(7)/y(37),T(1),1);
T(39) = getPowerDeriv(T(23),params(5)/(params(5)-1),1);

end
