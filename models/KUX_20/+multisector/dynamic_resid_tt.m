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

assert(length(T) >= 27);

T(1) = 1/params(5);
T(2) = params(10)^T(1);
T(3) = 1/y(7);
T(4) = (y(7)/y(21))^T(1);
T(5) = params(6)/params(7)*y(6);
T(6) = sqrt(params(6));
T(7) = params(12)^T(1);
T(8) = (y(7)/y(23))^T(1);
T(9) = params(14)^T(1);
T(10) = (y(7)/y(25))^T(1);
T(11) = params(16)^T(1);
T(12) = (y(7)/y(27))^T(1);
T(13) = params(18)^T(1);
T(14) = (y(7)/y(29))^T(1);
T(15) = params(20)^T(1);
T(16) = (y(7)/y(31))^T(1);
T(17) = params(22)^T(1);
T(18) = (y(7)/y(33))^T(1);
T(19) = params(24)^T(1);
T(20) = (y(7)/y(35))^T(1);
T(21) = params(26)^T(1);
T(22) = (y(7)/y(37))^T(1);
T(23) = T(2)*y(21)^(1-T(1))+T(7)*y(23)^(1-T(1))+T(9)*y(25)^(1-T(1))+T(11)*y(27)^(1-T(1))+T(13)*y(29)^(1-T(1))+T(15)*y(31)^(1-T(1))+T(17)*y(33)^(1-T(1))+T(19)*y(35)^(1-T(1))+T(21)*y(37)^(1-T(1));
T(24) = params(7)/T(6);
T(25) = y(10)*params(1)*9*T(24);
T(26) = params(10)*y(21)*params(9)+params(12)*y(23)*params(11)+params(14)*y(25)*params(13)+params(16)*y(27)*params(15)+params(18)*y(29)*params(17)+params(20)*y(31)*params(19)+params(22)*y(33)*params(21)+params(24)*y(35)*params(23)+params(26)*y(37)*params(25);
T(27) = log(T(24))-params(6)/2*(1/T(6))^2;

end
