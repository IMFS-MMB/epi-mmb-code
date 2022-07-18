function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = multisector.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(33, 1);
lhs = T(2)*T(3)*T(4);
rhs = T(5)+params(7)*params(10)*9*params(1)*y(10)*y(8)*params(9)/T(6);
residual(1) = lhs - rhs;
lhs = y(22);
rhs = y(21)/(steady_state(16))-1;
residual(2) = lhs - rhs;
lhs = T(3)*T(7)*T(8);
rhs = T(5)+params(7)*params(12)*9*params(1)*y(10)*y(8)*params(11)/T(6);
residual(3) = lhs - rhs;
lhs = y(24);
rhs = y(23)/(steady_state(18))-1;
residual(4) = lhs - rhs;
lhs = T(3)*T(9)*T(10);
rhs = T(5)+params(7)*params(14)*9*params(1)*y(10)*y(8)*params(13)/T(6);
residual(5) = lhs - rhs;
lhs = y(26);
rhs = y(25)/(steady_state(20))-1;
residual(6) = lhs - rhs;
lhs = T(3)*T(11)*T(12);
rhs = T(5)+params(7)*params(16)*9*params(1)*y(10)*y(8)*params(15)/T(6);
residual(7) = lhs - rhs;
lhs = y(28);
rhs = y(27)/(steady_state(22))-1;
residual(8) = lhs - rhs;
lhs = T(3)*T(13)*T(14);
rhs = T(5)+params(7)*params(18)*9*params(1)*y(10)*y(8)*params(17)/T(6);
residual(9) = lhs - rhs;
lhs = y(30);
rhs = y(29)/(steady_state(24))-1;
residual(10) = lhs - rhs;
lhs = T(3)*T(15)*T(16);
rhs = T(5)+params(7)*params(20)*9*params(1)*y(10)*y(8)*params(19)/T(6);
residual(11) = lhs - rhs;
lhs = y(32);
rhs = y(31)/(steady_state(26))-1;
residual(12) = lhs - rhs;
lhs = T(3)*T(17)*T(18);
rhs = T(5)+params(7)*params(22)*9*params(1)*y(10)*y(8)*params(21)/T(6);
residual(13) = lhs - rhs;
lhs = y(34);
rhs = y(33)/(steady_state(28))-1;
residual(14) = lhs - rhs;
lhs = T(3)*T(19)*T(20);
rhs = T(5)+params(7)*params(24)*9*params(1)*y(10)*y(8)*params(23)/T(6);
residual(15) = lhs - rhs;
lhs = y(36);
rhs = y(35)/(steady_state(30))-1;
residual(16) = lhs - rhs;
lhs = T(3)*T(21)*T(22);
rhs = T(5)+params(7)*params(26)*9*params(1)*y(10)*y(8)*params(25)/T(6);
residual(17) = lhs - rhs;
lhs = y(38);
rhs = y(37)/(steady_state(32))-1;
residual(18) = lhs - rhs;
lhs = y(37)+y(35)+y(33)+y(31)+y(29)+y(27)+y(25)+y(21)+y(23);
rhs = params(7)*y(6);
residual(19) = lhs - rhs;
lhs = y(7);
rhs = T(23)^(params(5)/(params(5)-1));
residual(20) = lhs - rhs;
lhs = y(9);
rhs = T(25)*T(26)+y(10)*params(2);
residual(21) = lhs - rhs;
lhs = y(8);
rhs = (-params(8))*(y(39)-y(40));
residual(22) = lhs - rhs;
lhs = y(17);
rhs = T(27)+params(8)*(y(39)*(1-params(4)-params(3))+params(3)*y(41));
residual(23) = lhs - rhs;
lhs = y(19);
rhs = T(27)+params(8)*y(41);
residual(24) = lhs - rhs;
lhs = y(18);
rhs = log(y(7))-params(6)/2*y(6)^2+params(8)*(y(40)*(1-y(9))+y(9)*y(39));
residual(25) = lhs - rhs;
lhs = y(11);
rhs = y(9)*y(12);
residual(26) = lhs - rhs;
residual(27) = y(10)+y(12)-y(3)-(1-params(3)-params(4))*y(1);
lhs = y(10);
rhs = (1-params(3)-params(4))*y(1)+y(2)+x(it_, 1);
residual(28) = lhs - rhs;
lhs = y(13);
rhs = y(4)+params(3)*y(1);
residual(29) = lhs - rhs;
lhs = y(14);
rhs = y(5)+params(4)*y(1);
residual(30) = lhs - rhs;
lhs = y(20);
rhs = y(14)-y(5);
residual(31) = lhs - rhs;
lhs = y(15);
rhs = (y(7)*y(12)+params(7)*(y(10)+y(13))/T(6))/T(24)-1;
residual(32) = lhs - rhs;
lhs = y(16);
rhs = T(6)*(y(6)*y(12)+(y(10)+y(13))/T(6))-1;
residual(33) = lhs - rhs;

end
