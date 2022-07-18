function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = multisector.static_resid_tt(T, y, x, params);
end
residual = zeros(33, 1);
lhs = T(2)*T(3)*T(4);
rhs = T(5)+params(7)*params(10)*9*params(1)*y(5)*y(3)*params(9)/T(6);
residual(1) = lhs - rhs;
lhs = y(17);
rhs = y(16)/(y(16))-1;
residual(2) = lhs - rhs;
lhs = T(3)*T(7)*T(8);
rhs = T(5)+params(7)*params(12)*9*params(1)*y(5)*y(3)*params(11)/T(6);
residual(3) = lhs - rhs;
lhs = y(19);
rhs = y(18)/(y(18))-1;
residual(4) = lhs - rhs;
lhs = T(3)*T(9)*T(10);
rhs = T(5)+params(7)*params(14)*9*params(1)*y(5)*y(3)*params(13)/T(6);
residual(5) = lhs - rhs;
lhs = y(21);
rhs = y(20)/(y(20))-1;
residual(6) = lhs - rhs;
lhs = T(3)*T(11)*T(12);
rhs = T(5)+params(7)*params(16)*9*params(1)*y(5)*y(3)*params(15)/T(6);
residual(7) = lhs - rhs;
lhs = y(23);
rhs = y(22)/(y(22))-1;
residual(8) = lhs - rhs;
lhs = T(3)*T(13)*T(14);
rhs = T(5)+params(7)*params(18)*9*params(1)*y(5)*y(3)*params(17)/T(6);
residual(9) = lhs - rhs;
lhs = y(25);
rhs = y(24)/(y(24))-1;
residual(10) = lhs - rhs;
lhs = T(3)*T(15)*T(16);
rhs = T(5)+params(7)*params(20)*9*params(1)*y(5)*y(3)*params(19)/T(6);
residual(11) = lhs - rhs;
lhs = y(27);
rhs = y(26)/(y(26))-1;
residual(12) = lhs - rhs;
lhs = T(3)*T(17)*T(18);
rhs = T(5)+params(7)*params(22)*9*params(1)*y(5)*y(3)*params(21)/T(6);
residual(13) = lhs - rhs;
lhs = y(29);
rhs = y(28)/(y(28))-1;
residual(14) = lhs - rhs;
lhs = T(3)*T(19)*T(20);
rhs = T(5)+params(7)*params(24)*9*params(1)*y(5)*y(3)*params(23)/T(6);
residual(15) = lhs - rhs;
lhs = y(31);
rhs = y(30)/(y(30))-1;
residual(16) = lhs - rhs;
lhs = T(3)*T(21)*T(22);
rhs = T(5)+params(7)*params(26)*9*params(1)*y(5)*y(3)*params(25)/T(6);
residual(17) = lhs - rhs;
lhs = y(33);
rhs = y(32)/(y(32))-1;
residual(18) = lhs - rhs;
lhs = y(32)+y(30)+y(28)+y(26)+y(24)+y(22)+y(20)+y(16)+y(18);
rhs = params(7)*y(1);
residual(19) = lhs - rhs;
lhs = y(2);
rhs = T(23)^(params(5)/(params(5)-1));
residual(20) = lhs - rhs;
lhs = y(4);
rhs = T(25)*T(26)+y(5)*params(2);
residual(21) = lhs - rhs;
lhs = y(3);
rhs = (-params(8))*(y(12)-y(13));
residual(22) = lhs - rhs;
lhs = y(12);
rhs = T(27)+params(8)*(y(12)*(1-params(4)-params(3))+params(3)*y(14));
residual(23) = lhs - rhs;
lhs = y(14);
rhs = T(27)+params(8)*y(14);
residual(24) = lhs - rhs;
lhs = y(13);
rhs = log(y(2))-params(6)/2*y(1)^2+params(8)*(y(13)*(1-y(4))+y(4)*y(12));
residual(25) = lhs - rhs;
lhs = y(6);
rhs = y(4)*y(7);
residual(26) = lhs - rhs;
residual(27) = y(5)-y(5)*(1-params(3)-params(4));
lhs = y(5);
rhs = y(6)+y(5)*(1-params(3)-params(4))+x(1);
residual(28) = lhs - rhs;
lhs = y(8);
rhs = y(8)+y(5)*params(3);
residual(29) = lhs - rhs;
lhs = y(9);
rhs = y(9)+y(5)*params(4);
residual(30) = lhs - rhs;
residual(31) = y(15);
lhs = y(10);
rhs = (y(2)*y(7)+params(7)*(y(5)+y(8))/T(6))/T(24)-1;
residual(32) = lhs - rhs;
lhs = y(11);
rhs = T(6)*(y(1)*y(7)+(y(5)+y(8))/T(6))-1;
residual(33) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
