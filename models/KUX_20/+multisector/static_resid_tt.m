function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
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

assert(length(T) >= 27);

T(1) = 1/params(5);
T(2) = params(10)^T(1);
T(3) = 1/y(2);
T(4) = (y(2)/y(16))^T(1);
T(5) = params(6)/params(7)*y(1);
T(6) = sqrt(params(6));
T(7) = params(12)^T(1);
T(8) = (y(2)/y(18))^T(1);
T(9) = params(14)^T(1);
T(10) = (y(2)/y(20))^T(1);
T(11) = params(16)^T(1);
T(12) = (y(2)/y(22))^T(1);
T(13) = params(18)^T(1);
T(14) = (y(2)/y(24))^T(1);
T(15) = params(20)^T(1);
T(16) = (y(2)/y(26))^T(1);
T(17) = params(22)^T(1);
T(18) = (y(2)/y(28))^T(1);
T(19) = params(24)^T(1);
T(20) = (y(2)/y(30))^T(1);
T(21) = params(26)^T(1);
T(22) = (y(2)/y(32))^T(1);
T(23) = T(2)*y(16)^(1-T(1))+T(7)*y(18)^(1-T(1))+T(9)*y(20)^(1-T(1))+T(11)*y(22)^(1-T(1))+T(13)*y(24)^(1-T(1))+T(15)*y(26)^(1-T(1))+T(17)*y(28)^(1-T(1))+T(19)*y(30)^(1-T(1))+T(21)*y(32)^(1-T(1));
T(24) = params(7)/T(6);
T(25) = y(5)*params(1)*9*T(24);
T(26) = params(10)*y(16)*params(9)+params(12)*y(18)*params(11)+params(14)*y(20)*params(13)+params(16)*y(22)*params(15)+params(18)*y(24)*params(17)+params(20)*y(26)*params(19)+params(22)*y(28)*params(21)+params(24)*y(30)*params(23)+params(26)*y(32)*params(25);
T(27) = log(T(24))-params(6)/2*(1/T(6))^2;

end
