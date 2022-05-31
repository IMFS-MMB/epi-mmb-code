function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = multisector.static_g1_tt(T, y, x, params);
end
g1 = zeros(33, 33);
g1(1,1)=T(28);
g1(1,2)=T(4)*T(2)*T(29)+T(2)*T(3)*1/y(16)*T(30);
g1(1,3)=(-(params(7)*params(10)*9*params(1)*y(5)*params(9)/T(6)));
g1(1,5)=(-(params(7)*params(10)*params(9)*9*params(1)*y(3)/T(6)));
g1(1,16)=T(2)*T(3)*T(30)*(-y(2))/(y(16)*y(16));
g1(2,16)=(-(((y(16))-y(16))/((y(16))*(y(16)))));
g1(2,17)=1;
g1(3,1)=T(28);
g1(3,2)=T(8)*T(7)*T(29)+T(3)*T(7)*1/y(18)*T(31);
g1(3,3)=(-(params(7)*params(12)*9*params(1)*y(5)*params(11)/T(6)));
g1(3,5)=(-(params(7)*params(12)*params(11)*9*params(1)*y(3)/T(6)));
g1(3,18)=T(3)*T(7)*T(31)*(-y(2))/(y(18)*y(18));
g1(4,18)=(-(((y(18))-y(18))/((y(18))*(y(18)))));
g1(4,19)=1;
g1(5,1)=T(28);
g1(5,2)=T(10)*T(9)*T(29)+T(3)*T(9)*1/y(20)*T(32);
g1(5,3)=(-(params(7)*params(14)*9*params(1)*y(5)*params(13)/T(6)));
g1(5,5)=(-(params(7)*params(14)*params(13)*9*params(1)*y(3)/T(6)));
g1(5,20)=T(3)*T(9)*T(32)*(-y(2))/(y(20)*y(20));
g1(6,20)=(-(((y(20))-y(20))/((y(20))*(y(20)))));
g1(6,21)=1;
g1(7,1)=T(28);
g1(7,2)=T(12)*T(11)*T(29)+T(3)*T(11)*1/y(22)*T(33);
g1(7,3)=(-(params(7)*params(16)*9*params(1)*y(5)*params(15)/T(6)));
g1(7,5)=(-(params(7)*params(16)*params(15)*9*params(1)*y(3)/T(6)));
g1(7,22)=T(3)*T(11)*T(33)*(-y(2))/(y(22)*y(22));
g1(8,22)=(-(((y(22))-y(22))/((y(22))*(y(22)))));
g1(8,23)=1;
g1(9,1)=T(28);
g1(9,2)=T(14)*T(13)*T(29)+T(3)*T(13)*1/y(24)*T(34);
g1(9,3)=(-(params(7)*params(18)*9*params(1)*y(5)*params(17)/T(6)));
g1(9,5)=(-(params(7)*params(18)*params(17)*9*params(1)*y(3)/T(6)));
g1(9,24)=T(3)*T(13)*T(34)*(-y(2))/(y(24)*y(24));
g1(10,24)=(-(((y(24))-y(24))/((y(24))*(y(24)))));
g1(10,25)=1;
g1(11,1)=T(28);
g1(11,2)=T(16)*T(15)*T(29)+T(3)*T(15)*1/y(26)*T(35);
g1(11,3)=(-(params(7)*params(20)*9*params(1)*y(5)*params(19)/T(6)));
g1(11,5)=(-(params(7)*params(20)*params(19)*9*params(1)*y(3)/T(6)));
g1(11,26)=T(3)*T(15)*T(35)*(-y(2))/(y(26)*y(26));
g1(12,26)=(-(((y(26))-y(26))/((y(26))*(y(26)))));
g1(12,27)=1;
g1(13,1)=T(28);
g1(13,2)=T(18)*T(17)*T(29)+T(3)*T(17)*1/y(28)*T(36);
g1(13,3)=(-(params(7)*params(22)*9*params(1)*y(5)*params(21)/T(6)));
g1(13,5)=(-(params(7)*params(22)*params(21)*9*params(1)*y(3)/T(6)));
g1(13,28)=T(3)*T(17)*T(36)*(-y(2))/(y(28)*y(28));
g1(14,28)=(-(((y(28))-y(28))/((y(28))*(y(28)))));
g1(14,29)=1;
g1(15,1)=T(28);
g1(15,2)=T(20)*T(19)*T(29)+T(3)*T(19)*1/y(30)*T(37);
g1(15,3)=(-(params(7)*params(24)*9*params(1)*y(5)*params(23)/T(6)));
g1(15,5)=(-(params(7)*params(24)*params(23)*9*params(1)*y(3)/T(6)));
g1(15,30)=T(3)*T(19)*T(37)*(-y(2))/(y(30)*y(30));
g1(16,30)=(-(((y(30))-y(30))/((y(30))*(y(30)))));
g1(16,31)=1;
g1(17,1)=T(28);
g1(17,2)=T(22)*T(21)*T(29)+T(3)*T(21)*1/y(32)*T(38);
g1(17,3)=(-(params(7)*params(26)*9*params(1)*y(5)*params(25)/T(6)));
g1(17,5)=(-(params(7)*params(26)*params(25)*9*params(1)*y(3)/T(6)));
g1(17,32)=T(3)*T(21)*T(38)*(-y(2))/(y(32)*y(32));
g1(18,32)=(-(((y(32))-y(32))/((y(32))*(y(32)))));
g1(18,33)=1;
g1(19,1)=(-params(7));
g1(19,16)=1;
g1(19,18)=1;
g1(19,20)=1;
g1(19,22)=1;
g1(19,24)=1;
g1(19,26)=1;
g1(19,28)=1;
g1(19,30)=1;
g1(19,32)=1;
g1(20,2)=1;
g1(20,16)=(-(T(2)*getPowerDeriv(y(16),1-T(1),1)*T(39)));
g1(20,18)=(-(T(39)*T(7)*getPowerDeriv(y(18),1-T(1),1)));
g1(20,20)=(-(T(39)*T(9)*getPowerDeriv(y(20),1-T(1),1)));
g1(20,22)=(-(T(39)*T(11)*getPowerDeriv(y(22),1-T(1),1)));
g1(20,24)=(-(T(39)*T(13)*getPowerDeriv(y(24),1-T(1),1)));
g1(20,26)=(-(T(39)*T(15)*getPowerDeriv(y(26),1-T(1),1)));
g1(20,28)=(-(T(39)*T(17)*getPowerDeriv(y(28),1-T(1),1)));
g1(20,30)=(-(T(39)*T(19)*getPowerDeriv(y(30),1-T(1),1)));
g1(20,32)=(-(T(39)*T(21)*getPowerDeriv(y(32),1-T(1),1)));
g1(21,4)=1;
g1(21,5)=(-(params(2)+params(1)*9*T(24)*T(26)));
g1(21,16)=(-(T(25)*params(10)*params(9)));
g1(21,18)=(-(T(25)*params(12)*params(11)));
g1(21,20)=(-(T(25)*params(14)*params(13)));
g1(21,22)=(-(T(25)*params(16)*params(15)));
g1(21,24)=(-(T(25)*params(18)*params(17)));
g1(21,26)=(-(T(25)*params(20)*params(19)));
g1(21,28)=(-(T(25)*params(22)*params(21)));
g1(21,30)=(-(T(25)*params(24)*params(23)));
g1(21,32)=(-(T(25)*params(26)*params(25)));
g1(22,3)=1;
g1(22,12)=params(8);
g1(22,13)=(-params(8));
g1(23,12)=1-params(8)*(1-params(4)-params(3));
g1(23,14)=(-(params(8)*params(3)));
g1(24,14)=1-params(8);
g1(25,1)=params(6)/2*2*y(1);
g1(25,2)=(-T(3));
g1(25,4)=(-(params(8)*(y(12)-y(13))));
g1(25,12)=(-(y(4)*params(8)));
g1(25,13)=1-params(8)*(1-y(4));
g1(26,4)=(-y(7));
g1(26,6)=1;
g1(26,7)=(-y(4));
g1(27,5)=1-(1-params(3)-params(4));
g1(28,5)=1-(1-params(3)-params(4));
g1(28,6)=(-1);
g1(29,5)=(-params(3));
g1(30,5)=(-params(4));
g1(31,15)=1;
g1(32,2)=(-(y(7)/T(24)));
g1(32,5)=(-1);
g1(32,7)=(-(y(2)/T(24)));
g1(32,8)=(-1);
g1(32,10)=1;
g1(33,1)=(-(T(6)*y(7)));
g1(33,5)=(-1);
g1(33,7)=(-(y(1)*T(6)));
g1(33,8)=(-1);
g1(33,11)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
