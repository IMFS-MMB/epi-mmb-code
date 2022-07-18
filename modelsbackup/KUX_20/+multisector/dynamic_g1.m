function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = multisector.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(33, 42);
g1(1,6)=T(28);
g1(1,7)=T(4)*T(2)*T(29)+T(2)*T(3)*1/y(21)*T(30);
g1(1,8)=(-(params(7)*params(10)*9*params(1)*y(10)*params(9)/T(6)));
g1(1,10)=(-(params(7)*params(10)*params(9)*9*params(1)*y(8)/T(6)));
g1(1,21)=T(2)*T(3)*T(30)*(-y(7))/(y(21)*y(21));
g1(2,21)=(-(1/(steady_state(16))));
g1(2,22)=1;
g1(3,6)=T(28);
g1(3,7)=T(8)*T(7)*T(29)+T(3)*T(7)*1/y(23)*T(31);
g1(3,8)=(-(params(7)*params(12)*9*params(1)*y(10)*params(11)/T(6)));
g1(3,10)=(-(params(7)*params(12)*params(11)*9*params(1)*y(8)/T(6)));
g1(3,23)=T(3)*T(7)*T(31)*(-y(7))/(y(23)*y(23));
g1(4,23)=(-(1/(steady_state(18))));
g1(4,24)=1;
g1(5,6)=T(28);
g1(5,7)=T(10)*T(9)*T(29)+T(3)*T(9)*1/y(25)*T(32);
g1(5,8)=(-(params(7)*params(14)*9*params(1)*y(10)*params(13)/T(6)));
g1(5,10)=(-(params(7)*params(14)*params(13)*9*params(1)*y(8)/T(6)));
g1(5,25)=T(3)*T(9)*T(32)*(-y(7))/(y(25)*y(25));
g1(6,25)=(-(1/(steady_state(20))));
g1(6,26)=1;
g1(7,6)=T(28);
g1(7,7)=T(12)*T(11)*T(29)+T(3)*T(11)*1/y(27)*T(33);
g1(7,8)=(-(params(7)*params(16)*9*params(1)*y(10)*params(15)/T(6)));
g1(7,10)=(-(params(7)*params(16)*params(15)*9*params(1)*y(8)/T(6)));
g1(7,27)=T(3)*T(11)*T(33)*(-y(7))/(y(27)*y(27));
g1(8,27)=(-(1/(steady_state(22))));
g1(8,28)=1;
g1(9,6)=T(28);
g1(9,7)=T(14)*T(13)*T(29)+T(3)*T(13)*1/y(29)*T(34);
g1(9,8)=(-(params(7)*params(18)*9*params(1)*y(10)*params(17)/T(6)));
g1(9,10)=(-(params(7)*params(18)*params(17)*9*params(1)*y(8)/T(6)));
g1(9,29)=T(3)*T(13)*T(34)*(-y(7))/(y(29)*y(29));
g1(10,29)=(-(1/(steady_state(24))));
g1(10,30)=1;
g1(11,6)=T(28);
g1(11,7)=T(16)*T(15)*T(29)+T(3)*T(15)*1/y(31)*T(35);
g1(11,8)=(-(params(7)*params(20)*9*params(1)*y(10)*params(19)/T(6)));
g1(11,10)=(-(params(7)*params(20)*params(19)*9*params(1)*y(8)/T(6)));
g1(11,31)=T(3)*T(15)*T(35)*(-y(7))/(y(31)*y(31));
g1(12,31)=(-(1/(steady_state(26))));
g1(12,32)=1;
g1(13,6)=T(28);
g1(13,7)=T(18)*T(17)*T(29)+T(3)*T(17)*1/y(33)*T(36);
g1(13,8)=(-(params(7)*params(22)*9*params(1)*y(10)*params(21)/T(6)));
g1(13,10)=(-(params(7)*params(22)*params(21)*9*params(1)*y(8)/T(6)));
g1(13,33)=T(3)*T(17)*T(36)*(-y(7))/(y(33)*y(33));
g1(14,33)=(-(1/(steady_state(28))));
g1(14,34)=1;
g1(15,6)=T(28);
g1(15,7)=T(20)*T(19)*T(29)+T(3)*T(19)*1/y(35)*T(37);
g1(15,8)=(-(params(7)*params(24)*9*params(1)*y(10)*params(23)/T(6)));
g1(15,10)=(-(params(7)*params(24)*params(23)*9*params(1)*y(8)/T(6)));
g1(15,35)=T(3)*T(19)*T(37)*(-y(7))/(y(35)*y(35));
g1(16,35)=(-(1/(steady_state(30))));
g1(16,36)=1;
g1(17,6)=T(28);
g1(17,7)=T(22)*T(21)*T(29)+T(3)*T(21)*1/y(37)*T(38);
g1(17,8)=(-(params(7)*params(26)*9*params(1)*y(10)*params(25)/T(6)));
g1(17,10)=(-(params(7)*params(26)*params(25)*9*params(1)*y(8)/T(6)));
g1(17,37)=T(3)*T(21)*T(38)*(-y(7))/(y(37)*y(37));
g1(18,37)=(-(1/(steady_state(32))));
g1(18,38)=1;
g1(19,6)=(-params(7));
g1(19,21)=1;
g1(19,23)=1;
g1(19,25)=1;
g1(19,27)=1;
g1(19,29)=1;
g1(19,31)=1;
g1(19,33)=1;
g1(19,35)=1;
g1(19,37)=1;
g1(20,7)=1;
g1(20,21)=(-(T(2)*getPowerDeriv(y(21),1-T(1),1)*T(39)));
g1(20,23)=(-(T(39)*T(7)*getPowerDeriv(y(23),1-T(1),1)));
g1(20,25)=(-(T(39)*T(9)*getPowerDeriv(y(25),1-T(1),1)));
g1(20,27)=(-(T(39)*T(11)*getPowerDeriv(y(27),1-T(1),1)));
g1(20,29)=(-(T(39)*T(13)*getPowerDeriv(y(29),1-T(1),1)));
g1(20,31)=(-(T(39)*T(15)*getPowerDeriv(y(31),1-T(1),1)));
g1(20,33)=(-(T(39)*T(17)*getPowerDeriv(y(33),1-T(1),1)));
g1(20,35)=(-(T(39)*T(19)*getPowerDeriv(y(35),1-T(1),1)));
g1(20,37)=(-(T(39)*T(21)*getPowerDeriv(y(37),1-T(1),1)));
g1(21,9)=1;
g1(21,10)=(-(params(2)+params(1)*9*T(24)*T(26)));
g1(21,21)=(-(T(25)*params(10)*params(9)));
g1(21,23)=(-(T(25)*params(12)*params(11)));
g1(21,25)=(-(T(25)*params(14)*params(13)));
g1(21,27)=(-(T(25)*params(16)*params(15)));
g1(21,29)=(-(T(25)*params(18)*params(17)));
g1(21,31)=(-(T(25)*params(20)*params(19)));
g1(21,33)=(-(T(25)*params(22)*params(21)));
g1(21,35)=(-(T(25)*params(24)*params(23)));
g1(21,37)=(-(T(25)*params(26)*params(25)));
g1(22,8)=1;
g1(22,39)=params(8);
g1(22,40)=(-params(8));
g1(23,17)=1;
g1(23,39)=(-(params(8)*(1-params(4)-params(3))));
g1(23,41)=(-(params(8)*params(3)));
g1(24,19)=1;
g1(24,41)=(-params(8));
g1(25,6)=params(6)/2*2*y(6);
g1(25,7)=(-T(3));
g1(25,9)=(-(params(8)*(y(39)-y(40))));
g1(25,39)=(-(y(9)*params(8)));
g1(25,18)=1;
g1(25,40)=(-(params(8)*(1-y(9))));
g1(26,9)=(-y(12));
g1(26,11)=1;
g1(26,12)=(-y(9));
g1(27,1)=(-(1-params(3)-params(4)));
g1(27,10)=1;
g1(27,3)=(-1);
g1(27,12)=1;
g1(28,1)=(-(1-params(3)-params(4)));
g1(28,10)=1;
g1(28,2)=(-1);
g1(28,42)=(-1);
g1(29,1)=(-params(3));
g1(29,4)=(-1);
g1(29,13)=1;
g1(30,1)=(-params(4));
g1(30,5)=(-1);
g1(30,14)=1;
g1(31,5)=1;
g1(31,14)=(-1);
g1(31,20)=1;
g1(32,7)=(-(y(12)/T(24)));
g1(32,10)=(-1);
g1(32,12)=(-(y(7)/T(24)));
g1(32,13)=(-1);
g1(32,15)=1;
g1(33,6)=(-(T(6)*y(12)));
g1(33,10)=(-1);
g1(33,12)=(-(y(6)*T(6)));
g1(33,13)=(-1);
g1(33,16)=1;

end
