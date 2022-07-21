function [residual, g1, g2, g3] = BKM_22_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(65, 1);
T(1) = (1-y(75))^params(26)*params(3);
T(2) = y(48)^(-params(10));
T(3) = y(52)^(-params(10));
T(4) = y(55)/params(16)*y(66);
T(5) = y(42)^(-params(10));
T(6) = (y(7)+y(1)+y(6))*T(5)+y(2)*T(2)+y(3)*T(3);
T(7) = params(12)*params(18)*y(97)^params(17);
T(8) = params(12)*params(18)*y(97)^(params(17)-1);
T(9) = params(18)*y(56)^params(17);
T(10) = (y(56)/(steady_state(31)))^params(20);
T(11) = (steady_state(41))/(steady_state(6))+x(it_, 8);
T(12) = y(66)/y(31)/T(11);
T(13) = T(12)^params(21);
T(14) = T(10)*T(13);
T(15) = exp(y(27))^params(22);
T(16) = T(14)*T(15);
T(17) = exp(y(27)-y(2))^params(23);
lhs = y(36);
rhs = min(3*params(4),params(4)*(1+y(2)/params(6)));
residual(1) = lhs - rhs;
lhs = y(37);
rhs = 0.3888888888888889-y(36);
residual(2) = lhs - rhs;
lhs = y(30);
rhs = y(1)*y(45);
residual(3) = lhs - rhs;
lhs = y(26);
rhs = y(1)-y(30);
residual(4) = lhs - rhs;
lhs = y(27);
rhs = y(2)+y(30)*params(8)-y(2)*(y(36)+y(37));
residual(5) = lhs - rhs;
lhs = y(32);
rhs = y(6)+y(30)*(1-params(8))-(y(36)+y(37))*y(6);
residual(6) = lhs - rhs;
lhs = y(28);
rhs = y(3)+y(2)*y(37)+params(9)*y(7);
residual(7) = lhs - rhs;
lhs = y(33);
rhs = (y(36)+y(37))*y(6)+y(7)*(1-params(9));
residual(8) = lhs - rhs;
lhs = y(29);
rhs = y(4)+y(36)*y(2);
residual(9) = lhs - rhs;
lhs = y(31);
rhs = y(5)-y(36)*y(2);
residual(10) = lhs - rhs;
lhs = y(45);
rhs = params(1)*y(42)*(y(2)*params(13)*y(48)+y(42)*y(6)*params(7))+params(2)*y(43)*(y(2)*params(15)*y(49)+y(6)*params(7)*y(43))+T(1)*(y(2)+y(6)*params(7));
residual(11) = lhs - rhs;
lhs = y(35);
rhs = ((1-params(9))*y(9)+(y(36)+y(37))*y(8))/(1-y(45));
residual(12) = lhs - rhs;
lhs = y(34);
rhs = y(8)*(1-y(37)-y(36))/(1-y(45))+y(45)*(1-params(8))/(1-y(45));
residual(13) = lhs - rhs;
lhs = y(38);
rhs = 1-y(39)-y(40);
residual(14) = lhs - rhs;
lhs = y(39);
rhs = y(45)*params(8)/(y(9)+1+y(8));
residual(15) = lhs - rhs;
lhs = y(40);
rhs = params(9)*y(9)/(y(9)+1+y(8));
residual(16) = lhs - rhs;
lhs = y(41);
rhs = log(y(42))+params(11)*log(1-y(43))+params(12)*((1-y(39)-y(40))*y(91)+y(39)*y(93)+y(40)*y(95));
residual(17) = lhs - rhs;
lhs = (y(7)+y(1)+y(6))*(y(42)+y(62)-y(43)*y(55)-y(61));
rhs = y(11)*(y(15)-params(8)*(y(15)-y(1))+y(18)+(1-params(9))*y(19))*y(14)/y(56);
residual(18) = lhs - rhs;
lhs = y(42)^(-1);
rhs = (1+y(75))*y(44)-(y(2)*params(13)*y(48)+y(42)*y(6)*params(7))*params(1)*y(46);
residual(19) = lhs - rhs;
lhs = params(8)*params(12)*(y(93)-y(91));
rhs = (y(9)+1+y(8))*y(46);
residual(20) = lhs - rhs;
lhs = params(11)/(1-y(43));
rhs = y(44)*y(55)*(1-y(78))+(y(2)*params(15)*y(49)+y(6)*params(7)*y(43))*params(2)*y(46);
residual(21) = lhs - rhs;
lhs = y(44);
rhs = (1-y(39)-y(40))*params(12)*y(92)*y(65)/y(97)+y(65)*y(39)*params(12)*y(94)/y(97)+y(65)*y(40)*params(12)*y(96)/y(97);
residual(22) = lhs - rhs;
lhs = y(47);
rhs = log(y(48))+params(11)*log(1-y(49))+params(12)*((1-y(37)-y(36))*y(93)+y(37)*y(95)+y(36)*params(19));
residual(23) = lhs - rhs;
lhs = y(2)*(y(48)+y(64)-y(49)*y(55)*params(14)-y(61));
rhs = y(14)*((1-y(37))*y(13)*y(16)+(y(15)-y(1))*params(8)*y(11))/y(56);
residual(24) = lhs - rhs;
lhs = T(2);
rhs = (1+y(76))*y(50);
residual(25) = lhs - rhs;
lhs = y(49);
rhs = max(1-params(11)/(y(55)*params(14)*y(50)),0);
residual(26) = lhs - rhs;
lhs = y(50);
rhs = y(65)*y(94)*(1-y(37)-y(36))*params(12)/y(97)/(1-y(36))+y(65)*y(96)*y(37)*params(12)/y(97)/(1-y(36));
residual(27) = lhs - rhs;
lhs = y(51);
rhs = log(y(52))+params(11)*log(1-y(53))+params(12)*y(95);
residual(28) = lhs - rhs;
lhs = y(3)*(y(52)+y(63)-y(55)*y(53)-y(61));
rhs = y(14)*(y(12)*y(17)+y(16)*y(37)*y(13)+y(19)*params(9)*y(11))/y(56);
residual(29) = lhs - rhs;
lhs = T(3);
rhs = (1+y(77))*y(54);
residual(30) = lhs - rhs;
lhs = params(11)/(1-y(53));
rhs = y(54)*y(55)*(1-y(80));
residual(31) = lhs - rhs;
lhs = y(54);
rhs = y(65)*params(12)*y(96)/y(97);
residual(32) = lhs - rhs;
lhs = y(57);
rhs = y(58)/y(59);
residual(33) = lhs - rhs;
lhs = y(58)/(1-params(12)*params(18));
rhs = T(4)*T(6)+T(7)*y(98)/(1-params(12)*params(18));
residual(34) = lhs - rhs;
lhs = y(59)/(1-params(12)*params(18));
rhs = y(66)*T(6)+T(8)*y(99)/(1-params(12)*params(18));
residual(35) = lhs - rhs;
lhs = 1;
rhs = params(18)*y(56)^(params(17)-1)+(1-params(18))*y(57)^(1-params(17));
residual(36) = lhs - rhs;
lhs = y(60);
rhs = T(9)*y(10)+(1-params(18))*y(57)^(-params(17));
residual(37) = lhs - rhs;
lhs = y(66)-y(55)*y(68);
rhs = y(61)*(y(7)+y(6)+y(3)+y(2)+y(1));
residual(38) = lhs - rhs;
lhs = y(65)/(steady_state(40));
rhs = T(16)*T(17)*exp(x(it_, 7));
residual(39) = lhs - rhs;
residual(40) = (y(7)+y(1)+y(6))*y(62)+y(3)*y(63)+y(2)*y(64);
lhs = y(67);
rhs = y(42)*(y(7)+y(1)+y(6))+y(2)*y(48)+y(3)*y(52);
residual(41) = lhs - rhs;
lhs = y(68);
rhs = y(43)*(y(7)+y(1)+y(6))+params(14)*y(2)*y(49)+y(3)*y(53);
residual(42) = lhs - rhs;
lhs = y(66)*y(60);
rhs = params(16)*y(68);
residual(43) = lhs - rhs;
lhs = y(69);
rhs = y(66)-y(67);
residual(44) = lhs - rhs;
lhs = y(70);
rhs = y(1);
residual(45) = lhs - rhs;
lhs = y(71);
rhs = y(2);
residual(46) = lhs - rhs;
lhs = y(72);
rhs = y(3);
residual(47) = lhs - rhs;
lhs = y(73);
rhs = y(6);
residual(48) = lhs - rhs;
lhs = y(74);
rhs = y(7);
residual(49) = lhs - rhs;
lhs = y(82);
rhs = y(67);
residual(50) = lhs - rhs;
lhs = y(83);
rhs = y(68);
residual(51) = lhs - rhs;
lhs = y(84);
rhs = y(66);
residual(52) = lhs - rhs;
lhs = y(89);
rhs = y(56);
residual(53) = lhs - rhs;
lhs = y(90);
rhs = y(65);
residual(54) = lhs - rhs;
lhs = y(85);
rhs = y(26);
residual(55) = lhs - rhs;
lhs = y(86);
rhs = y(27);
residual(56) = lhs - rhs;
lhs = y(87);
rhs = y(28);
residual(57) = lhs - rhs;
lhs = y(88);
rhs = y(29);
residual(58) = lhs - rhs;
lhs = y(75);
rhs = params(27)*y(20)+x(it_, 1);
residual(59) = lhs - rhs;
lhs = y(76);
rhs = params(28)*y(21)+x(it_, 2);
residual(60) = lhs - rhs;
lhs = y(77);
rhs = params(29)*y(22)+x(it_, 3);
residual(61) = lhs - rhs;
lhs = y(78);
rhs = params(30)*y(23)+x(it_, 4);
residual(62) = lhs - rhs;
lhs = y(79);
rhs = params(31)*y(24)+x(it_, 5);
residual(63) = lhs - rhs;
lhs = y(80);
rhs = params(32)*y(25)+x(it_, 6);
residual(64) = lhs - rhs;
lhs = y(81);
rhs = (log(y(42))+params(11)*log(1-y(43)))*(y(7)+y(1)+y(6))+y(2)*(log(y(48))+params(11)*log(1-y(49)))+y(3)*(log(y(52))+params(11)*log(1-y(53)))-y(4)*(log(y(52))-params(11)*log(1-y(53)));
residual(65) = lhs - rhs;
if nargout >= 2,
  g1 = zeros(65, 99);

  %
  % Jacobian matrix
  %

T(18) = getPowerDeriv(exp(y(27)-y(2)),params(23),1);
T(19) = getPowerDeriv(T(12),params(21),1);
T(20) = (y(7)+y(1)+y(6))*getPowerDeriv(y(42),(-params(10)),1);
T(21) = getPowerDeriv(y(48),(-params(10)),1);
T(22) = getPowerDeriv(y(52),(-params(10)),1);
 g1 = zeros(65, 108);
g1(1,2)=(-(params(4)*1/params(6)*(1-(params(4)*(1+y(2)/params(6))>3*params(4)))));
g1(1,36)=1;
g1(2,36)=1;
g1(2,37)=1;
g1(3,1)=(-y(45));
g1(3,30)=1;
g1(3,45)=(-y(1));
g1(4,1)=(-1);
g1(4,26)=1;
g1(4,30)=1;
g1(5,2)=(-(1-(y(36)+y(37))));
g1(5,27)=1;
g1(5,30)=(-params(8));
g1(5,36)=y(2);
g1(5,37)=y(2);
g1(6,30)=(-(1-params(8)));
g1(6,6)=(-(1-(y(36)+y(37))));
g1(6,32)=1;
g1(6,36)=y(6);
g1(6,37)=y(6);
g1(7,2)=(-y(37));
g1(7,3)=(-1);
g1(7,28)=1;
g1(7,7)=(-params(9));
g1(7,37)=(-y(2));
g1(8,6)=(-(y(36)+y(37)));
g1(8,7)=(-(1-params(9)));
g1(8,33)=1;
g1(8,36)=(-y(6));
g1(8,37)=(-y(6));
g1(9,2)=(-y(36));
g1(9,4)=(-1);
g1(9,29)=1;
g1(9,36)=(-y(2));
g1(10,2)=y(36);
g1(10,5)=(-1);
g1(10,31)=1;
g1(10,36)=y(2);
g1(11,2)=(-(T(1)+params(1)*y(42)*params(13)*y(48)+params(2)*y(43)*params(15)*y(49)));
g1(11,6)=(-(params(1)*y(42)*y(42)*params(7)+params(2)*y(43)*params(7)*y(43)+params(7)*T(1)));
g1(11,42)=(-(params(1)*(y(2)*params(13)*y(48)+y(42)*y(6)*params(7))+params(1)*y(42)*y(6)*params(7)));
g1(11,43)=(-(params(2)*(y(2)*params(15)*y(49)+y(6)*params(7)*y(43))+y(6)*params(7)*params(2)*y(43)));
g1(11,45)=1;
g1(11,48)=(-(params(1)*y(42)*y(2)*params(13)));
g1(11,49)=(-(params(2)*y(43)*y(2)*params(15)));
g1(11,75)=(-((y(2)+y(6)*params(7))*params(3)*(-(getPowerDeriv(1-y(75),params(26),1)))));
g1(12,8)=(-((y(36)+y(37))/(1-y(45))));
g1(12,9)=(-((1-params(9))/(1-y(45))));
g1(12,35)=1;
g1(12,36)=(-(y(8)/(1-y(45))));
g1(12,37)=(-(y(8)/(1-y(45))));
g1(12,45)=(-(((1-params(9))*y(9)+(y(36)+y(37))*y(8))/((1-y(45))*(1-y(45)))));
g1(13,8)=(-((1-y(37)-y(36))/(1-y(45))));
g1(13,34)=1;
g1(13,36)=(-((-y(8))/(1-y(45))));
g1(13,37)=(-((-y(8))/(1-y(45))));
g1(13,45)=(-(y(8)*(1-y(37)-y(36))/((1-y(45))*(1-y(45)))+(y(45)*(1-params(8))+(1-params(8))*(1-y(45)))/((1-y(45))*(1-y(45)))));
g1(14,38)=1;
g1(14,39)=1;
g1(14,40)=1;
g1(15,8)=(-((-(y(45)*params(8)))/((y(9)+1+y(8))*(y(9)+1+y(8)))));
g1(15,9)=(-((-(y(45)*params(8)))/((y(9)+1+y(8))*(y(9)+1+y(8)))));
g1(15,39)=1;
g1(15,45)=(-(params(8)/(y(9)+1+y(8))));
g1(16,8)=(-((-(params(9)*y(9)))/((y(9)+1+y(8))*(y(9)+1+y(8)))));
g1(16,9)=(-((params(9)*(y(9)+1+y(8))-params(9)*y(9))/((y(9)+1+y(8))*(y(9)+1+y(8)))));
g1(16,40)=1;
g1(17,39)=(-(params(12)*(y(93)-y(91))));
g1(17,40)=(-(params(12)*(y(95)-y(91))));
g1(17,41)=1;
g1(17,91)=(-((1-y(39)-y(40))*params(12)));
g1(17,42)=(-(1/y(42)));
g1(17,43)=(-(params(11)*(-1)/(1-y(43))));
g1(17,93)=(-(y(39)*params(12)));
g1(17,95)=(-(y(40)*params(12)));
g1(18,1)=y(42)+y(62)-y(43)*y(55)-y(61)-y(14)*params(8)*y(11)/y(56);
g1(18,6)=y(42)+y(62)-y(43)*y(55)-y(61);
g1(18,7)=y(42)+y(62)-y(43)*y(55)-y(61);
g1(18,42)=y(7)+y(1)+y(6);
g1(18,43)=(y(7)+y(1)+y(6))*(-y(55));
g1(18,55)=(y(7)+y(1)+y(6))*(-y(43));
g1(18,56)=(-((-(y(11)*(y(15)-params(8)*(y(15)-y(1))+y(18)+(1-params(9))*y(19))*y(14)))/(y(56)*y(56))));
g1(18,61)=(-(y(7)+y(1)+y(6)));
g1(18,11)=(-((y(15)-params(8)*(y(15)-y(1))+y(18)+(1-params(9))*y(19))*y(14)/y(56)));
g1(18,62)=y(7)+y(1)+y(6);
g1(18,14)=(-(y(11)*(y(15)-params(8)*(y(15)-y(1))+y(18)+(1-params(9))*y(19))/y(56)));
g1(18,15)=(-(y(14)*(1-params(8))*y(11)/y(56)));
g1(18,18)=(-(y(11)*y(14)/y(56)));
g1(18,19)=(-(y(14)*(1-params(9))*y(11)/y(56)));
g1(19,2)=params(1)*y(46)*params(13)*y(48);
g1(19,6)=params(1)*y(46)*y(42)*params(7);
g1(19,42)=getPowerDeriv(y(42),(-1),1)+y(6)*params(7)*params(1)*y(46);
g1(19,44)=(-(1+y(75)));
g1(19,46)=params(1)*(y(2)*params(13)*y(48)+y(42)*y(6)*params(7));
g1(19,48)=y(2)*params(13)*params(1)*y(46);
g1(19,75)=(-y(44));
g1(20,8)=(-y(46));
g1(20,9)=(-y(46));
g1(20,91)=(-(params(8)*params(12)));
g1(20,46)=(-(y(9)+1+y(8)));
g1(20,93)=params(8)*params(12);
g1(21,2)=(-(params(2)*y(46)*params(15)*y(49)));
g1(21,6)=(-(params(2)*y(46)*params(7)*y(43)));
g1(21,43)=params(11)/((1-y(43))*(1-y(43)))-y(6)*params(7)*params(2)*y(46);
g1(21,44)=(-(y(55)*(1-y(78))));
g1(21,46)=(-(params(2)*(y(2)*params(15)*y(49)+y(6)*params(7)*y(43))));
g1(21,49)=(-(y(2)*params(15)*params(2)*y(46)));
g1(21,55)=(-(y(44)*(1-y(78))));
g1(21,78)=(-(y(44)*(-y(55))));
g1(22,39)=(-(y(65)*y(92)*(-params(12))/y(97)+y(65)*params(12)*y(94)/y(97)));
g1(22,40)=(-(y(65)*params(12)*y(96)/y(97)+y(65)*y(92)*(-params(12))/y(97)));
g1(22,44)=1;
g1(22,92)=(-((1-y(39)-y(40))*params(12)*y(65)/y(97)));
g1(22,94)=(-(y(65)*y(39)*params(12)/y(97)));
g1(22,96)=(-(y(65)*y(40)*params(12)/y(97)));
g1(22,97)=(-((-((1-y(39)-y(40))*params(12)*y(92)*y(65)))/(y(97)*y(97))+(-(y(65)*y(39)*params(12)*y(94)))/(y(97)*y(97))+(-(y(65)*y(40)*params(12)*y(96)))/(y(97)*y(97))));
g1(22,65)=(-((1-y(39)-y(40))*params(12)*y(92)/y(97)+y(39)*params(12)*y(94)/y(97)+y(40)*params(12)*y(96)/y(97)));
g1(23,36)=(-(params(12)*(params(19)-y(93))));
g1(23,37)=(-(params(12)*(y(95)-y(93))));
g1(23,47)=1;
g1(23,93)=(-((1-y(37)-y(36))*params(12)));
g1(23,48)=(-(1/y(48)));
g1(23,49)=(-(params(11)*(-1)/(1-y(49))));
g1(23,95)=(-(y(37)*params(12)));
g1(24,1)=(-(y(14)*(-(params(8)*y(11)))/y(56)));
g1(24,2)=y(48)+y(64)-y(49)*y(55)*params(14)-y(61);
g1(24,37)=(-(y(14)*(-(y(13)*y(16)))/y(56)));
g1(24,48)=y(2);
g1(24,49)=y(2)*(-(y(55)*params(14)));
g1(24,55)=y(2)*(-(y(49)*params(14)));
g1(24,56)=(-((-(y(14)*((1-y(37))*y(13)*y(16)+(y(15)-y(1))*params(8)*y(11))))/(y(56)*y(56))));
g1(24,61)=(-y(2));
g1(24,11)=(-(params(8)*(y(15)-y(1))*y(14)/y(56)));
g1(24,13)=(-(y(14)*(1-y(37))*y(16)/y(56)));
g1(24,64)=y(2);
g1(24,14)=(-(((1-y(37))*y(13)*y(16)+(y(15)-y(1))*params(8)*y(11))/y(56)));
g1(24,15)=(-(y(14)*params(8)*y(11)/y(56)));
g1(24,16)=(-(y(14)*(1-y(37))*y(13)/y(56)));
g1(25,48)=T(21);
g1(25,50)=(-(1+y(76)));
g1(25,76)=(-y(50));
g1(26,49)=1;
g1(26,50)=(-((-((-(params(11)*y(55)*params(14)))/(y(55)*params(14)*y(50)*y(55)*params(14)*y(50))))*(1-params(11)/(y(55)*params(14)*y(50))>0)));
g1(26,55)=(-((1-params(11)/(y(55)*params(14)*y(50))>0)*(-((-(params(11)*params(14)*y(50)))/(y(55)*params(14)*y(50)*y(55)*params(14)*y(50))))));
g1(27,36)=(-((y(65)*y(94)*(1-y(37)-y(36))*params(12)/y(97)+(1-y(36))*y(65)*y(94)*(-params(12))/y(97))/((1-y(36))*(1-y(36)))+y(65)*y(96)*y(37)*params(12)/y(97)/((1-y(36))*(1-y(36)))));
g1(27,37)=(-(y(65)*y(94)*(-params(12))/y(97)/(1-y(36))+y(65)*params(12)*y(96)/y(97)/(1-y(36))));
g1(27,50)=1;
g1(27,94)=(-(y(65)*(1-y(37)-y(36))*params(12)/y(97)/(1-y(36))));
g1(27,96)=(-(y(65)*y(37)*params(12)/y(97)/(1-y(36))));
g1(27,97)=(-((-(y(65)*y(94)*(1-y(37)-y(36))*params(12)))/(y(97)*y(97))/(1-y(36))+(-(y(65)*y(96)*y(37)*params(12)))/(y(97)*y(97))/(1-y(36))));
g1(27,65)=(-(y(94)*(1-y(37)-y(36))*params(12)/y(97)/(1-y(36))+y(96)*y(37)*params(12)/y(97)/(1-y(36))));
g1(28,51)=1;
g1(28,95)=(-params(12));
g1(28,52)=(-(1/y(52)));
g1(28,53)=(-(params(11)*(-1)/(1-y(53))));
g1(29,3)=y(52)+y(63)-y(55)*y(53)-y(61);
g1(29,37)=(-(y(14)*y(13)*y(16)/y(56)));
g1(29,52)=y(3);
g1(29,53)=y(3)*(-y(55));
g1(29,55)=y(3)*(-y(53));
g1(29,56)=(-((-(y(14)*(y(12)*y(17)+y(16)*y(37)*y(13)+y(19)*params(9)*y(11))))/(y(56)*y(56))));
g1(29,61)=(-y(3));
g1(29,11)=(-(y(14)*params(9)*y(19)/y(56)));
g1(29,12)=(-(y(14)*y(17)/y(56)));
g1(29,63)=y(3);
g1(29,13)=(-(y(14)*y(37)*y(16)/y(56)));
g1(29,14)=(-((y(12)*y(17)+y(16)*y(37)*y(13)+y(19)*params(9)*y(11))/y(56)));
g1(29,16)=(-(y(14)*y(37)*y(13)/y(56)));
g1(29,17)=(-(y(14)*y(12)/y(56)));
g1(29,19)=(-(y(14)*params(9)*y(11)/y(56)));
g1(30,52)=T(22);
g1(30,54)=(-(1+y(77)));
g1(30,77)=(-y(54));
g1(31,53)=params(11)/((1-y(53))*(1-y(53)));
g1(31,54)=(-(y(55)*(1-y(80))));
g1(31,55)=(-(y(54)*(1-y(80))));
g1(31,80)=(-(y(54)*(-y(55))));
g1(32,54)=1;
g1(32,96)=(-(params(12)*y(65)/y(97)));
g1(32,97)=(-((-(y(65)*params(12)*y(96)))/(y(97)*y(97))));
g1(32,65)=(-(params(12)*y(96)/y(97)));
g1(33,57)=1;
g1(33,58)=(-(1/y(59)));
g1(33,59)=(-((-y(58))/(y(59)*y(59))));
g1(34,1)=(-(T(4)*T(5)));
g1(34,2)=(-(T(2)*T(4)));
g1(34,3)=(-(T(3)*T(4)));
g1(34,6)=(-(T(4)*T(5)));
g1(34,7)=(-(T(4)*T(5)));
g1(34,42)=(-(T(4)*T(20)));
g1(34,48)=(-(T(4)*y(2)*T(21)));
g1(34,52)=(-(T(4)*y(3)*T(22)));
g1(34,55)=(-(T(6)*y(66)*1/params(16)));
g1(34,97)=(-(y(98)*params(12)*params(18)*getPowerDeriv(y(97),params(17),1)/(1-params(12)*params(18))));
g1(34,58)=1/(1-params(12)*params(18));
g1(34,98)=(-(T(7)/(1-params(12)*params(18))));
g1(34,66)=(-(y(55)/params(16)*T(6)));
g1(35,1)=(-(y(66)*T(5)));
g1(35,2)=(-(T(2)*y(66)));
g1(35,3)=(-(T(3)*y(66)));
g1(35,6)=(-(y(66)*T(5)));
g1(35,7)=(-(y(66)*T(5)));
g1(35,42)=(-(y(66)*T(20)));
g1(35,48)=(-(y(66)*y(2)*T(21)));
g1(35,52)=(-(y(66)*y(3)*T(22)));
g1(35,97)=(-(y(99)*params(12)*params(18)*getPowerDeriv(y(97),params(17)-1,1)/(1-params(12)*params(18))));
g1(35,59)=1/(1-params(12)*params(18));
g1(35,99)=(-(T(8)/(1-params(12)*params(18))));
g1(35,66)=(-T(6));
g1(36,56)=(-(params(18)*getPowerDeriv(y(56),params(17)-1,1)));
g1(36,57)=(-((1-params(18))*getPowerDeriv(y(57),1-params(17),1)));
g1(37,56)=(-(y(10)*params(18)*getPowerDeriv(y(56),params(17),1)));
g1(37,57)=(-((1-params(18))*getPowerDeriv(y(57),(-params(17)),1)));
g1(37,10)=(-T(9));
g1(37,60)=1;
g1(38,1)=(-y(61));
g1(38,2)=(-y(61));
g1(38,3)=(-y(61));
g1(38,6)=(-y(61));
g1(38,7)=(-y(61));
g1(38,55)=(-y(68));
g1(38,61)=(-(y(7)+y(6)+y(3)+y(2)+y(1)));
g1(38,66)=1;
g1(38,68)=(-y(55));
g1(39,2)=(-(exp(x(it_, 7))*T(16)*(-exp(y(27)-y(2)))*T(18)));
g1(39,27)=(-(exp(x(it_, 7))*(T(17)*T(14)*exp(y(27))*getPowerDeriv(exp(y(27)),params(22),1)+T(16)*exp(y(27)-y(2))*T(18))));
g1(39,31)=(-(exp(x(it_, 7))*T(17)*T(15)*T(10)*(-y(66))/(y(31)*y(31))/T(11)*T(19)));
g1(39,56)=(-(exp(x(it_, 7))*T(17)*T(15)*T(13)*1/(steady_state(31))*getPowerDeriv(y(56)/(steady_state(31)),params(20),1)));
g1(39,65)=1/(steady_state(40));
g1(39,66)=(-(exp(x(it_, 7))*T(17)*T(15)*T(10)*T(19)*1/y(31)/T(11)));
g1(39,106)=(-(T(16)*T(17)*exp(x(it_, 7))));
g1(39,107)=(-(exp(x(it_, 7))*T(17)*T(15)*T(10)*T(19)*(-(y(66)/y(31)))/(T(11)*T(11))));
g1(40,1)=y(62);
g1(40,2)=y(64);
g1(40,3)=y(63);
g1(40,6)=y(62);
g1(40,7)=y(62);
g1(40,62)=y(7)+y(1)+y(6);
g1(40,63)=y(3);
g1(40,64)=y(2);
g1(41,1)=(-y(42));
g1(41,2)=(-y(48));
g1(41,3)=(-y(52));
g1(41,6)=(-y(42));
g1(41,7)=(-y(42));
g1(41,42)=(-(y(7)+y(1)+y(6)));
g1(41,48)=(-y(2));
g1(41,52)=(-y(3));
g1(41,67)=1;
g1(42,1)=(-y(43));
g1(42,2)=(-(y(49)*params(14)));
g1(42,3)=(-y(53));
g1(42,6)=(-y(43));
g1(42,7)=(-y(43));
g1(42,43)=(-(y(7)+y(1)+y(6)));
g1(42,49)=(-(y(2)*params(14)));
g1(42,53)=(-y(3));
g1(42,68)=1;
g1(43,60)=y(66);
g1(43,66)=y(60);
g1(43,68)=(-params(16));
g1(44,66)=(-1);
g1(44,67)=1;
g1(44,69)=1;
g1(45,1)=(-1);
g1(45,70)=1;
g1(46,2)=(-1);
g1(46,71)=1;
g1(47,3)=(-1);
g1(47,72)=1;
g1(48,6)=(-1);
g1(48,73)=1;
g1(49,7)=(-1);
g1(49,74)=1;
g1(50,67)=(-1);
g1(50,82)=1;
g1(51,68)=(-1);
g1(51,83)=1;
g1(52,66)=(-1);
g1(52,84)=1;
g1(53,56)=(-1);
g1(53,89)=1;
g1(54,65)=(-1);
g1(54,90)=1;
g1(55,26)=(-1);
g1(55,85)=1;
g1(56,27)=(-1);
g1(56,86)=1;
g1(57,28)=(-1);
g1(57,87)=1;
g1(58,29)=(-1);
g1(58,88)=1;
g1(59,20)=(-params(27));
g1(59,75)=1;
%g1(59,100)=(-1);
g1(60,21)=(-params(28));
g1(60,76)=1;
%g1(60,101)=(-1);
g1(61,22)=(-params(29));
g1(61,77)=1;
%g1(61,102)=(-1);
g1(62,23)=(-params(30));
g1(62,78)=1;
%g1(62,103)=(-1);
g1(63,24)=(-params(31));
g1(63,79)=1;
%g1(63,104)=(-1);
g1(64,25)=(-params(32));
g1(64,80)=1;
%g1(64,105)=(-1);
g1(65,1)=(-(log(y(42))+params(11)*log(1-y(43))));
g1(65,2)=(-(log(y(48))+params(11)*log(1-y(49))));
g1(65,3)=(-(log(y(52))+params(11)*log(1-y(53))));
g1(65,4)=log(y(52))-params(11)*log(1-y(53));
g1(65,6)=(-(log(y(42))+params(11)*log(1-y(43))));
g1(65,7)=(-(log(y(42))+params(11)*log(1-y(43))));
g1(65,42)=(-((y(7)+y(1)+y(6))*1/y(42)));
g1(65,43)=(-((y(7)+y(1)+y(6))*params(11)*(-1)/(1-y(43))));
g1(65,48)=(-(y(2)*1/y(48)));
g1(65,49)=(-(y(2)*params(11)*(-1)/(1-y(49))));
g1(65,52)=(-(y(3)*1/y(52)-y(4)*1/y(52)));
g1(65,53)=(-(y(3)*params(11)*(-1)/(1-y(53))-y(4)*(-(params(11)*(-1)/(1-y(53))))));
g1(65,81)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],56,9801);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],56,970299);
end
end
end
end
