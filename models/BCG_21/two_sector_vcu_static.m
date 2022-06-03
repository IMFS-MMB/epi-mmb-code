function [residual, g1, g2, g3] = two_sector_vcu_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 132, 1);

%
% Model equations
%

T71 = params(3)/(1+params(3));
T72 = (params(9)*y(1)/y(4))^T71;
T90 = (1-params(9))^T71;
T91 = 1/(1+params(3));
T94 = params(9)^T71;
T97 = T90*y(3)^T91+T94*y(4)^T91;
T100 = (y(20)*y(9))^params(1);
T103 = (y(18)*y(19))^(1-params(1));
lhs =y(10);
rhs =1/(y(6)-y(6)*params(2))-params(2)*params(6)/(y(6)-y(6)*params(2))*y(21)/y(21);
residual(1)= lhs-rhs;
lhs =y(10);
rhs =y(11)+y(12);
residual(2)= lhs-rhs;
residual(3) = y(10)*y(16)*y(20)-y(11)+y(11)*params(6)*(1-params(8));
lhs =y(9);
rhs =y(7)+(1-params(8))*y(9);
residual(4)= lhs-rhs;
lhs =y(16)*y(9);
rhs =params(11)*y(20)^params(10);
residual(5)= lhs-rhs;
lhs =y(14);
rhs =params(12)*y(13);
residual(6)= lhs-rhs;
lhs =y(3);
rhs =params(12)*(y(17)-params(13));
residual(7)= lhs-rhs;
lhs =y(16);
rhs =y(4)*params(1)*T72/(y(20)*y(9));
residual(8)= lhs-rhs;
lhs =y(15);
rhs =y(4)*T72*(1-params(1))/y(18);
residual(9)= lhs-rhs;
lhs =y(13);
rhs =(y(1)*(1-params(9))/y(3))^T71;
residual(10)= lhs-rhs;
lhs =y(1);
rhs =T97^(1+params(3));
residual(11)= lhs-rhs;
lhs =y(4);
rhs =T100*T103;
residual(12)= lhs-rhs;
lhs =y(1);
rhs =y(7)+y(5)+params(11)*y(20)^(1+params(10))/(1+params(10));
residual(13)= lhs-rhs;
lhs =y(17)-params(15);
rhs =(y(17)-params(15))*params(4)+y(22)+x(1);
residual(14)= lhs-rhs;
lhs =y(18)-params(16);
rhs =params(4)*(y(18)-params(16))+y(59)+x(2);
residual(15)= lhs-rhs;
residual(16) = y(12);
lhs =y(5);
rhs =y(6)*y(21);
residual(17)= lhs-rhs;
lhs =y(21)-1;
rhs =params(4)*(y(21)-1)+y(96);
residual(18)= lhs-rhs;
lhs =log(y(19));
rhs =log(y(19))*params(5)+x(3);
residual(19)= lhs-rhs;
lhs =y(7);
rhs =y(21)*y(8);
residual(20)= lhs-rhs;
lhs =y(1);
rhs =y(21)*y(2);
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(23);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(24);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(25);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =y(26);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =y(27);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(28);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(29);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =y(30);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =y(31);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =y(32);
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(33);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =y(34);
residual(33)= lhs-rhs;
lhs =y(34);
rhs =y(35);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =y(36);
residual(35)= lhs-rhs;
lhs =y(36);
rhs =y(37);
residual(36)= lhs-rhs;
lhs =y(37);
rhs =y(38);
residual(37)= lhs-rhs;
lhs =y(38);
rhs =y(39);
residual(38)= lhs-rhs;
lhs =y(39);
rhs =y(40);
residual(39)= lhs-rhs;
lhs =y(40);
rhs =y(41);
residual(40)= lhs-rhs;
lhs =y(41);
rhs =y(42);
residual(41)= lhs-rhs;
lhs =y(42);
rhs =y(43);
residual(42)= lhs-rhs;
lhs =y(43);
rhs =y(44);
residual(43)= lhs-rhs;
lhs =y(44);
rhs =y(45);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =y(46);
residual(45)= lhs-rhs;
lhs =y(46);
rhs =y(47);
residual(46)= lhs-rhs;
lhs =y(47);
rhs =y(48);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(49);
residual(48)= lhs-rhs;
lhs =y(49);
rhs =y(50);
residual(49)= lhs-rhs;
lhs =y(50);
rhs =y(51);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =y(52);
residual(51)= lhs-rhs;
lhs =y(52);
rhs =y(53);
residual(52)= lhs-rhs;
lhs =y(53);
rhs =y(54);
residual(53)= lhs-rhs;
lhs =y(54);
rhs =y(55);
residual(54)= lhs-rhs;
lhs =y(55);
rhs =y(56);
residual(55)= lhs-rhs;
lhs =y(56);
rhs =y(57);
residual(56)= lhs-rhs;
lhs =y(57);
rhs =y(58);
residual(57)= lhs-rhs;
lhs =y(58);
rhs =x(4);
residual(58)= lhs-rhs;
lhs =y(59);
rhs =y(60);
residual(59)= lhs-rhs;
lhs =y(60);
rhs =y(61);
residual(60)= lhs-rhs;
lhs =y(61);
rhs =y(62);
residual(61)= lhs-rhs;
lhs =y(62);
rhs =y(63);
residual(62)= lhs-rhs;
lhs =y(63);
rhs =y(64);
residual(63)= lhs-rhs;
lhs =y(64);
rhs =y(65);
residual(64)= lhs-rhs;
lhs =y(65);
rhs =y(66);
residual(65)= lhs-rhs;
lhs =y(66);
rhs =y(67);
residual(66)= lhs-rhs;
lhs =y(67);
rhs =y(68);
residual(67)= lhs-rhs;
lhs =y(68);
rhs =y(69);
residual(68)= lhs-rhs;
lhs =y(69);
rhs =y(70);
residual(69)= lhs-rhs;
lhs =y(70);
rhs =y(71);
residual(70)= lhs-rhs;
lhs =y(71);
rhs =y(72);
residual(71)= lhs-rhs;
lhs =y(72);
rhs =y(73);
residual(72)= lhs-rhs;
lhs =y(73);
rhs =y(74);
residual(73)= lhs-rhs;
lhs =y(74);
rhs =y(75);
residual(74)= lhs-rhs;
lhs =y(75);
rhs =y(76);
residual(75)= lhs-rhs;
lhs =y(76);
rhs =y(77);
residual(76)= lhs-rhs;
lhs =y(77);
rhs =y(78);
residual(77)= lhs-rhs;
lhs =y(78);
rhs =y(79);
residual(78)= lhs-rhs;
lhs =y(79);
rhs =y(80);
residual(79)= lhs-rhs;
lhs =y(80);
rhs =y(81);
residual(80)= lhs-rhs;
lhs =y(81);
rhs =y(82);
residual(81)= lhs-rhs;
lhs =y(82);
rhs =y(83);
residual(82)= lhs-rhs;
lhs =y(83);
rhs =y(84);
residual(83)= lhs-rhs;
lhs =y(84);
rhs =y(85);
residual(84)= lhs-rhs;
lhs =y(85);
rhs =y(86);
residual(85)= lhs-rhs;
lhs =y(86);
rhs =y(87);
residual(86)= lhs-rhs;
lhs =y(87);
rhs =y(88);
residual(87)= lhs-rhs;
lhs =y(88);
rhs =y(89);
residual(88)= lhs-rhs;
lhs =y(89);
rhs =y(90);
residual(89)= lhs-rhs;
lhs =y(90);
rhs =y(91);
residual(90)= lhs-rhs;
lhs =y(91);
rhs =y(92);
residual(91)= lhs-rhs;
lhs =y(92);
rhs =y(93);
residual(92)= lhs-rhs;
lhs =y(93);
rhs =y(94);
residual(93)= lhs-rhs;
lhs =y(94);
rhs =y(95);
residual(94)= lhs-rhs;
lhs =y(95);
rhs =x(5);
residual(95)= lhs-rhs;
lhs =y(96);
rhs =y(97);
residual(96)= lhs-rhs;
lhs =y(97);
rhs =y(98);
residual(97)= lhs-rhs;
lhs =y(98);
rhs =y(99);
residual(98)= lhs-rhs;
lhs =y(99);
rhs =y(100);
residual(99)= lhs-rhs;
lhs =y(100);
rhs =y(101);
residual(100)= lhs-rhs;
lhs =y(101);
rhs =y(102);
residual(101)= lhs-rhs;
lhs =y(102);
rhs =y(103);
residual(102)= lhs-rhs;
lhs =y(103);
rhs =y(104);
residual(103)= lhs-rhs;
lhs =y(104);
rhs =y(105);
residual(104)= lhs-rhs;
lhs =y(105);
rhs =y(106);
residual(105)= lhs-rhs;
lhs =y(106);
rhs =y(107);
residual(106)= lhs-rhs;
lhs =y(107);
rhs =y(108);
residual(107)= lhs-rhs;
lhs =y(108);
rhs =y(109);
residual(108)= lhs-rhs;
lhs =y(109);
rhs =y(110);
residual(109)= lhs-rhs;
lhs =y(110);
rhs =y(111);
residual(110)= lhs-rhs;
lhs =y(111);
rhs =y(112);
residual(111)= lhs-rhs;
lhs =y(112);
rhs =y(113);
residual(112)= lhs-rhs;
lhs =y(113);
rhs =y(114);
residual(113)= lhs-rhs;
lhs =y(114);
rhs =y(115);
residual(114)= lhs-rhs;
lhs =y(115);
rhs =y(116);
residual(115)= lhs-rhs;
lhs =y(116);
rhs =y(117);
residual(116)= lhs-rhs;
lhs =y(117);
rhs =y(118);
residual(117)= lhs-rhs;
lhs =y(118);
rhs =y(119);
residual(118)= lhs-rhs;
lhs =y(119);
rhs =y(120);
residual(119)= lhs-rhs;
lhs =y(120);
rhs =y(121);
residual(120)= lhs-rhs;
lhs =y(121);
rhs =y(122);
residual(121)= lhs-rhs;
lhs =y(122);
rhs =y(123);
residual(122)= lhs-rhs;
lhs =y(123);
rhs =y(124);
residual(123)= lhs-rhs;
lhs =y(124);
rhs =y(125);
residual(124)= lhs-rhs;
lhs =y(125);
rhs =y(126);
residual(125)= lhs-rhs;
lhs =y(126);
rhs =y(127);
residual(126)= lhs-rhs;
lhs =y(127);
rhs =y(128);
residual(127)= lhs-rhs;
lhs =y(128);
rhs =y(129);
residual(128)= lhs-rhs;
lhs =y(129);
rhs =y(130);
residual(129)= lhs-rhs;
lhs =y(130);
rhs =y(131);
residual(130)= lhs-rhs;
lhs =y(131);
rhs =y(132);
residual(131)= lhs-rhs;
lhs =y(132);
rhs =x(6);
residual(132)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(132, 132);

  %
  % Jacobian matrix
  %

T3 = (-1);
T374 = getPowerDeriv(params(9)*y(1)/y(4),T71,1);
T385 = getPowerDeriv(y(1)*(1-params(9))/y(3),T71,1);
T396 = getPowerDeriv(T97,1+params(3),1);
T434 = getPowerDeriv(y(20)*y(9),params(1),1);
T447 = getPowerDeriv(y(18)*y(19),1-params(1),1);
  g1(1,6)=(-((-(1-params(2)))/((y(6)-y(6)*params(2))*(y(6)-y(6)*params(2)))-y(21)*(-(params(2)*params(6)*(1-params(2))))/((y(6)-y(6)*params(2))*(y(6)-y(6)*params(2)))/y(21)));
  g1(1,10)=1;
  g1(2,10)=1;
  g1(2,11)=T3;
  g1(2,12)=T3;
  g1(3,10)=y(16)*y(20);
  g1(3,11)=T3+params(6)*(1-params(8));
  g1(3,16)=y(10)*y(20);
  g1(3,20)=y(10)*y(16);
  g1(4,7)=T3;
  g1(4,9)=1-(1-params(8));
  g1(5,9)=y(16);
  g1(5,16)=y(9);
  g1(5,20)=(-(params(11)*getPowerDeriv(y(20),params(10),1)));
  g1(6,13)=(-params(12));
  g1(6,14)=1;
  g1(7,3)=1;
  g1(7,17)=(-params(12));
  g1(8,1)=(-(y(4)*params(1)*params(9)/y(4)*T374/(y(20)*y(9))));
  g1(8,4)=(-((params(1)*T72+y(4)*params(1)*T374*(-(params(9)*y(1)))/(y(4)*y(4)))/(y(20)*y(9))));
  g1(8,9)=(-((-(y(20)*y(4)*params(1)*T72))/(y(20)*y(9)*y(20)*y(9))));
  g1(8,16)=1;
  g1(8,20)=(-((-(y(9)*y(4)*params(1)*T72))/(y(20)*y(9)*y(20)*y(9))));
  g1(9,1)=(-(y(4)*(1-params(1))*params(9)/y(4)*T374/y(18)));
  g1(9,4)=(-((T72*(1-params(1))+y(4)*(1-params(1))*T374*(-(params(9)*y(1)))/(y(4)*y(4)))/y(18)));
  g1(9,15)=1;
  g1(9,18)=(-((-(y(4)*T72*(1-params(1))))/(y(18)*y(18))));
  g1(10,1)=(-((1-params(9))/y(3)*T385));
  g1(10,3)=(-(T385*(-(y(1)*(1-params(9))))/(y(3)*y(3))));
  g1(10,13)=1;
  g1(11,1)=1;
  g1(11,3)=(-(T90*getPowerDeriv(y(3),T91,1)*T396));
  g1(11,4)=(-(T396*T94*getPowerDeriv(y(4),T91,1)));
  g1(12,4)=1;
  g1(12,9)=(-(T103*y(20)*T434));
  g1(12,18)=(-(T100*y(19)*T447));
  g1(12,19)=(-(T100*y(18)*T447));
  g1(12,20)=(-(T103*y(9)*T434));
  g1(13,1)=1;
  g1(13,5)=T3;
  g1(13,7)=T3;
  g1(13,20)=(-(params(11)*getPowerDeriv(y(20),1+params(10),1)/(1+params(10))));
  g1(14,17)=1-params(4);
  g1(14,22)=T3;
  g1(15,18)=1-params(4);
  g1(15,59)=T3;
  g1(16,12)=1;
  g1(17,5)=1;
  g1(17,6)=(-y(21));
  g1(17,21)=(-y(6));
  g1(18,21)=1-params(4);
  g1(18,96)=T3;
  g1(19,19)=1/y(19)-params(5)*1/y(19);
  g1(20,7)=1;
  g1(20,8)=(-y(21));
  g1(20,21)=(-y(8));
  g1(21,1)=1;
  g1(21,2)=(-y(21));
  g1(21,21)=(-y(2));
  g1(22,22)=1;
  g1(22,23)=T3;
  g1(23,23)=1;
  g1(23,24)=T3;
  g1(24,24)=1;
  g1(24,25)=T3;
  g1(25,25)=1;
  g1(25,26)=T3;
  g1(26,26)=1;
  g1(26,27)=T3;
  g1(27,27)=1;
  g1(27,28)=T3;
  g1(28,28)=1;
  g1(28,29)=T3;
  g1(29,29)=1;
  g1(29,30)=T3;
  g1(30,30)=1;
  g1(30,31)=T3;
  g1(31,31)=1;
  g1(31,32)=T3;
  g1(32,32)=1;
  g1(32,33)=T3;
  g1(33,33)=1;
  g1(33,34)=T3;
  g1(34,34)=1;
  g1(34,35)=T3;
  g1(35,35)=1;
  g1(35,36)=T3;
  g1(36,36)=1;
  g1(36,37)=T3;
  g1(37,37)=1;
  g1(37,38)=T3;
  g1(38,38)=1;
  g1(38,39)=T3;
  g1(39,39)=1;
  g1(39,40)=T3;
  g1(40,40)=1;
  g1(40,41)=T3;
  g1(41,41)=1;
  g1(41,42)=T3;
  g1(42,42)=1;
  g1(42,43)=T3;
  g1(43,43)=1;
  g1(43,44)=T3;
  g1(44,44)=1;
  g1(44,45)=T3;
  g1(45,45)=1;
  g1(45,46)=T3;
  g1(46,46)=1;
  g1(46,47)=T3;
  g1(47,47)=1;
  g1(47,48)=T3;
  g1(48,48)=1;
  g1(48,49)=T3;
  g1(49,49)=1;
  g1(49,50)=T3;
  g1(50,50)=1;
  g1(50,51)=T3;
  g1(51,51)=1;
  g1(51,52)=T3;
  g1(52,52)=1;
  g1(52,53)=T3;
  g1(53,53)=1;
  g1(53,54)=T3;
  g1(54,54)=1;
  g1(54,55)=T3;
  g1(55,55)=1;
  g1(55,56)=T3;
  g1(56,56)=1;
  g1(56,57)=T3;
  g1(57,57)=1;
  g1(57,58)=T3;
  g1(58,58)=1;
  g1(59,59)=1;
  g1(59,60)=T3;
  g1(60,60)=1;
  g1(60,61)=T3;
  g1(61,61)=1;
  g1(61,62)=T3;
  g1(62,62)=1;
  g1(62,63)=T3;
  g1(63,63)=1;
  g1(63,64)=T3;
  g1(64,64)=1;
  g1(64,65)=T3;
  g1(65,65)=1;
  g1(65,66)=T3;
  g1(66,66)=1;
  g1(66,67)=T3;
  g1(67,67)=1;
  g1(67,68)=T3;
  g1(68,68)=1;
  g1(68,69)=T3;
  g1(69,69)=1;
  g1(69,70)=T3;
  g1(70,70)=1;
  g1(70,71)=T3;
  g1(71,71)=1;
  g1(71,72)=T3;
  g1(72,72)=1;
  g1(72,73)=T3;
  g1(73,73)=1;
  g1(73,74)=T3;
  g1(74,74)=1;
  g1(74,75)=T3;
  g1(75,75)=1;
  g1(75,76)=T3;
  g1(76,76)=1;
  g1(76,77)=T3;
  g1(77,77)=1;
  g1(77,78)=T3;
  g1(78,78)=1;
  g1(78,79)=T3;
  g1(79,79)=1;
  g1(79,80)=T3;
  g1(80,80)=1;
  g1(80,81)=T3;
  g1(81,81)=1;
  g1(81,82)=T3;
  g1(82,82)=1;
  g1(82,83)=T3;
  g1(83,83)=1;
  g1(83,84)=T3;
  g1(84,84)=1;
  g1(84,85)=T3;
  g1(85,85)=1;
  g1(85,86)=T3;
  g1(86,86)=1;
  g1(86,87)=T3;
  g1(87,87)=1;
  g1(87,88)=T3;
  g1(88,88)=1;
  g1(88,89)=T3;
  g1(89,89)=1;
  g1(89,90)=T3;
  g1(90,90)=1;
  g1(90,91)=T3;
  g1(91,91)=1;
  g1(91,92)=T3;
  g1(92,92)=1;
  g1(92,93)=T3;
  g1(93,93)=1;
  g1(93,94)=T3;
  g1(94,94)=1;
  g1(94,95)=T3;
  g1(95,95)=1;
  g1(96,96)=1;
  g1(96,97)=T3;
  g1(97,97)=1;
  g1(97,98)=T3;
  g1(98,98)=1;
  g1(98,99)=T3;
  g1(99,99)=1;
  g1(99,100)=T3;
  g1(100,100)=1;
  g1(100,101)=T3;
  g1(101,101)=1;
  g1(101,102)=T3;
  g1(102,102)=1;
  g1(102,103)=T3;
  g1(103,103)=1;
  g1(103,104)=T3;
  g1(104,104)=1;
  g1(104,105)=T3;
  g1(105,105)=1;
  g1(105,106)=T3;
  g1(106,106)=1;
  g1(106,107)=T3;
  g1(107,107)=1;
  g1(107,108)=T3;
  g1(108,108)=1;
  g1(108,109)=T3;
  g1(109,109)=1;
  g1(109,110)=T3;
  g1(110,110)=1;
  g1(110,111)=T3;
  g1(111,111)=1;
  g1(111,112)=T3;
  g1(112,112)=1;
  g1(112,113)=T3;
  g1(113,113)=1;
  g1(113,114)=T3;
  g1(114,114)=1;
  g1(114,115)=T3;
  g1(115,115)=1;
  g1(115,116)=T3;
  g1(116,116)=1;
  g1(116,117)=T3;
  g1(117,117)=1;
  g1(117,118)=T3;
  g1(118,118)=1;
  g1(118,119)=T3;
  g1(119,119)=1;
  g1(119,120)=T3;
  g1(120,120)=1;
  g1(120,121)=T3;
  g1(121,121)=1;
  g1(121,122)=T3;
  g1(122,122)=1;
  g1(122,123)=T3;
  g1(123,123)=1;
  g1(123,124)=T3;
  g1(124,124)=1;
  g1(124,125)=T3;
  g1(125,125)=1;
  g1(125,126)=T3;
  g1(126,126)=1;
  g1(126,127)=T3;
  g1(127,127)=1;
  g1(127,128)=T3;
  g1(128,128)=1;
  g1(128,129)=T3;
  g1(129,129)=1;
  g1(129,130)=T3;
  g1(130,130)=1;
  g1(130,131)=T3;
  g1(131,131)=1;
  g1(131,132)=T3;
  g1(132,132)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],132,17424);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],132,2299968);
end
end
end
end