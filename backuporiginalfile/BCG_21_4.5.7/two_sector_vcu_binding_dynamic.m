function [residual, g1, g2, g3] = two_sector_vcu_binding_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(132, 1);
T20 = params(2)*params(6)/(y(248)-y(121)*params(2));
T41 = params(7)/2;
T43 = T41*(y(249)-y(122))^2;
T44 = y(122)^2;
T46 = params(7)*(y(249)-y(122))/y(122)+T43/T44;
T97 = params(3)/(1+params(3));
T98 = (params(9)*y(116)/y(119))^T97;
T116 = (1-params(9))^T97;
T117 = 1/(1+params(3));
T120 = params(9)^T97;
T123 = T116*y(118)^T117+T120*y(119)^T117;
T126 = (y(3)*y(135))^params(1);
T129 = (y(133)*y(134))^(1-params(1));
lhs =y(125);
rhs =1/(y(121)-params(2)*y(1))-T20*y(254)/y(136);
residual(1)= lhs-rhs;
lhs =y(125)*(1+params(7)*(y(122)-y(2))/y(2))-params(6)*y(250)*T46;
rhs =y(126)+y(127);
residual(2)= lhs-rhs;
residual(3) = y(250)*y(252)*y(253)-y(126)+params(6)*(1-params(8))*y(251);
lhs =y(124);
rhs =y(122)+(1-params(8))*y(3);
residual(4)= lhs-rhs;
lhs =y(3)*y(131);
rhs =params(11)*y(135)^params(10);
residual(5)= lhs-rhs;
lhs =y(129);
rhs =params(12)*y(128);
residual(6)= lhs-rhs;
lhs =y(118);
rhs =params(12)*(y(132)-params(13));
residual(7)= lhs-rhs;
lhs =y(131);
rhs =y(119)*params(1)*T98/(y(3)*y(135));
residual(8)= lhs-rhs;
lhs =y(130);
rhs =y(119)*T98*(1-params(1))/y(133);
residual(9)= lhs-rhs;
lhs =y(128);
rhs =(y(116)*(1-params(9))/y(118))^T97;
residual(10)= lhs-rhs;
lhs =y(116);
rhs =T123^(1+params(3));
residual(11)= lhs-rhs;
lhs =y(119);
rhs =T126*T129;
residual(12)= lhs-rhs;
lhs =y(116);
rhs =y(122)+y(120)+params(11)*y(135)^(1+params(10))/(1+params(10));
residual(13)= lhs-rhs;
lhs =y(132)-params(15);
rhs =params(4)*(y(4)-params(15))+y(137)+x(it_, 1);
residual(14)= lhs-rhs;
lhs =y(133)-params(16);
rhs =params(4)*(y(5)-params(16))+y(174)+x(it_, 2);
residual(15)= lhs-rhs;
residual(16) = y(122);
lhs =y(120);
rhs =y(121)*y(136);
residual(17)= lhs-rhs;
lhs =y(136)-1;
rhs =params(4)*(y(7)-1)+y(211);
residual(18)= lhs-rhs;
lhs =log(y(134));
rhs =params(5)*log(y(6))+x(it_, 3);
residual(19)= lhs-rhs;
lhs =y(122);
rhs =y(136)*y(123);
residual(20)= lhs-rhs;
lhs =y(116);
rhs =y(136)*y(117);
residual(21)= lhs-rhs;
lhs =y(137);
rhs =y(8);
residual(22)= lhs-rhs;
lhs =y(138);
rhs =y(9);
residual(23)= lhs-rhs;
lhs =y(139);
rhs =y(10);
residual(24)= lhs-rhs;
lhs =y(140);
rhs =y(11);
residual(25)= lhs-rhs;
lhs =y(141);
rhs =y(12);
residual(26)= lhs-rhs;
lhs =y(142);
rhs =y(13);
residual(27)= lhs-rhs;
lhs =y(143);
rhs =y(14);
residual(28)= lhs-rhs;
lhs =y(144);
rhs =y(15);
residual(29)= lhs-rhs;
lhs =y(145);
rhs =y(16);
residual(30)= lhs-rhs;
lhs =y(146);
rhs =y(17);
residual(31)= lhs-rhs;
lhs =y(147);
rhs =y(18);
residual(32)= lhs-rhs;
lhs =y(148);
rhs =y(19);
residual(33)= lhs-rhs;
lhs =y(149);
rhs =y(20);
residual(34)= lhs-rhs;
lhs =y(150);
rhs =y(21);
residual(35)= lhs-rhs;
lhs =y(151);
rhs =y(22);
residual(36)= lhs-rhs;
lhs =y(152);
rhs =y(23);
residual(37)= lhs-rhs;
lhs =y(153);
rhs =y(24);
residual(38)= lhs-rhs;
lhs =y(154);
rhs =y(25);
residual(39)= lhs-rhs;
lhs =y(155);
rhs =y(26);
residual(40)= lhs-rhs;
lhs =y(156);
rhs =y(27);
residual(41)= lhs-rhs;
lhs =y(157);
rhs =y(28);
residual(42)= lhs-rhs;
lhs =y(158);
rhs =y(29);
residual(43)= lhs-rhs;
lhs =y(159);
rhs =y(30);
residual(44)= lhs-rhs;
lhs =y(160);
rhs =y(31);
residual(45)= lhs-rhs;
lhs =y(161);
rhs =y(32);
residual(46)= lhs-rhs;
lhs =y(162);
rhs =y(33);
residual(47)= lhs-rhs;
lhs =y(163);
rhs =y(34);
residual(48)= lhs-rhs;
lhs =y(164);
rhs =y(35);
residual(49)= lhs-rhs;
lhs =y(165);
rhs =y(36);
residual(50)= lhs-rhs;
lhs =y(166);
rhs =y(37);
residual(51)= lhs-rhs;
lhs =y(167);
rhs =y(38);
residual(52)= lhs-rhs;
lhs =y(168);
rhs =y(39);
residual(53)= lhs-rhs;
lhs =y(169);
rhs =y(40);
residual(54)= lhs-rhs;
lhs =y(170);
rhs =y(41);
residual(55)= lhs-rhs;
lhs =y(171);
rhs =y(42);
residual(56)= lhs-rhs;
lhs =y(172);
rhs =y(43);
residual(57)= lhs-rhs;
lhs =y(173);
rhs =x(it_, 4);
residual(58)= lhs-rhs;
lhs =y(174);
rhs =y(44);
residual(59)= lhs-rhs;
lhs =y(175);
rhs =y(45);
residual(60)= lhs-rhs;
lhs =y(176);
rhs =y(46);
residual(61)= lhs-rhs;
lhs =y(177);
rhs =y(47);
residual(62)= lhs-rhs;
lhs =y(178);
rhs =y(48);
residual(63)= lhs-rhs;
lhs =y(179);
rhs =y(49);
residual(64)= lhs-rhs;
lhs =y(180);
rhs =y(50);
residual(65)= lhs-rhs;
lhs =y(181);
rhs =y(51);
residual(66)= lhs-rhs;
lhs =y(182);
rhs =y(52);
residual(67)= lhs-rhs;
lhs =y(183);
rhs =y(53);
residual(68)= lhs-rhs;
lhs =y(184);
rhs =y(54);
residual(69)= lhs-rhs;
lhs =y(185);
rhs =y(55);
residual(70)= lhs-rhs;
lhs =y(186);
rhs =y(56);
residual(71)= lhs-rhs;
lhs =y(187);
rhs =y(57);
residual(72)= lhs-rhs;
lhs =y(188);
rhs =y(58);
residual(73)= lhs-rhs;
lhs =y(189);
rhs =y(59);
residual(74)= lhs-rhs;
lhs =y(190);
rhs =y(60);
residual(75)= lhs-rhs;
lhs =y(191);
rhs =y(61);
residual(76)= lhs-rhs;
lhs =y(192);
rhs =y(62);
residual(77)= lhs-rhs;
lhs =y(193);
rhs =y(63);
residual(78)= lhs-rhs;
lhs =y(194);
rhs =y(64);
residual(79)= lhs-rhs;
lhs =y(195);
rhs =y(65);
residual(80)= lhs-rhs;
lhs =y(196);
rhs =y(66);
residual(81)= lhs-rhs;
lhs =y(197);
rhs =y(67);
residual(82)= lhs-rhs;
lhs =y(198);
rhs =y(68);
residual(83)= lhs-rhs;
lhs =y(199);
rhs =y(69);
residual(84)= lhs-rhs;
lhs =y(200);
rhs =y(70);
residual(85)= lhs-rhs;
lhs =y(201);
rhs =y(71);
residual(86)= lhs-rhs;
lhs =y(202);
rhs =y(72);
residual(87)= lhs-rhs;
lhs =y(203);
rhs =y(73);
residual(88)= lhs-rhs;
lhs =y(204);
rhs =y(74);
residual(89)= lhs-rhs;
lhs =y(205);
rhs =y(75);
residual(90)= lhs-rhs;
lhs =y(206);
rhs =y(76);
residual(91)= lhs-rhs;
lhs =y(207);
rhs =y(77);
residual(92)= lhs-rhs;
lhs =y(208);
rhs =y(78);
residual(93)= lhs-rhs;
lhs =y(209);
rhs =y(79);
residual(94)= lhs-rhs;
lhs =y(210);
rhs =x(it_, 5);
residual(95)= lhs-rhs;
lhs =y(211);
rhs =y(80);
residual(96)= lhs-rhs;
lhs =y(212);
rhs =y(81);
residual(97)= lhs-rhs;
lhs =y(213);
rhs =y(82);
residual(98)= lhs-rhs;
lhs =y(214);
rhs =y(83);
residual(99)= lhs-rhs;
lhs =y(215);
rhs =y(84);
residual(100)= lhs-rhs;
lhs =y(216);
rhs =y(85);
residual(101)= lhs-rhs;
lhs =y(217);
rhs =y(86);
residual(102)= lhs-rhs;
lhs =y(218);
rhs =y(87);
residual(103)= lhs-rhs;
lhs =y(219);
rhs =y(88);
residual(104)= lhs-rhs;
lhs =y(220);
rhs =y(89);
residual(105)= lhs-rhs;
lhs =y(221);
rhs =y(90);
residual(106)= lhs-rhs;
lhs =y(222);
rhs =y(91);
residual(107)= lhs-rhs;
lhs =y(223);
rhs =y(92);
residual(108)= lhs-rhs;
lhs =y(224);
rhs =y(93);
residual(109)= lhs-rhs;
lhs =y(225);
rhs =y(94);
residual(110)= lhs-rhs;
lhs =y(226);
rhs =y(95);
residual(111)= lhs-rhs;
lhs =y(227);
rhs =y(96);
residual(112)= lhs-rhs;
lhs =y(228);
rhs =y(97);
residual(113)= lhs-rhs;
lhs =y(229);
rhs =y(98);
residual(114)= lhs-rhs;
lhs =y(230);
rhs =y(99);
residual(115)= lhs-rhs;
lhs =y(231);
rhs =y(100);
residual(116)= lhs-rhs;
lhs =y(232);
rhs =y(101);
residual(117)= lhs-rhs;
lhs =y(233);
rhs =y(102);
residual(118)= lhs-rhs;
lhs =y(234);
rhs =y(103);
residual(119)= lhs-rhs;
lhs =y(235);
rhs =y(104);
residual(120)= lhs-rhs;
lhs =y(236);
rhs =y(105);
residual(121)= lhs-rhs;
lhs =y(237);
rhs =y(106);
residual(122)= lhs-rhs;
lhs =y(238);
rhs =y(107);
residual(123)= lhs-rhs;
lhs =y(239);
rhs =y(108);
residual(124)= lhs-rhs;
lhs =y(240);
rhs =y(109);
residual(125)= lhs-rhs;
lhs =y(241);
rhs =y(110);
residual(126)= lhs-rhs;
lhs =y(242);
rhs =y(111);
residual(127)= lhs-rhs;
lhs =y(243);
rhs =y(112);
residual(128)= lhs-rhs;
lhs =y(244);
rhs =y(113);
residual(129)= lhs-rhs;
lhs =y(245);
rhs =y(114);
residual(130)= lhs-rhs;
lhs =y(246);
rhs =y(115);
residual(131)= lhs-rhs;
lhs =y(247);
rhs =x(it_, 6);
residual(132)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(132, 260);

  %
  % Jacobian matrix
  %

T3 = (-1);
T516 = getPowerDeriv(params(9)*y(116)/y(119),T97,1);
T519 = y(119)*params(1)*params(9)/y(119)*T516;
T527 = getPowerDeriv(y(116)*(1-params(9))/y(118),T97,1);
T533 = (-(y(116)*(1-params(9))))/(y(118)*y(118));
T537 = T116*getPowerDeriv(y(118),T117,1);
T538 = getPowerDeriv(T123,1+params(3),1);
T543 = (-(params(9)*y(116)))/(y(119)*y(119));
T545 = params(1)*T516*T543;
T556 = T120*getPowerDeriv(y(119),T117,1);
T568 = y(254)*(-(params(2)*params(6)*(-params(2))))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)));
T598 = (y(122)*(-params(7))-params(7)*(y(249)-y(122)))/(y(122)*y(122))+(T44*T41*(-(2*(y(249)-y(122))))-T43*2*y(122))/(T44*T44);
T604 = params(7)/y(122)+T41*2*(y(249)-y(122))/T44;
T613 = getPowerDeriv(y(3)*y(135),params(1),1);
T627 = getPowerDeriv(y(133)*y(134),1-params(1),1);
  g1(1,1)=(-(params(2)/((y(121)-params(2)*y(1))*(y(121)-params(2)*y(1)))));
  g1(1,121)=(-(T3/((y(121)-params(2)*y(1))*(y(121)-params(2)*y(1)))-T568/y(136)));
  g1(1,248)=y(254)*(-(params(2)*params(6)))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136);
  g1(1,125)=1;
  g1(1,136)=(-(T20*y(254)))/(y(136)*y(136));
  g1(1,254)=T20/y(136);
  g1(2,2)=y(125)*(y(2)*(-params(7))-params(7)*(y(122)-y(2)))/(y(2)*y(2));
  g1(2,122)=y(125)*params(7)/y(2)-params(6)*y(250)*T598;
  g1(2,249)=(-(params(6)*y(250)*T604));
  g1(2,125)=1+params(7)*(y(122)-y(2))/y(2);
  g1(2,250)=(-(params(6)*T46));
  g1(2,126)=T3;
  g1(2,127)=T3;
  g1(3,250)=y(252)*y(253);
  g1(3,126)=T3;
  g1(3,251)=params(6)*(1-params(8));
  g1(3,252)=y(250)*y(253);
  g1(3,253)=y(250)*y(252);
  g1(4,122)=T3;
  g1(4,3)=(-(1-params(8)));
  g1(4,124)=1;
  g1(5,3)=y(131);
  g1(5,131)=y(3);
  g1(5,135)=(-(params(11)*getPowerDeriv(y(135),params(10),1)));
  g1(6,128)=(-params(12));
  g1(6,129)=1;
  g1(7,118)=1;
  g1(7,132)=(-params(12));
  g1(8,116)=(-(T519/(y(3)*y(135))));
  g1(8,119)=(-((params(1)*T98+y(119)*T545)/(y(3)*y(135))));
  g1(8,3)=(-((-(y(135)*y(119)*params(1)*T98))/(y(3)*y(135)*y(3)*y(135))));
  g1(8,131)=1;
  g1(8,135)=(-((-(y(3)*y(119)*params(1)*T98))/(y(3)*y(135)*y(3)*y(135))));
  g1(9,116)=(-(y(119)*(1-params(1))*params(9)/y(119)*T516/y(133)));
  g1(9,119)=(-((T98*(1-params(1))+y(119)*(1-params(1))*T516*T543)/y(133)));
  g1(9,130)=1;
  g1(9,133)=(-((-(y(119)*T98*(1-params(1))))/(y(133)*y(133))));
  g1(10,116)=(-((1-params(9))/y(118)*T527));
  g1(10,118)=(-(T527*T533));
  g1(10,128)=1;
  g1(11,116)=1;
  g1(11,118)=(-(T537*T538));
  g1(11,119)=(-(T538*T556));
  g1(12,119)=1;
  g1(12,3)=(-(T129*y(135)*T613));
  g1(12,133)=(-(T126*y(134)*T627));
  g1(12,134)=(-(T126*y(133)*T627));
  g1(12,135)=(-(T129*y(3)*T613));
  g1(13,116)=1;
  g1(13,120)=T3;
  g1(13,122)=T3;
  g1(13,135)=(-(params(11)*getPowerDeriv(y(135),1+params(10),1)/(1+params(10))));
  g1(14,4)=(-params(4));
  g1(14,132)=1;
  g1(14,137)=T3;
  g1(14,255)=T3;
  g1(15,5)=(-params(4));
  g1(15,133)=1;
  g1(15,174)=T3;
  g1(15,256)=T3;
  g1(16,122)=1;
  g1(17,120)=1;
  g1(17,121)=(-y(136));
  g1(17,136)=(-y(121));
  g1(18,7)=(-params(4));
  g1(18,136)=1;
  g1(18,211)=T3;
  g1(19,6)=(-(params(5)*1/y(6)));
  g1(19,134)=1/y(134);
  g1(19,257)=T3;
  g1(20,122)=1;
  g1(20,123)=(-y(136));
  g1(20,136)=(-y(123));
  g1(21,116)=1;
  g1(21,117)=(-y(136));
  g1(21,136)=(-y(117));
  g1(22,137)=1;
  g1(22,8)=T3;
  g1(23,138)=1;
  g1(23,9)=T3;
  g1(24,139)=1;
  g1(24,10)=T3;
  g1(25,140)=1;
  g1(25,11)=T3;
  g1(26,141)=1;
  g1(26,12)=T3;
  g1(27,142)=1;
  g1(27,13)=T3;
  g1(28,143)=1;
  g1(28,14)=T3;
  g1(29,144)=1;
  g1(29,15)=T3;
  g1(30,145)=1;
  g1(30,16)=T3;
  g1(31,146)=1;
  g1(31,17)=T3;
  g1(32,147)=1;
  g1(32,18)=T3;
  g1(33,148)=1;
  g1(33,19)=T3;
  g1(34,149)=1;
  g1(34,20)=T3;
  g1(35,150)=1;
  g1(35,21)=T3;
  g1(36,151)=1;
  g1(36,22)=T3;
  g1(37,152)=1;
  g1(37,23)=T3;
  g1(38,153)=1;
  g1(38,24)=T3;
  g1(39,154)=1;
  g1(39,25)=T3;
  g1(40,155)=1;
  g1(40,26)=T3;
  g1(41,156)=1;
  g1(41,27)=T3;
  g1(42,157)=1;
  g1(42,28)=T3;
  g1(43,158)=1;
  g1(43,29)=T3;
  g1(44,159)=1;
  g1(44,30)=T3;
  g1(45,160)=1;
  g1(45,31)=T3;
  g1(46,161)=1;
  g1(46,32)=T3;
  g1(47,162)=1;
  g1(47,33)=T3;
  g1(48,163)=1;
  g1(48,34)=T3;
  g1(49,164)=1;
  g1(49,35)=T3;
  g1(50,165)=1;
  g1(50,36)=T3;
  g1(51,166)=1;
  g1(51,37)=T3;
  g1(52,167)=1;
  g1(52,38)=T3;
  g1(53,168)=1;
  g1(53,39)=T3;
  g1(54,169)=1;
  g1(54,40)=T3;
  g1(55,170)=1;
  g1(55,41)=T3;
  g1(56,171)=1;
  g1(56,42)=T3;
  g1(57,172)=1;
  g1(57,43)=T3;
  g1(58,173)=1;
  g1(58,258)=T3;
  g1(59,174)=1;
  g1(59,44)=T3;
  g1(60,175)=1;
  g1(60,45)=T3;
  g1(61,176)=1;
  g1(61,46)=T3;
  g1(62,177)=1;
  g1(62,47)=T3;
  g1(63,178)=1;
  g1(63,48)=T3;
  g1(64,179)=1;
  g1(64,49)=T3;
  g1(65,180)=1;
  g1(65,50)=T3;
  g1(66,181)=1;
  g1(66,51)=T3;
  g1(67,182)=1;
  g1(67,52)=T3;
  g1(68,183)=1;
  g1(68,53)=T3;
  g1(69,184)=1;
  g1(69,54)=T3;
  g1(70,185)=1;
  g1(70,55)=T3;
  g1(71,186)=1;
  g1(71,56)=T3;
  g1(72,187)=1;
  g1(72,57)=T3;
  g1(73,188)=1;
  g1(73,58)=T3;
  g1(74,189)=1;
  g1(74,59)=T3;
  g1(75,190)=1;
  g1(75,60)=T3;
  g1(76,191)=1;
  g1(76,61)=T3;
  g1(77,192)=1;
  g1(77,62)=T3;
  g1(78,193)=1;
  g1(78,63)=T3;
  g1(79,194)=1;
  g1(79,64)=T3;
  g1(80,195)=1;
  g1(80,65)=T3;
  g1(81,196)=1;
  g1(81,66)=T3;
  g1(82,197)=1;
  g1(82,67)=T3;
  g1(83,198)=1;
  g1(83,68)=T3;
  g1(84,199)=1;
  g1(84,69)=T3;
  g1(85,200)=1;
  g1(85,70)=T3;
  g1(86,201)=1;
  g1(86,71)=T3;
  g1(87,202)=1;
  g1(87,72)=T3;
  g1(88,203)=1;
  g1(88,73)=T3;
  g1(89,204)=1;
  g1(89,74)=T3;
  g1(90,205)=1;
  g1(90,75)=T3;
  g1(91,206)=1;
  g1(91,76)=T3;
  g1(92,207)=1;
  g1(92,77)=T3;
  g1(93,208)=1;
  g1(93,78)=T3;
  g1(94,209)=1;
  g1(94,79)=T3;
  g1(95,210)=1;
  g1(95,259)=T3;
  g1(96,211)=1;
  g1(96,80)=T3;
  g1(97,212)=1;
  g1(97,81)=T3;
  g1(98,213)=1;
  g1(98,82)=T3;
  g1(99,214)=1;
  g1(99,83)=T3;
  g1(100,215)=1;
  g1(100,84)=T3;
  g1(101,216)=1;
  g1(101,85)=T3;
  g1(102,217)=1;
  g1(102,86)=T3;
  g1(103,218)=1;
  g1(103,87)=T3;
  g1(104,219)=1;
  g1(104,88)=T3;
  g1(105,220)=1;
  g1(105,89)=T3;
  g1(106,221)=1;
  g1(106,90)=T3;
  g1(107,222)=1;
  g1(107,91)=T3;
  g1(108,223)=1;
  g1(108,92)=T3;
  g1(109,224)=1;
  g1(109,93)=T3;
  g1(110,225)=1;
  g1(110,94)=T3;
  g1(111,226)=1;
  g1(111,95)=T3;
  g1(112,227)=1;
  g1(112,96)=T3;
  g1(113,228)=1;
  g1(113,97)=T3;
  g1(114,229)=1;
  g1(114,98)=T3;
  g1(115,230)=1;
  g1(115,99)=T3;
  g1(116,231)=1;
  g1(116,100)=T3;
  g1(117,232)=1;
  g1(117,101)=T3;
  g1(118,233)=1;
  g1(118,102)=T3;
  g1(119,234)=1;
  g1(119,103)=T3;
  g1(120,235)=1;
  g1(120,104)=T3;
  g1(121,236)=1;
  g1(121,105)=T3;
  g1(122,237)=1;
  g1(122,106)=T3;
  g1(123,238)=1;
  g1(123,107)=T3;
  g1(124,239)=1;
  g1(124,108)=T3;
  g1(125,240)=1;
  g1(125,109)=T3;
  g1(126,241)=1;
  g1(126,110)=T3;
  g1(127,242)=1;
  g1(127,111)=T3;
  g1(128,243)=1;
  g1(128,112)=T3;
  g1(129,244)=1;
  g1(129,113)=T3;
  g1(130,245)=1;
  g1(130,114)=T3;
  g1(131,246)=1;
  g1(131,115)=T3;
  g1(132,247)=1;
  g1(132,260)=T3;

if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(100,3);
T758 = getPowerDeriv(params(9)*y(116)/y(119),T97,2);
T759 = params(9)/y(119)*T758;
T783 = T543*T543*T758+T516*(-((-(params(9)*y(116)))*(y(119)+y(119))))/(y(119)*y(119)*y(119)*y(119));
T852 = getPowerDeriv(y(116)*(1-params(9))/y(118),T97,2);
T853 = (1-params(9))/y(118)*T852;
T874 = getPowerDeriv(T123,1+params(3),2);
T875 = T537*T874;
T889 = getPowerDeriv(y(3)*y(135),params(1),2);
T890 = y(135)*T889;
T896 = getPowerDeriv(y(133)*y(134),1-params(1),2);
T897 = y(134)*T896;
  v2(1,1)=1;
  v2(1,2)=1;
  v2(1,3)=(-((-(params(2)*((y(121)-params(2)*y(1))*(-params(2))+(y(121)-params(2)*y(1))*(-params(2)))))/((y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1)))));
  v2(2,1)=1;
  v2(2,2)=31201;
  v2(2,3)=(-(((y(121)-params(2)*y(1))*(-params(2))+(y(121)-params(2)*y(1))*(-params(2)))/((y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1)))));
  v2(3,1)=1;
  v2(3,2)=121;
  v2(3,3)=  v2(2,3);
  v2(4,1)=1;
  v2(4,2)=31321;
  v2(4,3)=(-((y(121)-params(2)*y(1)+y(121)-params(2)*y(1))/((y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1))*(y(121)-params(2)*y(1)))-y(254)*(-((-(params(2)*params(6)*(-params(2))))*((y(248)-y(121)*params(2))*(-params(2))+(y(248)-y(121)*params(2))*(-params(2)))))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136)));
  v2(5,1)=1;
  v2(5,2)=64341;
  v2(5,3)=y(254)*(-((-(params(2)*params(6)))*((y(248)-y(121)*params(2))*(-params(2))+(y(248)-y(121)*params(2))*(-params(2)))))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136);
  v2(6,1)=1;
  v2(6,2)=31448;
  v2(6,3)=  v2(5,3);
  v2(7,1)=1;
  v2(7,2)=64468;
  v2(7,3)=y(254)*(-((-(params(2)*params(6)))*(y(248)-y(121)*params(2)+y(248)-y(121)*params(2))))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136);
  v2(8,1)=1;
  v2(8,2)=35221;
  v2(8,3)=(-T568)/(y(136)*y(136));
  v2(9,1)=1;
  v2(9,2)=31336;
  v2(9,3)=  v2(8,3);
  v2(10,1)=1;
  v2(10,2)=35348;
  v2(10,3)=(-(y(254)*(-(params(2)*params(6)))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))))/(y(136)*y(136));
  v2(11,1)=1;
  v2(11,2)=64356;
  v2(11,3)=  v2(10,3);
  v2(12,1)=1;
  v2(12,2)=35236;
  v2(12,3)=(-((-(T20*y(254)))*(y(136)+y(136))))/(y(136)*y(136)*y(136)*y(136));
  v2(13,1)=1;
  v2(13,2)=65901;
  v2(13,3)=(-(params(2)*params(6)*(-params(2))))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136);
  v2(14,1)=1;
  v2(14,2)=31454;
  v2(14,3)=  v2(13,3);
  v2(15,1)=1;
  v2(15,2)=66028;
  v2(15,3)=(-(params(2)*params(6)))/((y(248)-y(121)*params(2))*(y(248)-y(121)*params(2)))/y(136);
  v2(16,1)=1;
  v2(16,2)=64474;
  v2(16,3)=  v2(15,3);
  v2(17,1)=1;
  v2(17,2)=65916;
  v2(17,3)=(-T20)/(y(136)*y(136));
  v2(18,1)=1;
  v2(18,2)=35354;
  v2(18,3)=  v2(17,3);
  v2(19,1)=2;
  v2(19,2)=262;
  v2(19,3)=y(125)*(-((y(2)*(-params(7))-params(7)*(y(122)-y(2)))*(y(2)+y(2))))/(y(2)*y(2)*y(2)*y(2));
  v2(20,1)=2;
  v2(20,2)=31462;
  v2(20,3)=y(125)*(-params(7))/(y(2)*y(2));
  v2(21,1)=2;
  v2(21,2)=382;
  v2(21,3)=  v2(20,3);
  v2(22,1)=2;
  v2(22,2)=31582;
  v2(22,3)=(-(params(6)*y(250)*((-((y(122)*(-params(7))-params(7)*(y(249)-y(122)))*(y(122)+y(122))))/(y(122)*y(122)*y(122)*y(122))+(T44*T44*(T41*(-(2*(y(249)-y(122))))*2*y(122)+T44*2*T41-(T41*(-(2*(y(249)-y(122))))*2*y(122)+2*T43))-(T44*T41*(-(2*(y(249)-y(122))))-T43*2*y(122))*(T44*2*y(122)+T44*2*y(122)))/(T44*T44*T44*T44))));
  v2(23,1)=2;
  v2(23,2)=64602;
  v2(23,3)=(-(params(6)*y(250)*((-params(7))/(y(122)*y(122))+(T44*T41*(-2)-2*y(122)*T41*2*(y(249)-y(122)))/(T44*T44))));
  v2(24,1)=2;
  v2(24,2)=31709;
  v2(24,3)=  v2(23,3);
  v2(25,1)=2;
  v2(25,2)=64729;
  v2(25,3)=(-(params(6)*y(250)*2*T41/T44));
  v2(26,1)=2;
  v2(26,2)=32242;
  v2(26,3)=(y(2)*(-params(7))-params(7)*(y(122)-y(2)))/(y(2)*y(2));
  v2(27,1)=2;
  v2(27,2)=385;
  v2(27,3)=  v2(26,3);
  v2(28,1)=2;
  v2(28,2)=32362;
  v2(28,3)=params(7)/y(2);
  v2(29,1)=2;
  v2(29,2)=31585;
  v2(29,3)=  v2(28,3);
  v2(30,1)=2;
  v2(30,2)=64862;
  v2(30,3)=(-(params(6)*T598));
  v2(31,1)=2;
  v2(31,2)=31710;
  v2(31,3)=  v2(30,3);
  v2(32,1)=2;
  v2(32,2)=64989;
  v2(32,3)=(-(params(6)*T604));
  v2(33,1)=2;
  v2(33,2)=64730;
  v2(33,3)=  v2(32,3);
  v2(34,1)=3;
  v2(34,2)=65510;
  v2(34,3)=y(253);
  v2(35,1)=3;
  v2(35,2)=64992;
  v2(35,3)=  v2(34,3);
  v2(36,1)=3;
  v2(36,2)=65770;
  v2(36,3)=y(252);
  v2(37,1)=3;
  v2(37,2)=64993;
  v2(37,3)=  v2(36,3);
  v2(38,1)=3;
  v2(38,2)=65772;
  v2(38,3)=y(250);
  v2(39,1)=3;
  v2(39,2)=65513;
  v2(39,3)=  v2(38,3);
  v2(40,1)=5;
  v2(40,2)=33803;
  v2(40,3)=1;
  v2(41,1)=5;
  v2(41,2)=651;
  v2(41,3)=  v2(40,3);
  v2(42,1)=5;
  v2(42,2)=34975;
  v2(42,3)=(-(params(11)*getPowerDeriv(y(135),params(10),2)));
  v2(43,1)=8;
  v2(43,2)=30016;
  v2(43,3)=(-(y(119)*params(1)*params(9)/y(119)*T759/(y(3)*y(135))));
  v2(44,1)=8;
  v2(44,2)=30796;
  v2(44,3)=(-((params(1)*params(9)/y(119)*T516+y(119)*params(1)*(T543*T759+T516*(-params(9))/(y(119)*y(119))))/(y(3)*y(135))));
  v2(45,1)=8;
  v2(45,2)=30019;
  v2(45,3)=  v2(44,3);
  v2(46,1)=8;
  v2(46,2)=30799;
  v2(46,3)=(-((T545+T545+y(119)*params(1)*T783)/(y(3)*y(135))));
  v2(47,1)=8;
  v2(47,2)=636;
  v2(47,3)=(-((-(y(135)*T519))/(y(3)*y(135)*y(3)*y(135))));
  v2(48,1)=8;
  v2(48,2)=29903;
  v2(48,3)=  v2(47,3);
  v2(49,1)=8;
  v2(49,2)=639;
  v2(49,3)=(-((-(y(135)*(params(1)*T98+y(119)*T545)))/(y(3)*y(135)*y(3)*y(135))));
  v2(50,1)=8;
  v2(50,2)=30683;
  v2(50,3)=  v2(49,3);
  v2(51,1)=8;
  v2(51,2)=523;
  v2(51,3)=(-((-((-(y(135)*y(119)*params(1)*T98))*(y(135)*y(3)*y(135)+y(135)*y(3)*y(135))))/(y(3)*y(135)*y(3)*y(135)*y(3)*y(135)*y(3)*y(135))));
  v2(52,1)=8;
  v2(52,2)=34956;
  v2(52,3)=(-((-(y(3)*T519))/(y(3)*y(135)*y(3)*y(135))));
  v2(53,1)=8;
  v2(53,2)=30035;
  v2(53,3)=  v2(52,3);
  v2(54,1)=8;
  v2(54,2)=34959;
  v2(54,3)=(-((-(y(3)*(params(1)*T98+y(119)*T545)))/(y(3)*y(135)*y(3)*y(135))));
  v2(55,1)=8;
  v2(55,2)=30815;
  v2(55,3)=  v2(54,3);
  v2(56,1)=8;
  v2(56,2)=34843;
  v2(56,3)=(-((y(3)*y(135)*y(3)*y(135)*(-(y(119)*params(1)*T98))-(-(y(3)*y(119)*params(1)*T98))*(y(135)*y(3)*y(135)+y(135)*y(3)*y(135)))/(y(3)*y(135)*y(3)*y(135)*y(3)*y(135)*y(3)*y(135))));
  v2(57,1)=8;
  v2(57,2)=655;
  v2(57,3)=  v2(56,3);
  v2(58,1)=8;
  v2(58,2)=34975;
  v2(58,3)=(-((-((-(y(3)*y(119)*params(1)*T98))*(y(3)*y(3)*y(135)+y(3)*y(3)*y(135))))/(y(3)*y(135)*y(3)*y(135)*y(3)*y(135)*y(3)*y(135))));
  v2(59,1)=9;
  v2(59,2)=30016;
  v2(59,3)=(-(y(119)*(1-params(1))*params(9)/y(119)*T759/y(133)));
  v2(60,1)=9;
  v2(60,2)=30796;
  v2(60,3)=(-(((1-params(1))*params(9)/y(119)*T516+y(119)*(1-params(1))*(T543*T759+T516*(-params(9))/(y(119)*y(119))))/y(133)));
  v2(61,1)=9;
  v2(61,2)=30019;
  v2(61,3)=  v2(60,3);
  v2(62,1)=9;
  v2(62,2)=30799;
  v2(62,3)=(-(((1-params(1))*T516*T543+(1-params(1))*T516*T543+y(119)*(1-params(1))*T783)/y(133)));
  v2(63,1)=9;
  v2(63,2)=34436;
  v2(63,3)=(-((-(y(119)*(1-params(1))*params(9)/y(119)*T516))/(y(133)*y(133))));
  v2(64,1)=9;
  v2(64,2)=30033;
  v2(64,3)=  v2(63,3);
  v2(65,1)=9;
  v2(65,2)=34439;
  v2(65,3)=(-((-(T98*(1-params(1))+y(119)*(1-params(1))*T516*T543))/(y(133)*y(133))));
  v2(66,1)=9;
  v2(66,2)=30813;
  v2(66,3)=  v2(65,3);
  v2(67,1)=9;
  v2(67,2)=34453;
  v2(67,3)=(-((-((-(y(119)*T98*(1-params(1))))*(y(133)+y(133))))/(y(133)*y(133)*y(133)*y(133))));
  v2(68,1)=10;
  v2(68,2)=30016;
  v2(68,3)=(-((1-params(9))/y(118)*T853));
  v2(69,1)=10;
  v2(69,2)=30536;
  v2(69,3)=(-(T533*T853+T527*(-(1-params(9)))/(y(118)*y(118))));
  v2(70,1)=10;
  v2(70,2)=30018;
  v2(70,3)=  v2(69,3);
  v2(71,1)=10;
  v2(71,2)=30538;
  v2(71,3)=(-(T533*T533*T852+T527*(-((-(y(116)*(1-params(9))))*(y(118)+y(118))))/(y(118)*y(118)*y(118)*y(118))));
  v2(72,1)=11;
  v2(72,2)=30538;
  v2(72,3)=(-(T538*T116*getPowerDeriv(y(118),T117,2)+T537*T875));
  v2(73,1)=11;
  v2(73,2)=30798;
  v2(73,3)=(-(T556*T875));
  v2(74,1)=11;
  v2(74,2)=30539;
  v2(74,3)=  v2(73,3);
  v2(75,1)=11;
  v2(75,2)=30799;
  v2(75,3)=(-(T556*T556*T874+T538*T120*getPowerDeriv(y(119),T117,2)));
  v2(76,1)=12;
  v2(76,2)=523;
  v2(76,3)=(-(T129*y(135)*T890));
  v2(77,1)=12;
  v2(77,2)=34323;
  v2(77,3)=(-(y(135)*T613*y(134)*T627));
  v2(78,1)=12;
  v2(78,2)=653;
  v2(78,3)=  v2(77,3);
  v2(79,1)=12;
  v2(79,2)=34453;
  v2(79,3)=(-(T126*y(134)*T897));
  v2(80,1)=12;
  v2(80,2)=34583;
  v2(80,3)=(-(y(135)*T613*y(133)*T627));
  v2(81,1)=12;
  v2(81,2)=654;
  v2(81,3)=  v2(80,3);
  v2(82,1)=12;
  v2(82,2)=34713;
  v2(82,3)=(-(T126*(T627+y(133)*T897)));
  v2(83,1)=12;
  v2(83,2)=34454;
  v2(83,3)=  v2(82,3);
  v2(84,1)=12;
  v2(84,2)=34714;
  v2(84,3)=(-(T126*y(133)*y(133)*T896));
  v2(85,1)=12;
  v2(85,2)=34843;
  v2(85,3)=(-(T129*(T613+y(3)*T890)));
  v2(86,1)=12;
  v2(86,2)=655;
  v2(86,3)=  v2(85,3);
  v2(87,1)=12;
  v2(87,2)=34973;
  v2(87,3)=(-(y(134)*T627*y(3)*T613));
  v2(88,1)=12;
  v2(88,2)=34455;
  v2(88,3)=  v2(87,3);
  v2(89,1)=12;
  v2(89,2)=34974;
  v2(89,3)=(-(y(133)*T627*y(3)*T613));
  v2(90,1)=12;
  v2(90,2)=34715;
  v2(90,3)=  v2(89,3);
  v2(91,1)=12;
  v2(91,2)=34975;
  v2(91,3)=(-(T129*y(3)*y(3)*T889));
  v2(92,1)=13;
  v2(92,2)=34975;
  v2(92,3)=(-(params(11)*getPowerDeriv(y(135),1+params(10),2)/(1+params(10))));
  v2(93,1)=17;
  v2(93,2)=35221;
  v2(93,3)=T3;
  v2(94,1)=17;
  v2(94,2)=31336;
  v2(94,3)=  v2(93,3);
  v2(95,1)=19;
  v2(95,2)=1306;
  v2(95,3)=(-(params(5)*T3/(y(6)*y(6))));
  v2(96,1)=19;
  v2(96,2)=34714;
  v2(96,3)=T3/(y(134)*y(134));
  v2(97,1)=20;
  v2(97,2)=35223;
  v2(97,3)=T3;
  v2(98,1)=20;
  v2(98,2)=31856;
  v2(98,3)=  v2(97,3);
  v2(99,1)=21;
  v2(99,2)=35217;
  v2(99,3)=T3;
  v2(100,1)=21;
  v2(100,2)=30296;
  v2(100,3)=  v2(99,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),132,67600);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],132,17576000);
end
end
end
end