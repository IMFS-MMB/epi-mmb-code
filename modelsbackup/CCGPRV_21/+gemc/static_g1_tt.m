function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 529);

T = gemc.static_resid_tt(T, y, x, params);

T(374) = (-1);
T(375) = y(221)*y(34)*y(221)*y(34);
T(376) = (-((y(287)-y(202))*y(221)))/T(375)-(-(params(157)*(T(22)*y(287)-T(21)*T(5)*y(202))*y(221)*T(24)))/(y(34)*y(221)*T(24)*y(34)*y(221)*T(24));
T(377) = getPowerDeriv(T(25),(-params(322)),1);
T(378) = (-(y(221)*T(27)))/T(375)-(-(y(221)*T(24)*T(28)))/(y(34)*y(221)*T(24)*y(34)*y(221)*T(24));
T(379) = getPowerDeriv(T(29),(-params(322)),1);
T(380) = getPowerDeriv(T(25)/T(32),1-params(322),1);
T(381) = getPowerDeriv(T(34),1-params(266),1);
T(382) = getPowerDeriv(y(288)/(y(221)*y(34)),(1-params(322))*params(266),1);
T(383) = (-(y(288)*y(221)))/T(375)*T(382);
T(384) = y(34)*y(221)*y(38)*y(34)*y(221)*y(38);
T(385) = y(34)*y(221)*T(24)*y(38)*y(34)*y(221)*T(24)*y(38);
T(386) = getPowerDeriv(T(37),params(323),1);
T(387) = getPowerDeriv(T(29)/T(42),1-params(322),1);
T(388) = getPowerDeriv(T(44),1-params(266),1);
T(389) = getPowerDeriv(T(46),params(324),1);
T(390) = getPowerDeriv(y(34)*y(221)*y(227),params(74)/(params(74)+params(76)-1),1);
T(391) = getPowerDeriv(T(52),params(9),1);
T(392) = getPowerDeriv(T(54),params(8),1);
T(393) = y(268)*params(8)*T(51)*y(221)*y(227)*T(390)*T(391)*T(392);
T(394) = getPowerDeriv(T(81),1-params(322),1);
T(395) = getPowerDeriv(T(83),1-params(266),1);
T(396) = getPowerDeriv(T(87),1-params(322),1);
T(397) = getPowerDeriv(T(89),1-params(266),1);
T(398) = getPowerDeriv(T(92),1-params(145),1);
T(399) = getPowerDeriv(y(226)*(y(223)-exp(y(154))*y(221)*y(34)*y(42)*params(118)),params(74),1);
T(400) = getPowerDeriv(T(109)/(y(221)*y(34)),1-params(74),1);
T(401) = (y(223)-exp(y(154))*y(221)*y(34)*params(118))*(y(223)-exp(y(154))*y(221)*y(34)*params(118));
T(402) = getPowerDeriv(y(42)*y(227)*(y(35)-params(118)),params(74),1);
T(403) = 4*(1+y(120))^3;
T(404) = getPowerDeriv(y(268)*params(285)*y(121)/((1-y(187))*(params(285)-1)),1-params(309),1);
T(405) = getPowerDeriv(T(109)*y(129),1-params(74),1);
T(406) = (-exp(y(182)))/(exp(y(182))*exp(y(182)))*getPowerDeriv(1/exp(y(182)),params(287)-1,1);
T(407) = T(67)*exp(params(58)-params(56)+params(147)*(1-params(6)))*(-(y(258)*T(2)*T(7)))/T(70);
T(408) = getPowerDeriv(T(72),params(144),1);
T(409) = getPowerDeriv(T(72),params(145),1);
T(410) = getPowerDeriv(T(65),1-params(144),1);
T(411) = getPowerDeriv(y(285)-(y(202)+(1-params(173))*y(204)),1-params(182),1);
T(412) = T(374)/(y(221)*y(34))-params(157)*(-T(22))/(y(34)*y(221)*T(24));
T(413) = (-((y(286))/(y(287))))/(y(221)*y(34))-params(158)*(-(T(21)*T(5)*(y(286))/(y(287))))/(y(34)*y(221)*T(24));
T(414) = (-(T(379)*T(413)/(y(221)*y(34))));
T(415) = (-(T(49)*T(388)*T(41)*T(387)*T(413)/T(42)));
T(416) = y(258)*T(76)*T(74)*params(142)*params(26)*T(310)/(1+y(120))*(T(2)*(y(205)+y(223))-(T(2)*y(223)+T(2)*y(205)))/((y(205)+y(223))*(y(205)+y(223)))/T(64);
T(417) = 1/T(109)-T(144)*(T(326)*exp(params(58)+params(149))*((1-params(268))*(1-params(269))*(T(107)-exp(params(58)+params(149))*T(146)*T(147))+(1-params(268))*params(269)*(1-T(107)*exp(params(58)+params(149)))+params(268)*(T(1)*exp(params(149))-exp(params(58)+params(149))))/y(293)-params(138)*(1-T(107)*exp(params(58)+params(149)))/T(109));
T(418) = getPowerDeriv(y(216),params(276)-1,1);
T(419) = getPowerDeriv(y(217),params(277)-1,1);
T(420) = getPowerDeriv(y(218),params(282)-1,1);
T(421) = getPowerDeriv(y(219),params(283)-1,1);
T(422) = (-((y(287)-y(202))*y(34)))/T(375)-(-(params(157)*(T(22)*y(287)-T(21)*T(5)*y(202))*y(34)*T(24)))/(y(34)*y(221)*T(24)*y(34)*y(221)*T(24));
T(423) = (-(y(34)*T(27)))/T(375)-(-(T(28)*y(34)*T(24)))/(y(34)*y(221)*T(24)*y(34)*y(221)*T(24));
T(424) = T(382)*(-(y(288)*y(34)))/T(375);
T(425) = 1/(y(34)*y(221)*y(38))-T(24)*params(159)/(y(34)*y(221)*T(24)*y(38));
T(426) = (-(1/T(190)));
T(427) = 1/T(11);
T(428) = T(67)*exp(params(58)-params(56)+params(147)*(1-params(6)))*(1-y(188))*T(2)*T(7)/T(70);
T(429) = 1/y(268);
T(430) = T(71)*getPowerDeriv(y(261),1-params(6),1)+T(67)*(-(exp(params(58)-params(56)+params(147)*(1-params(6)))*(1-y(188))*y(258)*T(2)*T(7)*T(68)*exp(params(147))*T(6)*getPowerDeriv(y(261)*exp(params(147))*T(6),1-params(6),1)))/(T(70)*T(70));
T(431) = (-T(11))/((y(261))*(y(261)));
T(432) = getPowerDeriv(y(274)/y(262),(-params(276)),1);
T(433) = getPowerDeriv(y(270)/y(262),(-params(276)),1);
T(434) = getPowerDeriv(y(274)/y(263),(-params(277)),1);
T(435) = getPowerDeriv(y(270)/y(263),(-params(277)),1);
T(436) = (y(289)-T(188)*T(6)*exp(params(148)))/T(190)-T(196)*T(188)*T(6)*exp(params(148))/T(190);
T(437) = getPowerDeriv(y(274)/y(264),(-params(283)),1);
T(438) = getPowerDeriv(y(270)/y(264),(-params(283)),1);
T(439) = (y(291)-T(107)*y(291)*T(6)*exp(params(149)))/T(190)-T(196)*T(107)*y(291)*T(6)*exp(params(149))/T(190);
T(440) = 1/y(267);
T(441) = (T(122)*T(4)/exp(params(149))/(T(4)*y(268))-T(319)*T(429))/(T(122)*T(122));
T(442) = (-(y(268)*y(321)/T(109)))/(y(267)*y(267));
T(443) = getPowerDeriv(y(274)/y(267),(-params(282)),1);
T(444) = getPowerDeriv(y(270)/y(267),(-params(282)),1);
T(445) = getPowerDeriv(y(268),params(6),1);
T(446) = T(67)*(-(exp(params(58)-params(56)+params(147)*(1-params(6)))*(1-y(188))*y(258)*T(2)*T(7)*T(69)*T(6)*getPowerDeriv(T(6)*y(268),params(6),1)))/(T(70)*T(70));
T(447) = getPowerDeriv(y(268)/y(277),(-params(284)),1);
T(448) = (-y(267))/(y(268)*y(268));
T(449) = T(128)*(-((params(278)-1)*params(368)/2*y(593)*y(594)))/(y(268)*y(268));
T(450) = (T(122)*(-(T(4)*T(318)))/(T(4)*y(268)*T(4)*y(268))-T(319)*T(448))/(T(122)*T(122));
T(451) = y(321)/T(109)/y(267);
T(452) = getPowerDeriv(y(270)/y(274),(-params(288)),1);
T(453) = getPowerDeriv(y(270)/y(276),(-params(287)),1);
T(454) = getPowerDeriv(y(594)/y(270),(-params(281)),1);
T(455) = T(226)*T(226);
T(456) = getPowerDeriv(y(273)/y(277),(-params(284)),1);
T(457) = getPowerDeriv(y(277)/y(274),(-params(288)),1);
T(458) = (y(274)*T(6)*y(594)/(T(6)*y(594))-y(594)/(T(6)*y(594))*T(6)*y(274))/(y(274)*y(274));
T(459) = T(126)*T(458)*2*T(127);
T(460) = getPowerDeriv(y(274)/y(276),(-params(287)),1);
T(461) = getPowerDeriv(y(597)*y(274),params(523),1);
T(462) = T(11)*T(11);
T(463) = getPowerDeriv(T(11),1-params(74),1);
T(464) = T(190)*T(190);
T(465) = T(193)*T(193);
T(466) = T(306)*(y(280)*y(44)*T(3)-y(44)*T(3)*y(280))/(y(280)*y(280));
T(467) = getPowerDeriv((y(285))/(y(287)),1-params(322),1);
T(468) = getPowerDeriv((y(285))/(y(286)),1-params(322),1);
T(469) = (1-(y(202)+y(204))/(y(287)))/(y(221)*y(34))-params(158)*(T(22)-T(21)*T(5)*(y(202)+y(204))/(y(287)))/(y(34)*y(221)*T(24));
T(470) = 1/(y(221)*y(34));
T(471) = T(470)-T(22)*params(157)/(y(34)*y(221)*T(24));
T(472) = (-((-((y(202)+y(204))*(y(286))))/((y(287))*(y(287)))))/(y(221)*y(34))-params(158)*(-(T(21)*T(5)*(-((y(202)+y(204))*(y(286))))/((y(287))*(y(287)))))/(y(34)*y(221)*T(24));
T(473) = T(382)*T(470);
T(474) = T(5)*1/exp(params(148))*y(263)*T(6)*exp(params(148));
T(475) = (y(263)-T(474))/T(190)-T(196)*T(474)/T(190);
T(476) = (-((T(107)*y(290)+T(107)*y(291)-T(107)*(y(290)+y(291)))/((T(107)*y(290)+T(107)*y(291))*(T(107)*y(290)+T(107)*y(291)))/((y(290)+y(291))/(T(107)*y(290)+T(107)*y(291)))));
T(477) = (y(264)-T(107)*y(264)*T(6)*exp(params(149)))/T(190)-T(196)*T(107)*y(264)*T(6)*exp(params(149))/T(190);
T(478) = y(268)*(-(T(107)*y(321)))/(T(109)*T(109))/y(267);
T(479) = T(107)*getPowerDeriv(y(294)*T(107),1-params(76),1);
T(480) = getPowerDeriv(y(295),1-params(74),1);
T(481) = y(268)*1/T(109)/y(267);
T(482) = 1/y(329)/T(119);
T(483) = ((y(328))-y(328))/((y(328))*(y(328)))*getPowerDeriv(T(50),1-params(9),1);
T(484) = y(268)*params(8)*T(392)*T(53)*T(483);
T(485) = getPowerDeriv(T(51)*T(57),params(8),1);
T(486) = y(268)*params(8)*T(485)*T(51)*getPowerDeriv(y(329),params(9),1);
T(487) = (-y(326))/(y(329)*y(329))/T(119);
T(488) = getPowerDeriv(T(364),1-params(393),1);
T(489) = getPowerDeriv(T(248)*y(427),1-params(373),1);
T(490) = 1/T(248)-T(250)*(T(355)*exp(params(58)+params(405))*((1-params(413))*(1-params(414))*(T(5)*T(242)-exp(params(58)+params(405))*T(147)*T(252))+(1-params(413))*params(414)*(1-T(5)*T(242)*exp(params(58)+params(405)))+params(413)*(T(1)*exp(params(405))-exp(params(58)+params(405))))/y(502)-params(402)*(1-T(5)*T(242)*exp(params(58)+params(405)))/T(248));
T(491) = (y(513)*y(530)+y(521)*y(516)*y(462))*(y(513)*y(530)+y(521)*y(516)*y(462));
T(492) = getPowerDeriv(T(238),(-params(507)),1);
T(493) = (-exp(y(499)))/(exp(y(499))*exp(y(499)))*getPowerDeriv(1/exp(y(499)),params(499)-1,1);
T(494) = 1/y(508);
T(495) = (T(241)*T(4)/exp(params(405))/(T(4)*y(508))-T(348)*T(494))/(T(241)*T(241));
T(496) = (-(y(508)*y(526)/T(248)))/(y(500)*y(500));
T(497) = getPowerDeriv(y(514)/y(500),(-params(498)),1);
T(498) = getPowerDeriv(y(513)/y(500),(-params(498)),1);
T(499) = y(508)*(-(T(5)*T(242)*y(526)))/(T(248)*T(248))/y(500);
T(500) = getPowerDeriv(y(503),params(498)-1,1);
T(501) = getPowerDeriv(y(504),params(498)-1,1);
T(502) = getPowerDeriv(y(506)/y(505),params(392),1);
T(503) = getPowerDeriv(T(259),1-params(507),1);
T(504) = getPowerDeriv(params(381)*y(505)*y(538),params(373),1);
T(505) = 1/y(505);
T(506) = getPowerDeriv(y(506)*params(381)*y(537),params(373),1);
T(507) = getPowerDeriv(T(255),params(393),1);
T(508) = 1/y(549);
T(509) = (T(241)*(-(T(4)*T(347)))/(T(4)*y(508)*T(4)*y(508))-T(348)*(-y(500))/(y(508)*y(508)))/(T(241)*T(241));
T(510) = y(526)/T(248)/y(500);
T(511) = T(275)*(-(y(605)*y(603)*params(520)*(params(500)-1)/2))/(y(508)*y(508));
T(512) = getPowerDeriv(T(297),1/(1-params(486)),1);
T(513) = getPowerDeriv(y(513)/y(549),(-params(498)),1);
T(514) = getPowerDeriv(y(513)/y(517),(-params(499)),1);
T(515) = getPowerDeriv(y(605)/y(513),(-params(503)),1);
T(516) = getPowerDeriv(y(587)*y(514),params(371),1);
T(517) = (y(514)*T(6)*y(605)/(T(6)*y(605))-y(605)/(T(6)*y(605))*T(6)*y(514))/(y(514)*y(514));
T(518) = T(273)*T(517)*2*T(274);
T(519) = getPowerDeriv(y(514)/y(549),(-params(498)),1);
T(520) = getPowerDeriv(y(514)/y(517),(-params(499)),1);
T(521) = (-((-y(520))/(y(518)*y(518))/(y(520)/y(518))));
T(522) = (-(1/y(518)/(y(520)/y(518))));
T(523) = y(508)*1/T(248)/y(500);
T(524) = T(343)*(y(547)*y(457)*T(234)-y(457)*T(234)*y(547))/(y(547)*y(547));
T(525) = T(505)-params(425)*T(236)/(T(24)*y(505));
T(526) = T(128)*(params(278)-1)*params(368)/2*y(594)/y(268);
T(527) = T(128)*(params(278)-1)*params(368)/2*y(593)/y(268);
T(528) = T(275)*y(605)*params(520)*(params(500)-1)/2/y(508);
T(529) = T(275)*y(603)*params(520)*(params(500)-1)/2/y(508);

end
