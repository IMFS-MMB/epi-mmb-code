function [nzij_pred, nzij_current, nzij_fwrd] = dynamic_g1_nz()
% Returns the coordinates of non-zero elements in the Jacobian, in column-major order, for each lead/lag (only for endogenous)
  nzij_pred = zeros(116, 2, 'int32');
  nzij_pred(1,1)=1; nzij_pred(1,2)=5;
  nzij_pred(2,1)=4; nzij_pred(2,2)=8;
  nzij_pred(3,1)=5; nzij_pred(3,2)=8;
  nzij_pred(4,1)=6; nzij_pred(4,2)=8;
  nzij_pred(5,1)=9; nzij_pred(5,2)=8;
  nzij_pred(6,1)=11; nzij_pred(6,2)=14;
  nzij_pred(7,1)=15; nzij_pred(7,2)=15;
  nzij_pred(8,1)=14; nzij_pred(8,2)=17;
  nzij_pred(9,1)=18; nzij_pred(9,2)=19;
  nzij_pred(10,1)=19; nzij_pred(10,2)=20;
  nzij_pred(11,1)=20; nzij_pred(11,2)=21;
  nzij_pred(12,1)=21; nzij_pred(12,2)=22;
  nzij_pred(13,1)=22; nzij_pred(13,2)=23;
  nzij_pred(14,1)=23; nzij_pred(14,2)=24;
  nzij_pred(15,1)=24; nzij_pred(15,2)=25;
  nzij_pred(16,1)=25; nzij_pred(16,2)=26;
  nzij_pred(17,1)=26; nzij_pred(17,2)=27;
  nzij_pred(18,1)=27; nzij_pred(18,2)=28;
  nzij_pred(19,1)=28; nzij_pred(19,2)=29;
  nzij_pred(20,1)=29; nzij_pred(20,2)=30;
  nzij_pred(21,1)=30; nzij_pred(21,2)=31;
  nzij_pred(22,1)=31; nzij_pred(22,2)=32;
  nzij_pred(23,1)=32; nzij_pred(23,2)=33;
  nzij_pred(24,1)=33; nzij_pred(24,2)=34;
  nzij_pred(25,1)=34; nzij_pred(25,2)=35;
  nzij_pred(26,1)=35; nzij_pred(26,2)=36;
  nzij_pred(27,1)=36; nzij_pred(27,2)=37;
  nzij_pred(28,1)=37; nzij_pred(28,2)=38;
  nzij_pred(29,1)=38; nzij_pred(29,2)=39;
  nzij_pred(30,1)=39; nzij_pred(30,2)=40;
  nzij_pred(31,1)=40; nzij_pred(31,2)=41;
  nzij_pred(32,1)=41; nzij_pred(32,2)=42;
  nzij_pred(33,1)=42; nzij_pred(33,2)=43;
  nzij_pred(34,1)=43; nzij_pred(34,2)=44;
  nzij_pred(35,1)=44; nzij_pred(35,2)=45;
  nzij_pred(36,1)=45; nzij_pred(36,2)=46;
  nzij_pred(37,1)=46; nzij_pred(37,2)=47;
  nzij_pred(38,1)=47; nzij_pred(38,2)=48;
  nzij_pred(39,1)=48; nzij_pred(39,2)=49;
  nzij_pred(40,1)=49; nzij_pred(40,2)=50;
  nzij_pred(41,1)=50; nzij_pred(41,2)=51;
  nzij_pred(42,1)=51; nzij_pred(42,2)=52;
  nzij_pred(43,1)=52; nzij_pred(43,2)=53;
  nzij_pred(44,1)=53; nzij_pred(44,2)=54;
  nzij_pred(45,1)=55; nzij_pred(45,2)=56;
  nzij_pred(46,1)=56; nzij_pred(46,2)=57;
  nzij_pred(47,1)=57; nzij_pred(47,2)=58;
  nzij_pred(48,1)=58; nzij_pred(48,2)=59;
  nzij_pred(49,1)=59; nzij_pred(49,2)=60;
  nzij_pred(50,1)=60; nzij_pred(50,2)=61;
  nzij_pred(51,1)=61; nzij_pred(51,2)=62;
  nzij_pred(52,1)=62; nzij_pred(52,2)=63;
  nzij_pred(53,1)=63; nzij_pred(53,2)=64;
  nzij_pred(54,1)=64; nzij_pred(54,2)=65;
  nzij_pred(55,1)=65; nzij_pred(55,2)=66;
  nzij_pred(56,1)=66; nzij_pred(56,2)=67;
  nzij_pred(57,1)=67; nzij_pred(57,2)=68;
  nzij_pred(58,1)=68; nzij_pred(58,2)=69;
  nzij_pred(59,1)=69; nzij_pred(59,2)=70;
  nzij_pred(60,1)=70; nzij_pred(60,2)=71;
  nzij_pred(61,1)=71; nzij_pred(61,2)=72;
  nzij_pred(62,1)=72; nzij_pred(62,2)=73;
  nzij_pred(63,1)=73; nzij_pred(63,2)=74;
  nzij_pred(64,1)=74; nzij_pred(64,2)=75;
  nzij_pred(65,1)=75; nzij_pred(65,2)=76;
  nzij_pred(66,1)=76; nzij_pred(66,2)=77;
  nzij_pred(67,1)=77; nzij_pred(67,2)=78;
  nzij_pred(68,1)=78; nzij_pred(68,2)=79;
  nzij_pred(69,1)=79; nzij_pred(69,2)=80;
  nzij_pred(70,1)=80; nzij_pred(70,2)=81;
  nzij_pred(71,1)=81; nzij_pred(71,2)=82;
  nzij_pred(72,1)=82; nzij_pred(72,2)=83;
  nzij_pred(73,1)=83; nzij_pred(73,2)=84;
  nzij_pred(74,1)=84; nzij_pred(74,2)=85;
  nzij_pred(75,1)=85; nzij_pred(75,2)=86;
  nzij_pred(76,1)=86; nzij_pred(76,2)=87;
  nzij_pred(77,1)=87; nzij_pred(77,2)=88;
  nzij_pred(78,1)=88; nzij_pred(78,2)=89;
  nzij_pred(79,1)=89; nzij_pred(79,2)=90;
  nzij_pred(80,1)=90; nzij_pred(80,2)=91;
  nzij_pred(81,1)=92; nzij_pred(81,2)=93;
  nzij_pred(82,1)=93; nzij_pred(82,2)=94;
  nzij_pred(83,1)=94; nzij_pred(83,2)=95;
  nzij_pred(84,1)=95; nzij_pred(84,2)=96;
  nzij_pred(85,1)=96; nzij_pred(85,2)=97;
  nzij_pred(86,1)=97; nzij_pred(86,2)=98;
  nzij_pred(87,1)=98; nzij_pred(87,2)=99;
  nzij_pred(88,1)=99; nzij_pred(88,2)=100;
  nzij_pred(89,1)=100; nzij_pred(89,2)=101;
  nzij_pred(90,1)=101; nzij_pred(90,2)=102;
  nzij_pred(91,1)=102; nzij_pred(91,2)=103;
  nzij_pred(92,1)=103; nzij_pred(92,2)=104;
  nzij_pred(93,1)=104; nzij_pred(93,2)=105;
  nzij_pred(94,1)=105; nzij_pred(94,2)=106;
  nzij_pred(95,1)=106; nzij_pred(95,2)=107;
  nzij_pred(96,1)=107; nzij_pred(96,2)=108;
  nzij_pred(97,1)=108; nzij_pred(97,2)=109;
  nzij_pred(98,1)=109; nzij_pred(98,2)=110;
  nzij_pred(99,1)=110; nzij_pred(99,2)=111;
  nzij_pred(100,1)=111; nzij_pred(100,2)=112;
  nzij_pred(101,1)=112; nzij_pred(101,2)=113;
  nzij_pred(102,1)=113; nzij_pred(102,2)=114;
  nzij_pred(103,1)=114; nzij_pred(103,2)=115;
  nzij_pred(104,1)=115; nzij_pred(104,2)=116;
  nzij_pred(105,1)=116; nzij_pred(105,2)=117;
  nzij_pred(106,1)=117; nzij_pred(106,2)=118;
  nzij_pred(107,1)=118; nzij_pred(107,2)=119;
  nzij_pred(108,1)=119; nzij_pred(108,2)=120;
  nzij_pred(109,1)=120; nzij_pred(109,2)=121;
  nzij_pred(110,1)=121; nzij_pred(110,2)=122;
  nzij_pred(111,1)=122; nzij_pred(111,2)=123;
  nzij_pred(112,1)=123; nzij_pred(112,2)=124;
  nzij_pred(113,1)=124; nzij_pred(113,2)=125;
  nzij_pred(114,1)=125; nzij_pred(114,2)=126;
  nzij_pred(115,1)=126; nzij_pred(115,2)=127;
  nzij_pred(116,1)=127; nzij_pred(116,2)=128;
  nzij_current = zeros(150, 2, 'int32');
  nzij_current(1,1)=8; nzij_current(1,2)=1;
  nzij_current(2,1)=10; nzij_current(2,2)=1;
  nzij_current(3,1)=17; nzij_current(3,2)=1;
  nzij_current(4,1)=17; nzij_current(4,2)=2;
  nzij_current(5,1)=6; nzij_current(5,2)=3;
  nzij_current(6,1)=7; nzij_current(6,2)=3;
  nzij_current(7,1)=8; nzij_current(7,2)=3;
  nzij_current(8,1)=9; nzij_current(8,2)=3;
  nzij_current(9,1)=10; nzij_current(9,2)=4;
  nzij_current(10,1)=13; nzij_current(10,2)=4;
  nzij_current(11,1)=1; nzij_current(11,2)=5;
  nzij_current(12,1)=13; nzij_current(12,2)=5;
  nzij_current(13,1)=12; nzij_current(13,2)=6;
  nzij_current(14,1)=16; nzij_current(14,2)=7;
  nzij_current(15,1)=4; nzij_current(15,2)=8;
  nzij_current(16,1)=1; nzij_current(16,2)=9;
  nzij_current(17,1)=2; nzij_current(17,2)=9;
  nzij_current(18,1)=2; nzij_current(18,2)=10;
  nzij_current(19,1)=3; nzij_current(19,2)=10;
  nzij_current(20,1)=2; nzij_current(20,2)=11;
  nzij_current(21,1)=7; nzij_current(21,2)=12;
  nzij_current(22,1)=5; nzij_current(22,2)=13;
  nzij_current(23,1)=6; nzij_current(23,2)=13;
  nzij_current(24,1)=7; nzij_current(24,2)=14;
  nzij_current(25,1)=9; nzij_current(25,2)=14;
  nzij_current(26,1)=11; nzij_current(26,2)=14;
  nzij_current(27,1)=9; nzij_current(27,2)=15;
  nzij_current(28,1)=15; nzij_current(28,2)=15;
  nzij_current(29,1)=5; nzij_current(29,2)=16;
  nzij_current(30,1)=6; nzij_current(30,2)=16;
  nzij_current(31,1)=9; nzij_current(31,2)=16;
  nzij_current(32,1)=10; nzij_current(32,2)=16;
  nzij_current(33,1)=1; nzij_current(33,2)=17;
  nzij_current(34,1)=13; nzij_current(34,2)=17;
  nzij_current(35,1)=14; nzij_current(35,2)=17;
  nzij_current(36,1)=16; nzij_current(36,2)=17;
  nzij_current(37,1)=17; nzij_current(37,2)=17;
  nzij_current(38,1)=18; nzij_current(38,2)=18;
  nzij_current(39,1)=19; nzij_current(39,2)=19;
  nzij_current(40,1)=20; nzij_current(40,2)=20;
  nzij_current(41,1)=21; nzij_current(41,2)=21;
  nzij_current(42,1)=22; nzij_current(42,2)=22;
  nzij_current(43,1)=23; nzij_current(43,2)=23;
  nzij_current(44,1)=24; nzij_current(44,2)=24;
  nzij_current(45,1)=25; nzij_current(45,2)=25;
  nzij_current(46,1)=26; nzij_current(46,2)=26;
  nzij_current(47,1)=27; nzij_current(47,2)=27;
  nzij_current(48,1)=28; nzij_current(48,2)=28;
  nzij_current(49,1)=29; nzij_current(49,2)=29;
  nzij_current(50,1)=30; nzij_current(50,2)=30;
  nzij_current(51,1)=31; nzij_current(51,2)=31;
  nzij_current(52,1)=32; nzij_current(52,2)=32;
  nzij_current(53,1)=33; nzij_current(53,2)=33;
  nzij_current(54,1)=34; nzij_current(54,2)=34;
  nzij_current(55,1)=35; nzij_current(55,2)=35;
  nzij_current(56,1)=36; nzij_current(56,2)=36;
  nzij_current(57,1)=37; nzij_current(57,2)=37;
  nzij_current(58,1)=38; nzij_current(58,2)=38;
  nzij_current(59,1)=39; nzij_current(59,2)=39;
  nzij_current(60,1)=40; nzij_current(60,2)=40;
  nzij_current(61,1)=41; nzij_current(61,2)=41;
  nzij_current(62,1)=42; nzij_current(62,2)=42;
  nzij_current(63,1)=43; nzij_current(63,2)=43;
  nzij_current(64,1)=44; nzij_current(64,2)=44;
  nzij_current(65,1)=45; nzij_current(65,2)=45;
  nzij_current(66,1)=46; nzij_current(66,2)=46;
  nzij_current(67,1)=47; nzij_current(67,2)=47;
  nzij_current(68,1)=48; nzij_current(68,2)=48;
  nzij_current(69,1)=49; nzij_current(69,2)=49;
  nzij_current(70,1)=50; nzij_current(70,2)=50;
  nzij_current(71,1)=51; nzij_current(71,2)=51;
  nzij_current(72,1)=52; nzij_current(72,2)=52;
  nzij_current(73,1)=53; nzij_current(73,2)=53;
  nzij_current(74,1)=54; nzij_current(74,2)=54;
  nzij_current(75,1)=11; nzij_current(75,2)=55;
  nzij_current(76,1)=55; nzij_current(76,2)=55;
  nzij_current(77,1)=56; nzij_current(77,2)=56;
  nzij_current(78,1)=57; nzij_current(78,2)=57;
  nzij_current(79,1)=58; nzij_current(79,2)=58;
  nzij_current(80,1)=59; nzij_current(80,2)=59;
  nzij_current(81,1)=60; nzij_current(81,2)=60;
  nzij_current(82,1)=61; nzij_current(82,2)=61;
  nzij_current(83,1)=62; nzij_current(83,2)=62;
  nzij_current(84,1)=63; nzij_current(84,2)=63;
  nzij_current(85,1)=64; nzij_current(85,2)=64;
  nzij_current(86,1)=65; nzij_current(86,2)=65;
  nzij_current(87,1)=66; nzij_current(87,2)=66;
  nzij_current(88,1)=67; nzij_current(88,2)=67;
  nzij_current(89,1)=68; nzij_current(89,2)=68;
  nzij_current(90,1)=69; nzij_current(90,2)=69;
  nzij_current(91,1)=70; nzij_current(91,2)=70;
  nzij_current(92,1)=71; nzij_current(92,2)=71;
  nzij_current(93,1)=72; nzij_current(93,2)=72;
  nzij_current(94,1)=73; nzij_current(94,2)=73;
  nzij_current(95,1)=74; nzij_current(95,2)=74;
  nzij_current(96,1)=75; nzij_current(96,2)=75;
  nzij_current(97,1)=76; nzij_current(97,2)=76;
  nzij_current(98,1)=77; nzij_current(98,2)=77;
  nzij_current(99,1)=78; nzij_current(99,2)=78;
  nzij_current(100,1)=79; nzij_current(100,2)=79;
  nzij_current(101,1)=80; nzij_current(101,2)=80;
  nzij_current(102,1)=81; nzij_current(102,2)=81;
  nzij_current(103,1)=82; nzij_current(103,2)=82;
  nzij_current(104,1)=83; nzij_current(104,2)=83;
  nzij_current(105,1)=84; nzij_current(105,2)=84;
  nzij_current(106,1)=85; nzij_current(106,2)=85;
  nzij_current(107,1)=86; nzij_current(107,2)=86;
  nzij_current(108,1)=87; nzij_current(108,2)=87;
  nzij_current(109,1)=88; nzij_current(109,2)=88;
  nzij_current(110,1)=89; nzij_current(110,2)=89;
  nzij_current(111,1)=90; nzij_current(111,2)=90;
  nzij_current(112,1)=91; nzij_current(112,2)=91;
  nzij_current(113,1)=14; nzij_current(113,2)=92;
  nzij_current(114,1)=92; nzij_current(114,2)=92;
  nzij_current(115,1)=93; nzij_current(115,2)=93;
  nzij_current(116,1)=94; nzij_current(116,2)=94;
  nzij_current(117,1)=95; nzij_current(117,2)=95;
  nzij_current(118,1)=96; nzij_current(118,2)=96;
  nzij_current(119,1)=97; nzij_current(119,2)=97;
  nzij_current(120,1)=98; nzij_current(120,2)=98;
  nzij_current(121,1)=99; nzij_current(121,2)=99;
  nzij_current(122,1)=100; nzij_current(122,2)=100;
  nzij_current(123,1)=101; nzij_current(123,2)=101;
  nzij_current(124,1)=102; nzij_current(124,2)=102;
  nzij_current(125,1)=103; nzij_current(125,2)=103;
  nzij_current(126,1)=104; nzij_current(126,2)=104;
  nzij_current(127,1)=105; nzij_current(127,2)=105;
  nzij_current(128,1)=106; nzij_current(128,2)=106;
  nzij_current(129,1)=107; nzij_current(129,2)=107;
  nzij_current(130,1)=108; nzij_current(130,2)=108;
  nzij_current(131,1)=109; nzij_current(131,2)=109;
  nzij_current(132,1)=110; nzij_current(132,2)=110;
  nzij_current(133,1)=111; nzij_current(133,2)=111;
  nzij_current(134,1)=112; nzij_current(134,2)=112;
  nzij_current(135,1)=113; nzij_current(135,2)=113;
  nzij_current(136,1)=114; nzij_current(136,2)=114;
  nzij_current(137,1)=115; nzij_current(137,2)=115;
  nzij_current(138,1)=116; nzij_current(138,2)=116;
  nzij_current(139,1)=117; nzij_current(139,2)=117;
  nzij_current(140,1)=118; nzij_current(140,2)=118;
  nzij_current(141,1)=119; nzij_current(141,2)=119;
  nzij_current(142,1)=120; nzij_current(142,2)=120;
  nzij_current(143,1)=121; nzij_current(143,2)=121;
  nzij_current(144,1)=122; nzij_current(144,2)=122;
  nzij_current(145,1)=123; nzij_current(145,2)=123;
  nzij_current(146,1)=124; nzij_current(146,2)=124;
  nzij_current(147,1)=125; nzij_current(147,2)=125;
  nzij_current(148,1)=126; nzij_current(148,2)=126;
  nzij_current(149,1)=127; nzij_current(149,2)=127;
  nzij_current(150,1)=128; nzij_current(150,2)=128;
  nzij_fwrd = zeros(6, 2, 'int32');
  nzij_fwrd(1,1)=1; nzij_fwrd(1,2)=5;
  nzij_fwrd(2,1)=3; nzij_fwrd(2,2)=9;
  nzij_fwrd(3,1)=3; nzij_fwrd(3,2)=10;
  nzij_fwrd(4,1)=3; nzij_fwrd(4,2)=13;
  nzij_fwrd(5,1)=3; nzij_fwrd(5,2)=16;
  nzij_fwrd(6,1)=1; nzij_fwrd(6,2)=17;
end