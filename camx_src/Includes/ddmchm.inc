c----CAMx v7.32 250801
c
c     DDMCHM.INC sets species pointers for the DDM chemistry
c     These equivalences must be consistent with the internal
c     species lists defined in data statements in READCHM
c
c      Copyright 1996 - 2025
c     Ramboll
c
c     Modifications:
c      01/20/21 -gwilson-  Added "blank" pointers for traditional tracer species
c
      integer lmap(NSPNAM)
c
      common /lname/ lmap
c
      integer      lapo2,   lautx,   lbzo2
      integer      lc2o3,    lcro,   lcxo3
      integer      lddo2,   lepx2,    lho2
      integer         li,     lio,   liso2
      integer      lmeo2,      lo,    lo1d
      integer        loh,    loio,   lopo3
      integer       lro2,    lror,    lto2
      integer      ltpo2,   lxlo2,    lxo2
      integer      lxo2h,   lxo2n,   lxpar
      integer      lxprp,     lbr,    lbro
      integer        lcl,    lclo,   lhco3
      integer      larpx,   lddca,   ledoh
      integer      lethr,   lhact,   lhket
      integer      libta,   lipoh,   lnpoh
      integer      lpdoh,    lroh,   lsxd5
      integer      ltprd,   lcgly,   laacd
      integer      lacet,   lald2,   laldx
      integer      lbenz,   lbrcl,   lcat1
      integer      lclad,   lclao,     lco
      integer      lcres,   lcron,    ldms
      integer      lech4,   lepox,    leth
      integer      letha,   lethy,   letoh
      integer      lfacd,   lform,    lgly
      integer      lglyd,   lh2o2,    lhcl
      integer        lhi,   lhio3,   lhno3
      integer       lhoi,   lhono,   lhpld
      integer        li2,   li2o2,    libr
      integer       licl,   lino2,   lino3
      integer      lintr,   liole,   lisop
      integer      lispd,   lispx,   lixoy
      integer       lket,   lmeoh,   lmepx
      integer      lmgly,   ln2o5,    lnh3
      integer        lno,    lno2,    lno3
      integer      lntr1,   lntr2,     lo3
      integer       lole,   lopan,   lopen
      integer      lpacd,    lpan,   lpanx
      integer       lpar,    lpna,   lprpa
      integer      lrooh,    lso2,    lsqt
      integer      lsulf,   lterp,    ltol
      integer      lxopn,    lxyl,   lccrs
      integer      lacg1,   lacg2,   lbcg1
      integer      lbcg2,   lcprm,   lfcrs
      integer      lfprm,   livoa,   lsvoa
      integer        lna,    lpal,    lpca
      integer       lpcl,    lpec,    lpfe
      integer      lph2o,     lpk,    lpmg
      integer       lpmn,   lpnh4,   lpno3
      integer       lhoa,    lpsi,   lpso4
      integer       lpti,   laoa1,   laoa2
      integer      lboa1,   lboa2,   laoa0
      integer      lboa0,    lbr2,   lbrn2
      integer      lbrn3,   lch3i,    lcl2
      integer      lcln2,   lcln3,   lfmbr
      integer      lfmcl,    lhbr,   lhobr
      integer      lhocl,   lialk,    lmb2
      integer      lmb2c,    lmb3,    lmbc
      integer      lmbc2,    lmi2,    lmib
      integer       lmic,   lssbr,   lsscl
      integer      lssn3,    lhg0,    lhg2
      integer     lhgiip, lhgiipc,    lhgp
      integer      livob,   livod,   livog
      integer      lpap0,   lpap1,   lpap2
      integer      lpap3,   lpap4,   lpas0
      integer      lpas1,   lpas2,   lpas3
      integer      lpas4,   lpbs0,   lpbs1
      integer      lpbs2,   lpbs3,   lpbs4
      integer      lpcp0,   lpcp1,   lpcp2
      integer      lpcp3,   lpcp4,   lpfp0
      integer      lpfp1,   lpfp2,   lpfp3
      integer      lpfp4,   lvap1,   lvap2
      integer      lvap3,   lvap4,   lvas1
      integer      lvas2,   lvas3,   lvas4
      integer      lvbs1,   lvbs2,   lvbs3
      integer      lvbs4,   lvcp1,   lvcp2
      integer      lvcp3,   lvcp4,   lvfp1
      integer      lvfp2,   lvfp3,   lvfp4
      integer      laco3,   lactp,   ladcn
      integer      laddc,   lapip,   lbal1
      integer      lbal2,   lbalp,   lbenp
      integer       lcho,   letep,   lethp
      integer      lhc3p,   lhc5p,   lhc8p
      integer      lketp,   llimp,   lmacp
      integer       lmcp,   lmcto,   lmctp
      integer      lmekp,    lmo2,   lmvkp
      integer      lolip,   lolnd,   lolnn
      integer      loltp,   lorap,   lper1
      integer      lper2,   ltlp1,   ltolp
      integer       ltr2,   lualp,    lxy2
      integer      lxyl1,   lxylp,   lxyo2
      integer      lxyop,    lacd,    lace
      integer       lact,    lald,    lapi
      integer       lben,    lcsl,   ldcb1
      integer      ldcb2,   ldcb3,   ldien
      integer       leoh,    lepx,    lete
      integer      leteg,    lhc3,    lhc5
      integer       lhc8,   lhno4,   lishp
      integer       liso,   lison,    llim
      integer      lmahp,    lmct,    lmoh
      integer      lnald,    loli,    lolt
      integer      lonit,    lop1,    lop2
      integer      lora1,   lora2,    lpaa
      integer      lphen,    lppn,   luald
      integer       lxym,    lxyo,    lxyp
      integer      lbzc3,    lbzo,   lmac3
      integer      lmco3,    lo3p,   lrco3
      integer      lro2c,   lro2x,   ltbuo
      integer      lxace,   lxacr,   lxaf1
      integer      lxaf2,   lxaf3,   lxbac
      integer      lxbal,   lxcch,    lxco
      integer      lxgld,   lxgly,   lxhch
      integer      lxho2,   lxipr,   lxma3
      integer      lxmac,   lxmc3,   lxmek
      integer      lxmeo,   lxmgl,   lxmvk
      integer      lxno2,    lxoh,   lxpd2
      integer      lxrc3,   lxrch,   lxrn3
      integer      lxtbu,   ly6px,   lyapx
      integer      lyrpx,   lzrn3,   lacro
      integer      lacye,   lafg1,   lafg2
      integer      lafg3,   lalk1,   lalk2
      integer      lalk3,   lalk4,   lalk5
      integer      lapin,   laro1,   laro2
      integer      lb124,   lbacl,   lbald
      integer      lbd13,   lccho,   lco3h
      integer      lcooh,   lethe,   lhcho
      integer      liprd,   lmacr,    lmek
      integer      lmpan,    lmvk,   lmxyl
      integer      lnphe,   lole1,   lole2
      integer      loxyl,   lpan2,   lpbzn
      integer      lprd2,   lprpe,   lpxyl
      integer      lr6px,   lrapx,   lrcho
      integer      lrno3,   lro3h,   lsesq
      integer      ltolu,     lxn,   ltr01
      integer      ltr02,   ltr03,   ltr04
      integer      ltr05,   ltr06,   ltr07
      integer      ltr08,   ltr09
c
      equivalence
     &   (lmap(  1),   lapo2),(lmap(  2),   lautx),(lmap(  3),   lbzo2),
     &   (lmap(  4),   lc2o3),(lmap(  5),    lcro),(lmap(  6),   lcxo3),
     &   (lmap(  7),   lddo2),(lmap(  8),   lepx2),(lmap(  9),    lho2),
     &   (lmap( 10),      li),(lmap( 11),     lio),(lmap( 12),   liso2),
     &   (lmap( 13),   lmeo2),(lmap( 14),      lo),(lmap( 15),    lo1d),
     &   (lmap( 16),     loh),(lmap( 17),    loio),(lmap( 18),   lopo3),
     &   (lmap( 19),    lro2),(lmap( 20),    lror),(lmap( 21),    lto2),
     &   (lmap( 22),   ltpo2),(lmap( 23),   lxlo2),(lmap( 24),    lxo2),
     &   (lmap( 25),   lxo2h),(lmap( 26),   lxo2n),(lmap( 27),   lxpar),
     &   (lmap( 28),   lxprp),(lmap( 29),     lbr),(lmap( 30),    lbro),
     &   (lmap( 31),     lcl),(lmap( 32),    lclo),(lmap( 33),   lhco3),
     &   (lmap( 34),   larpx),(lmap( 35),   lddca),(lmap( 36),   ledoh),
     &   (lmap( 37),   lethr),(lmap( 38),   lhact),(lmap( 39),   lhket),
     &   (lmap( 40),   libta),(lmap( 41),   lipoh),(lmap( 42),   lnpoh),
     &   (lmap( 43),   lpdoh),(lmap( 44),    lroh),(lmap( 45),   lsxd5),
     &   (lmap( 46),   ltprd),(lmap( 47),   lcgly),(lmap( 48),   laacd),
     &   (lmap( 49),   lacet),(lmap( 50),   lald2),(lmap( 51),   laldx),
     &   (lmap( 52),   lbenz),(lmap( 53),   lbrcl),(lmap( 54),   lcat1),
     &   (lmap( 55),   lclad),(lmap( 56),   lclao),(lmap( 57),     lco),
     &   (lmap( 58),   lcres),(lmap( 59),   lcron),(lmap( 60),    ldms),
     &   (lmap( 61),   lech4),(lmap( 62),   lepox),(lmap( 63),    leth),
     &   (lmap( 64),   letha),(lmap( 65),   lethy),(lmap( 66),   letoh),
     &   (lmap( 67),   lfacd),(lmap( 68),   lform),(lmap( 69),    lgly),
     &   (lmap( 70),   lglyd),(lmap( 71),   lh2o2),(lmap( 72),    lhcl),
     &   (lmap( 73),     lhi),(lmap( 74),   lhio3),(lmap( 75),   lhno3),
     &   (lmap( 76),    lhoi),(lmap( 77),   lhono),(lmap( 78),   lhpld),
     &   (lmap( 79),     li2),(lmap( 80),   li2o2),(lmap( 81),    libr),
     &   (lmap( 82),    licl),(lmap( 83),   lino2),(lmap( 84),   lino3),
     &   (lmap( 85),   lintr),(lmap( 86),   liole),(lmap( 87),   lisop),
     &   (lmap( 88),   lispd),(lmap( 89),   lispx),(lmap( 90),   lixoy),
     &   (lmap( 91),    lket),(lmap( 92),   lmeoh),(lmap( 93),   lmepx),
     &   (lmap( 94),   lmgly),(lmap( 95),   ln2o5),(lmap( 96),    lnh3),
     &   (lmap( 97),     lno),(lmap( 98),    lno2),(lmap( 99),    lno3),
     &   (lmap(100),   lntr1),(lmap(101),   lntr2),(lmap(102),     lo3),
     &   (lmap(103),    lole),(lmap(104),   lopan),(lmap(105),   lopen),
     &   (lmap(106),   lpacd),(lmap(107),    lpan),(lmap(108),   lpanx),
     &   (lmap(109),    lpar),(lmap(110),    lpna),(lmap(111),   lprpa),
     &   (lmap(112),   lrooh),(lmap(113),    lso2),(lmap(114),    lsqt),
     &   (lmap(115),   lsulf),(lmap(116),   lterp),(lmap(117),    ltol),
     &   (lmap(118),   lxopn),(lmap(119),    lxyl),(lmap(120),   lccrs),
     &   (lmap(121),   lacg1),(lmap(122),   lacg2),(lmap(123),   lbcg1),
     &   (lmap(124),   lbcg2),(lmap(125),   lcprm),(lmap(126),   lfcrs),
     &   (lmap(127),   lfprm),(lmap(128),   livoa),(lmap(129),   lsvoa),
     &   (lmap(130),     lna),(lmap(131),    lpal),(lmap(132),    lpca),
     &   (lmap(133),    lpcl),(lmap(134),    lpec),(lmap(135),    lpfe),
     &   (lmap(136),   lph2o),(lmap(137),     lpk),(lmap(138),    lpmg),
     &   (lmap(139),    lpmn),(lmap(140),   lpnh4),(lmap(141),   lpno3),
     &   (lmap(142),    lhoa),(lmap(143),    lpsi),(lmap(144),   lpso4),
     &   (lmap(145),    lpti),(lmap(146),   laoa1),(lmap(147),   laoa2),
     &   (lmap(148),   lboa1),(lmap(149),   lboa2),(lmap(150),   laoa0),
     &   (lmap(151),   lboa0),(lmap(152),    lbr2),(lmap(153),   lbrn2),
     &   (lmap(154),   lbrn3),(lmap(155),   lch3i),(lmap(156),    lcl2),
     &   (lmap(157),   lcln2),(lmap(158),   lcln3),(lmap(159),   lfmbr),
     &   (lmap(160),   lfmcl),(lmap(161),    lhbr),(lmap(162),   lhobr),
     &   (lmap(163),   lhocl),(lmap(164),   lialk),(lmap(165),    lmb2),
     &   (lmap(166),   lmb2c),(lmap(167),    lmb3),(lmap(168),    lmbc),
     &   (lmap(169),   lmbc2),(lmap(170),    lmi2),(lmap(171),    lmib),
     &   (lmap(172),    lmic),(lmap(173),   lssbr),(lmap(174),   lsscl),
     &   (lmap(175),   lssn3),(lmap(176),    lhg0),(lmap(177),    lhg2),
     &   (lmap(178),  lhgiip),(lmap(179), lhgiipc),(lmap(180),    lhgp),
     &   (lmap(181),   livob),(lmap(182),   livod),(lmap(183),   livog),
     &   (lmap(184),   lpap0),(lmap(185),   lpap1),(lmap(186),   lpap2),
     &   (lmap(187),   lpap3),(lmap(188),   lpap4),(lmap(189),   lpas0),
     &   (lmap(190),   lpas1),(lmap(191),   lpas2),(lmap(192),   lpas3),
     &   (lmap(193),   lpas4),(lmap(194),   lpbs0),(lmap(195),   lpbs1),
     &   (lmap(196),   lpbs2),(lmap(197),   lpbs3),(lmap(198),   lpbs4),
     &   (lmap(199),   lpcp0),(lmap(200),   lpcp1),(lmap(201),   lpcp2),
     &   (lmap(202),   lpcp3),(lmap(203),   lpcp4),(lmap(204),   lpfp0),
     &   (lmap(205),   lpfp1),(lmap(206),   lpfp2),(lmap(207),   lpfp3),
     &   (lmap(208),   lpfp4),(lmap(209),   lvap1),(lmap(210),   lvap2),
     &   (lmap(211),   lvap3),(lmap(212),   lvap4),(lmap(213),   lvas1),
     &   (lmap(214),   lvas2),(lmap(215),   lvas3),(lmap(216),   lvas4),
     &   (lmap(217),   lvbs1),(lmap(218),   lvbs2),(lmap(219),   lvbs3),
     &   (lmap(220),   lvbs4),(lmap(221),   lvcp1),(lmap(222),   lvcp2),
     &   (lmap(223),   lvcp3),(lmap(224),   lvcp4),(lmap(225),   lvfp1),
     &   (lmap(226),   lvfp2),(lmap(227),   lvfp3),(lmap(228),   lvfp4),
     &   (lmap(229),   laco3),(lmap(230),   lactp),(lmap(231),   ladcn),
     &   (lmap(232),   laddc),(lmap(233),   lapip),(lmap(234),   lbal1),
     &   (lmap(235),   lbal2),(lmap(236),   lbalp),(lmap(237),   lbenp),
     &   (lmap(238),    lcho),(lmap(239),   letep),(lmap(240),   lethp),
     &   (lmap(241),   lhc3p),(lmap(242),   lhc5p),(lmap(243),   lhc8p),
     &   (lmap(244),   lketp),(lmap(245),   llimp),(lmap(246),   lmacp),
     &   (lmap(247),    lmcp),(lmap(248),   lmcto),(lmap(249),   lmctp),
     &   (lmap(250),   lmekp),(lmap(251),    lmo2),(lmap(252),   lmvkp),
     &   (lmap(253),   lolip),(lmap(254),   lolnd),(lmap(255),   lolnn),
     &   (lmap(256),   loltp),(lmap(257),   lorap),(lmap(258),   lper1),
     &   (lmap(259),   lper2),(lmap(260),   ltlp1),(lmap(261),   ltolp),
     &   (lmap(262),    ltr2),(lmap(263),   lualp),(lmap(264),    lxy2),
     &   (lmap(265),   lxyl1),(lmap(266),   lxylp),(lmap(267),   lxyo2),
     &   (lmap(268),   lxyop),(lmap(269),    lacd),(lmap(270),    lace)
      equivalence
     &   (lmap(271),    lact),(lmap(272),    lald),(lmap(273),    lapi),
     &   (lmap(274),    lben),(lmap(275),    lcsl),(lmap(276),   ldcb1),
     &   (lmap(277),   ldcb2),(lmap(278),   ldcb3),(lmap(279),   ldien),
     &   (lmap(280),    leoh),(lmap(281),    lepx),(lmap(282),    lete),
     &   (lmap(283),   leteg),(lmap(284),    lhc3),(lmap(285),    lhc5),
     &   (lmap(286),    lhc8),(lmap(287),   lhno4),(lmap(288),   lishp),
     &   (lmap(289),    liso),(lmap(290),   lison),(lmap(291),    llim),
     &   (lmap(292),   lmahp),(lmap(293),    lmct),(lmap(294),    lmoh),
     &   (lmap(295),   lnald),(lmap(296),    loli),(lmap(297),    lolt),
     &   (lmap(298),   lonit),(lmap(299),    lop1),(lmap(300),    lop2),
     &   (lmap(301),   lora1),(lmap(302),   lora2),(lmap(303),    lpaa),
     &   (lmap(304),   lphen),(lmap(305),    lppn),(lmap(306),   luald),
     &   (lmap(307),    lxym),(lmap(308),    lxyo),(lmap(309),    lxyp),
     &   (lmap(310),   lbzc3),(lmap(311),    lbzo),(lmap(312),   lmac3),
     &   (lmap(313),   lmco3),(lmap(314),    lo3p),(lmap(315),   lrco3),
     &   (lmap(316),   lro2c),(lmap(317),   lro2x),(lmap(318),   ltbuo),
     &   (lmap(319),   lxace),(lmap(320),   lxacr),(lmap(321),   lxaf1),
     &   (lmap(322),   lxaf2),(lmap(323),   lxaf3),(lmap(324),   lxbac),
     &   (lmap(325),   lxbal),(lmap(326),   lxcch),(lmap(327),    lxco),
     &   (lmap(328),   lxgld),(lmap(329),   lxgly),(lmap(330),   lxhch),
     &   (lmap(331),   lxho2),(lmap(332),   lxipr),(lmap(333),   lxma3),
     &   (lmap(334),   lxmac),(lmap(335),   lxmc3),(lmap(336),   lxmek),
     &   (lmap(337),   lxmeo),(lmap(338),   lxmgl),(lmap(339),   lxmvk),
     &   (lmap(340),   lxno2),(lmap(341),    lxoh),(lmap(342),   lxpd2),
     &   (lmap(343),   lxrc3),(lmap(344),   lxrch),(lmap(345),   lxrn3),
     &   (lmap(346),   lxtbu),(lmap(347),   ly6px),(lmap(348),   lyapx),
     &   (lmap(349),   lyrpx),(lmap(350),   lzrn3),(lmap(351),   lacro),
     &   (lmap(352),   lacye),(lmap(353),   lafg1),(lmap(354),   lafg2),
     &   (lmap(355),   lafg3),(lmap(356),   lalk1),(lmap(357),   lalk2),
     &   (lmap(358),   lalk3),(lmap(359),   lalk4),(lmap(360),   lalk5),
     &   (lmap(361),   lapin),(lmap(362),   laro1),(lmap(363),   laro2),
     &   (lmap(364),   lb124),(lmap(365),   lbacl),(lmap(366),   lbald),
     &   (lmap(367),   lbd13),(lmap(368),   lccho),(lmap(369),   lco3h),
     &   (lmap(370),   lcooh),(lmap(371),   lethe),(lmap(372),   lhcho),
     &   (lmap(373),   liprd),(lmap(374),   lmacr),(lmap(375),    lmek),
     &   (lmap(376),   lmpan),(lmap(377),    lmvk),(lmap(378),   lmxyl),
     &   (lmap(379),   lnphe),(lmap(380),   lole1),(lmap(381),   lole2),
     &   (lmap(382),   loxyl),(lmap(383),   lpan2),(lmap(384),   lpbzn),
     &   (lmap(385),   lprd2),(lmap(386),   lprpe),(lmap(387),   lpxyl),
     &   (lmap(388),   lr6px),(lmap(389),   lrapx),(lmap(390),   lrcho),
     &   (lmap(391),   lrno3),(lmap(392),   lro3h),(lmap(393),   lsesq),
     &   (lmap(394),   ltolu),(lmap(395),     lxn),(lmap(396),   ltr01),
     &   (lmap(397),   ltr02),(lmap(398),   ltr03),(lmap(399),   ltr04),
     &   (lmap(400),   ltr05),(lmap(401),   ltr06),(lmap(402),   ltr07),
     &   (lmap(403),   ltr08),(lmap(404),   ltr09)
