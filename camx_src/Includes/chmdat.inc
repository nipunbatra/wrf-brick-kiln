c----CAMx v7.32 250801
c  
c     CHMDAT.INC contains all chemistry variables 
c                            
c      Copyright 1996 - 2025
c     Ramboll 
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c      10/14/04    Modifications for mechanism 10
c      11/30/06    Created by taking Static variables out of chmstry.com
c      12/29/06 -bkoo-     Expanded species list for the updated SOA scheme
c      01/08/07 -bkoo-     Expanded species list for Mechanism 6 (CB05)
c                          Now ipigrxn is a scalar pointer to NO self-reaction
c      07/04/07 -bkoo-     Added pointer to hydrolysis of N2O5
c      04/28/08 -gyarwood- Added EBI solver
c      04/12/12 -cemery-   Added T and P adjustments to photolysis rates
c      07/20/12 -ou-       Added Mechanism 1 (CB06 with Iodine chemistry)
c      10/08/12 -jjung-    Added Mechanism 2 (CB6r1)
c      04/30/13 -cemery-   Added surface model
c      06/28/13 -bkoo-     Added indexes for new GREASD PiG chemistry
c      10/14/13 -bkoo-     Added pointer to hydrolysis of organic nitrate
c      10/24/13 -bkoo-     Revised Mechanism 2 (CB6r1 -> CB6r2)
c      03/18/14 -bkoo-     Added species for benzene SOA
c      04/07/14 -cemery-   Added top con file
c      06/27/14 -bkoo-     Added Mechanism 3 (CB6r2h)
c      11/11/14 -bkoo-     Added Mechanism 4 (CB6r3)
c      12/07/14 -bkoo-     Added species for VBS
c      01/08/16 -bkoo-     Updated for SAPRC99->SAPRC07TC
c                                      Revised VBS
c                                      Removal of IEH solver
c      02/29/16 -bkoo-     Added pointer to heterogeneous rxn of N2O5 + HCL
c      07/20/16 -bkoo-     Added pointer to heterogeneous hydrolysis of INTR
c      09/02/16 -bkoo-     Revised heterogeneous rxn pointers
c      11/28/16 -bkoo-     Revised VBS precursors
c      11/22/17 -bkoo-     Added PFE/PMN/PK/PCA/PMG
c      01/12/18 -bkoo-     Removed BNZA/TOLA/XYLA/ISP/TRP
c      01/06/19 -cemery-   Added PAL/PSI/PTI
c      01/09/19 -cemery-   Added DMS to CB6r4 (remains Mech 4)
c      11/17/20 -cemery-   Added definition of TUV photolysis file labels
c      01/20/21 -gwilson-  Added "blank" pointers for traditional tracer species
c      03/24/21 -gy-       Added CB7 (Mech 7)
c      08/10/21 -gy-       Added RACM2s21 (Mech 9)
c 
c-----------------------------------------------------------------------
c     Parameters for some of the switches:
c
c     CDCMC  -- string for requesting the CMC chemistry solver
c     CDEBI  -- string for requesting the EBI chemistry solver
c     CDIEH  -- string for requesting the IEH chemistry solver
c     CDLSOD -- string for requesting the LSODE chemistry solver
c     IDCMC  -- code for using the CMC chemistry solver
c     IDEBI  -- code for using the EBI chemistry solver
c     IDIEH  -- code for using the IEH chemistry solver
c     IDLSOD -- code for using the LSODE chemistry solver
c
c-----------------------------------------------------------------------
c
      character*10 CDCMC
      character*10 CDEBI
      character*10 CDIEH
      character*10 CDLSOD
      integer      IDCMC
      integer      IDEBI
      integer      IDIEH
      integer      IDLSOD
c
      parameter( CDCMC  = "CMC       " )
      parameter( CDEBI  = "EBI       " )
      parameter( CDIEH  = "IEH       " )
      parameter( CDLSOD = "LSODE     " )
      parameter( IDCMC  = 1 )
      parameter( IDEBI  = 2 )
      parameter( IDIEH  = 3 )
      parameter( IDLSOD = 4 )
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    ngas    --  number of gas species being modeled
c    naero   --  number of aersol species being modeled
c    nspec   --  total number of modeled species
c    nrad    --  number of radical species being modeled
c    nreact  --  number of chemical reactions
c    idmech  --  the code which determines which chemical mechanism is used
c    idsolv  --  the code which determines which chemstry solver to use
c    navspc  --  number of species to write to output average file
c    ndepspc --  number of species to write to deposition file
c    nicspc  --  number of species in the initial conditions file
c    nbcspc  --  number of species in the boundary conditions file
c    ntcspc  --  number of species in the top concentrations file
c    nemspc  --  number of species in the ALL emissions files
c    nsmspc  --  number of species in surface model
c    nsmrxn  --  number of reactions in surface model
c
c-----------------------------------------------------------------------
c
       integer      ngas
       integer      naero
       integer      nspec
       integer      nrad
       integer      nreact
       integer      idmech
       integer      idsolv
       integer      navspc
       integer      ndepspc
       integer      nicspc
       integer      nbcspc
       integer      ntcspc
       integer      nemspc
       integer      nsmspc
       integer      nsmrxn
c
      common /chm1/ ngas, naero, nspec, nrad, nreact,
     &                         idmech, idsolv, navspc, ndepspc, nicspc,
     &                           nbcspc, ntcspc, nsmspc, nsmrxn, nemspc
c
c-----------------------------------------------------------------------
c
c    aeropt --  keyword for aerosol scheme ('NONE', 'INERT', 'CF')
c    leqsam --  logical flag to use EQSAM chemistry
c
c-----------------------------------------------------------------------
c
      character*10 aeropt
      logical leqsam
c
      common /chm2/ aeropt
      common /chm2dat/ leqsam
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     nphot1   -- number of primary photolysis reactions
c     nphot2   -- number of secondary (scaled) photolysis reactions
c
c-----------------------------------------------------------------------
c
      integer nphot1
      integer nphot2
      character*20 tuvlabl(10)
c
      common /photmap/ nphot1, nphot2, tuvlabl
c 
c-----------------------------------------------------------------------
c     Mercury chemistry variables 
c     Variables for vertical profiles of Cl2, Br and BrO for Hg chemistry
c
c     NHTHAL   --  number of heights for halogen concentrations 
c     hthal    --  heights (km) for the top of each concentration range
c     cl2day  --  profile of Cl2 (ppm) with height during day over ocean
c     cl2nite --  profile of Cl2 (ppm) with height at night over ocean
c     brlprof  --  profile of Br (ppm) with height over land
c     brwprof  --  profile of Br (ppm) with height over ocean
c     brolprof --  profile of BrO (ppm) with height over land
c     browprof --  profile of BrO (ppm) with height over ocean
c
c     MWHG     --  Hg molecular weight 
c-----------------------------------------------------------------------
c
      integer, parameter :: NHTHAL = 6
      real, parameter :: MWHG = 200.59
c
      real, dimension(NHTHAL) :: hthal,brlprof,brwprof,brolprof,
     &                           browprof,cl2day,cl2nite
c
      common /halprofs/ hthal, brlprof, brwprof, brolprof, browprof,
     &                  cl2day, cl2nite
c 
c-----------------------------------------------------------------------
c     Parameters for mapping input species to internal model order:
c
c      NSPNAM  -- number of recognizable state species
c
      integer   NSPNAM
c
      parameter ( NSPNAM  =   404 )
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     kmap     -- mapping of species on chemistry parameters file to
c                 internal order
c
c-----------------------------------------------------------------------
c
      integer kmap(NSPNAM)
c
      common /kname/ kmap
c
      integer      kapo2,   kautx,   kbzo2
      integer      kc2o3,    kcro,   kcxo3
      integer      kddo2,   kepx2,    kho2
      integer         ki,     kio,   kiso2
      integer      kmeo2,      ko,    ko1d
      integer        koh,    koio,   kopo3
      integer       kro2,    kror,    kto2
      integer      ktpo2,   kxlo2,    kxo2
      integer      kxo2h,   kxo2n,   kxpar
      integer      kxprp,     kbr,    kbro
      integer        kcl,    kclo,   khco3
      integer      karpx,   kddca,   kedoh
      integer      kethr,   khact,   khket
      integer      kibta,   kipoh,   knpoh
      integer      kpdoh,    kroh,   ksxd5
      integer      ktprd,   kcgly,   kaacd
      integer      kacet,   kald2,   kaldx
      integer      kbenz,   kbrcl,   kcat1
      integer      kclad,   kclao,     kco
      integer      kcres,   kcron,    kdms
      integer      kech4,   kepox,    keth
      integer      ketha,   kethy,   ketoh
      integer      kfacd,   kform,    kgly
      integer      kglyd,   kh2o2,    khcl
      integer        khi,   khio3,   khno3
      integer       khoi,   khono,   khpld
      integer        ki2,   ki2o2,    kibr
      integer       kicl,   kino2,   kino3
      integer      kintr,   kiole,   kisop
      integer      kispd,   kispx,   kixoy
      integer       kket,   kmeoh,   kmepx
      integer      kmgly,   kn2o5,    knh3
      integer        kno,    kno2,    kno3
      integer      kntr1,   kntr2,     ko3
      integer       kole,   kopan,   kopen
      integer      kpacd,    kpan,   kpanx
      integer       kpar,    kpna,   kprpa
      integer      krooh,    kso2,    ksqt
      integer      ksulf,   kterp,    ktol
      integer      kxopn,    kxyl,   kccrs
      integer      kacg1,   kacg2,   kbcg1
      integer      kbcg2,   kcprm,   kfcrs
      integer      kfprm,   kivoa,   ksvoa
      integer        kna,    kpal,    kpca
      integer       kpcl,    kpec,    kpfe
      integer      kph2o,     kpk,    kpmg
      integer       kpmn,   kpnh4,   kpno3
      integer       khoa,    kpsi,   kpso4
      integer       kpti,   kaoa1,   kaoa2
      integer      kboa1,   kboa2,   kaoa0
      integer      kboa0,    kbr2,   kbrn2
      integer      kbrn3,   kch3i,    kcl2
      integer      kcln2,   kcln3,   kfmbr
      integer      kfmcl,    khbr,   khobr
      integer      khocl,   kialk,    kmb2
      integer      kmb2c,    kmb3,    kmbc
      integer      kmbc2,    kmi2,    kmib
      integer       kmic,   kssbr,   ksscl
      integer      kssn3,    khg0,    khg2
      integer     khgiip, khgiipc,    khgp
      integer      kivob,   kivod,   kivog
      integer      kpap0,   kpap1,   kpap2
      integer      kpap3,   kpap4,   kpas0
      integer      kpas1,   kpas2,   kpas3
      integer      kpas4,   kpbs0,   kpbs1
      integer      kpbs2,   kpbs3,   kpbs4
      integer      kpcp0,   kpcp1,   kpcp2
      integer      kpcp3,   kpcp4,   kpfp0
      integer      kpfp1,   kpfp2,   kpfp3
      integer      kpfp4,   kvap1,   kvap2
      integer      kvap3,   kvap4,   kvas1
      integer      kvas2,   kvas3,   kvas4
      integer      kvbs1,   kvbs2,   kvbs3
      integer      kvbs4,   kvcp1,   kvcp2
      integer      kvcp3,   kvcp4,   kvfp1
      integer      kvfp2,   kvfp3,   kvfp4
      integer      kaco3,   kactp,   kadcn
      integer      kaddc,   kapip,   kbal1
      integer      kbal2,   kbalp,   kbenp
      integer       kcho,   ketep,   kethp
      integer      khc3p,   khc5p,   khc8p
      integer      kketp,   klimp,   kmacp
      integer       kmcp,   kmcto,   kmctp
      integer      kmekp,    kmo2,   kmvkp
      integer      kolip,   kolnd,   kolnn
      integer      koltp,   korap,   kper1
      integer      kper2,   ktlp1,   ktolp
      integer       ktr2,   kualp,    kxy2
      integer      kxyl1,   kxylp,   kxyo2
      integer      kxyop,    kacd,    kace
      integer       kact,    kald,    kapi
      integer       kben,    kcsl,   kdcb1
      integer      kdcb2,   kdcb3,   kdien
      integer       keoh,    kepx,    kete
      integer      keteg,    khc3,    khc5
      integer       khc8,   khno4,   kishp
      integer       kiso,   kison,    klim
      integer      kmahp,    kmct,    kmoh
      integer      knald,    koli,    kolt
      integer      konit,    kop1,    kop2
      integer      kora1,   kora2,    kpaa
      integer      kphen,    kppn,   kuald
      integer       kxym,    kxyo,    kxyp
      integer      kbzc3,    kbzo,   kmac3
      integer      kmco3,    ko3p,   krco3
      integer      kro2c,   kro2x,   ktbuo
      integer      kxace,   kxacr,   kxaf1
      integer      kxaf2,   kxaf3,   kxbac
      integer      kxbal,   kxcch,    kxco
      integer      kxgld,   kxgly,   kxhch
      integer      kxho2,   kxipr,   kxma3
      integer      kxmac,   kxmc3,   kxmek
      integer      kxmeo,   kxmgl,   kxmvk
      integer      kxno2,    kxoh,   kxpd2
      integer      kxrc3,   kxrch,   kxrn3
      integer      kxtbu,   ky6px,   kyapx
      integer      kyrpx,   kzrn3,   kacro
      integer      kacye,   kafg1,   kafg2
      integer      kafg3,   kalk1,   kalk2
      integer      kalk3,   kalk4,   kalk5
      integer      kapin,   karo1,   karo2
      integer      kb124,   kbacl,   kbald
      integer      kbd13,   kccho,   kco3h
      integer      kcooh,   kethe,   khcho
      integer      kiprd,   kmacr,    kmek
      integer      kmpan,    kmvk,   kmxyl
      integer      knphe,   kole1,   kole2
      integer      koxyl,   kpan2,   kpbzn
      integer      kprd2,   kprpe,   kpxyl
      integer      kr6px,   krapx,   krcho
      integer      krno3,   kro3h,   ksesq
      integer      ktolu,     kxn,   ktr01
      integer      ktr02,   ktr03,   ktr04
      integer      ktr05,   ktr06,   ktr07
      integer      ktr08,   ktr09
c
      equivalence
     &   (kmap(  1),   kapo2),(kmap(  2),   kautx),(kmap(  3),   kbzo2),
     &   (kmap(  4),   kc2o3),(kmap(  5),    kcro),(kmap(  6),   kcxo3),
     &   (kmap(  7),   kddo2),(kmap(  8),   kepx2),(kmap(  9),    kho2),
     &   (kmap( 10),      ki),(kmap( 11),     kio),(kmap( 12),   kiso2),
     &   (kmap( 13),   kmeo2),(kmap( 14),      ko),(kmap( 15),    ko1d),
     &   (kmap( 16),     koh),(kmap( 17),    koio),(kmap( 18),   kopo3),
     &   (kmap( 19),    kro2),(kmap( 20),    kror),(kmap( 21),    kto2),
     &   (kmap( 22),   ktpo2),(kmap( 23),   kxlo2),(kmap( 24),    kxo2),
     &   (kmap( 25),   kxo2h),(kmap( 26),   kxo2n),(kmap( 27),   kxpar),
     &   (kmap( 28),   kxprp),(kmap( 29),     kbr),(kmap( 30),    kbro),
     &   (kmap( 31),     kcl),(kmap( 32),    kclo),(kmap( 33),   khco3),
     &   (kmap( 34),   karpx),(kmap( 35),   kddca),(kmap( 36),   kedoh),
     &   (kmap( 37),   kethr),(kmap( 38),   khact),(kmap( 39),   khket),
     &   (kmap( 40),   kibta),(kmap( 41),   kipoh),(kmap( 42),   knpoh),
     &   (kmap( 43),   kpdoh),(kmap( 44),    kroh),(kmap( 45),   ksxd5),
     &   (kmap( 46),   ktprd),(kmap( 47),   kcgly),(kmap( 48),   kaacd),
     &   (kmap( 49),   kacet),(kmap( 50),   kald2),(kmap( 51),   kaldx),
     &   (kmap( 52),   kbenz),(kmap( 53),   kbrcl),(kmap( 54),   kcat1),
     &   (kmap( 55),   kclad),(kmap( 56),   kclao),(kmap( 57),     kco),
     &   (kmap( 58),   kcres),(kmap( 59),   kcron),(kmap( 60),    kdms),
     &   (kmap( 61),   kech4),(kmap( 62),   kepox),(kmap( 63),    keth),
     &   (kmap( 64),   ketha),(kmap( 65),   kethy),(kmap( 66),   ketoh),
     &   (kmap( 67),   kfacd),(kmap( 68),   kform),(kmap( 69),    kgly),
     &   (kmap( 70),   kglyd),(kmap( 71),   kh2o2),(kmap( 72),    khcl),
     &   (kmap( 73),     khi),(kmap( 74),   khio3),(kmap( 75),   khno3),
     &   (kmap( 76),    khoi),(kmap( 77),   khono),(kmap( 78),   khpld),
     &   (kmap( 79),     ki2),(kmap( 80),   ki2o2),(kmap( 81),    kibr),
     &   (kmap( 82),    kicl),(kmap( 83),   kino2),(kmap( 84),   kino3),
     &   (kmap( 85),   kintr),(kmap( 86),   kiole),(kmap( 87),   kisop),
     &   (kmap( 88),   kispd),(kmap( 89),   kispx),(kmap( 90),   kixoy),
     &   (kmap( 91),    kket),(kmap( 92),   kmeoh),(kmap( 93),   kmepx),
     &   (kmap( 94),   kmgly),(kmap( 95),   kn2o5),(kmap( 96),    knh3),
     &   (kmap( 97),     kno),(kmap( 98),    kno2),(kmap( 99),    kno3),
     &   (kmap(100),   kntr1),(kmap(101),   kntr2),(kmap(102),     ko3),
     &   (kmap(103),    kole),(kmap(104),   kopan),(kmap(105),   kopen),
     &   (kmap(106),   kpacd),(kmap(107),    kpan),(kmap(108),   kpanx),
     &   (kmap(109),    kpar),(kmap(110),    kpna),(kmap(111),   kprpa),
     &   (kmap(112),   krooh),(kmap(113),    kso2),(kmap(114),    ksqt),
     &   (kmap(115),   ksulf),(kmap(116),   kterp),(kmap(117),    ktol),
     &   (kmap(118),   kxopn),(kmap(119),    kxyl),(kmap(120),   kccrs),
     &   (kmap(121),   kacg1),(kmap(122),   kacg2),(kmap(123),   kbcg1),
     &   (kmap(124),   kbcg2),(kmap(125),   kcprm),(kmap(126),   kfcrs),
     &   (kmap(127),   kfprm),(kmap(128),   kivoa),(kmap(129),   ksvoa),
     &   (kmap(130),     kna),(kmap(131),    kpal),(kmap(132),    kpca),
     &   (kmap(133),    kpcl),(kmap(134),    kpec),(kmap(135),    kpfe),
     &   (kmap(136),   kph2o),(kmap(137),     kpk),(kmap(138),    kpmg),
     &   (kmap(139),    kpmn),(kmap(140),   kpnh4),(kmap(141),   kpno3),
     &   (kmap(142),    khoa),(kmap(143),    kpsi),(kmap(144),   kpso4),
     &   (kmap(145),    kpti),(kmap(146),   kaoa1),(kmap(147),   kaoa2),
     &   (kmap(148),   kboa1),(kmap(149),   kboa2),(kmap(150),   kaoa0),
     &   (kmap(151),   kboa0),(kmap(152),    kbr2),(kmap(153),   kbrn2),
     &   (kmap(154),   kbrn3),(kmap(155),   kch3i),(kmap(156),    kcl2),
     &   (kmap(157),   kcln2),(kmap(158),   kcln3),(kmap(159),   kfmbr),
     &   (kmap(160),   kfmcl),(kmap(161),    khbr),(kmap(162),   khobr),
     &   (kmap(163),   khocl),(kmap(164),   kialk),(kmap(165),    kmb2),
     &   (kmap(166),   kmb2c),(kmap(167),    kmb3),(kmap(168),    kmbc),
     &   (kmap(169),   kmbc2),(kmap(170),    kmi2),(kmap(171),    kmib),
     &   (kmap(172),    kmic),(kmap(173),   kssbr),(kmap(174),   ksscl),
     &   (kmap(175),   kssn3),(kmap(176),    khg0),(kmap(177),    khg2),
     &   (kmap(178),  khgiip),(kmap(179), khgiipc),(kmap(180),    khgp),
     &   (kmap(181),   kivob),(kmap(182),   kivod),(kmap(183),   kivog),
     &   (kmap(184),   kpap0),(kmap(185),   kpap1),(kmap(186),   kpap2),
     &   (kmap(187),   kpap3),(kmap(188),   kpap4),(kmap(189),   kpas0),
     &   (kmap(190),   kpas1),(kmap(191),   kpas2),(kmap(192),   kpas3),
     &   (kmap(193),   kpas4),(kmap(194),   kpbs0),(kmap(195),   kpbs1),
     &   (kmap(196),   kpbs2),(kmap(197),   kpbs3),(kmap(198),   kpbs4),
     &   (kmap(199),   kpcp0),(kmap(200),   kpcp1),(kmap(201),   kpcp2),
     &   (kmap(202),   kpcp3),(kmap(203),   kpcp4),(kmap(204),   kpfp0),
     &   (kmap(205),   kpfp1),(kmap(206),   kpfp2),(kmap(207),   kpfp3),
     &   (kmap(208),   kpfp4),(kmap(209),   kvap1),(kmap(210),   kvap2),
     &   (kmap(211),   kvap3),(kmap(212),   kvap4),(kmap(213),   kvas1),
     &   (kmap(214),   kvas2),(kmap(215),   kvas3),(kmap(216),   kvas4),
     &   (kmap(217),   kvbs1),(kmap(218),   kvbs2),(kmap(219),   kvbs3),
     &   (kmap(220),   kvbs4),(kmap(221),   kvcp1),(kmap(222),   kvcp2),
     &   (kmap(223),   kvcp3),(kmap(224),   kvcp4),(kmap(225),   kvfp1),
     &   (kmap(226),   kvfp2),(kmap(227),   kvfp3),(kmap(228),   kvfp4),
     &   (kmap(229),   kaco3),(kmap(230),   kactp),(kmap(231),   kadcn),
     &   (kmap(232),   kaddc),(kmap(233),   kapip),(kmap(234),   kbal1),
     &   (kmap(235),   kbal2),(kmap(236),   kbalp),(kmap(237),   kbenp),
     &   (kmap(238),    kcho),(kmap(239),   ketep),(kmap(240),   kethp),
     &   (kmap(241),   khc3p),(kmap(242),   khc5p),(kmap(243),   khc8p),
     &   (kmap(244),   kketp),(kmap(245),   klimp),(kmap(246),   kmacp),
     &   (kmap(247),    kmcp),(kmap(248),   kmcto),(kmap(249),   kmctp),
     &   (kmap(250),   kmekp),(kmap(251),    kmo2),(kmap(252),   kmvkp),
     &   (kmap(253),   kolip),(kmap(254),   kolnd),(kmap(255),   kolnn),
     &   (kmap(256),   koltp),(kmap(257),   korap),(kmap(258),   kper1),
     &   (kmap(259),   kper2),(kmap(260),   ktlp1),(kmap(261),   ktolp),
     &   (kmap(262),    ktr2),(kmap(263),   kualp),(kmap(264),    kxy2),
     &   (kmap(265),   kxyl1),(kmap(266),   kxylp),(kmap(267),   kxyo2),
     &   (kmap(268),   kxyop),(kmap(269),    kacd),(kmap(270),    kace)
      equivalence
     &   (kmap(271),    kact),(kmap(272),    kald),(kmap(273),    kapi),
     &   (kmap(274),    kben),(kmap(275),    kcsl),(kmap(276),   kdcb1),
     &   (kmap(277),   kdcb2),(kmap(278),   kdcb3),(kmap(279),   kdien),
     &   (kmap(280),    keoh),(kmap(281),    kepx),(kmap(282),    kete),
     &   (kmap(283),   keteg),(kmap(284),    khc3),(kmap(285),    khc5),
     &   (kmap(286),    khc8),(kmap(287),   khno4),(kmap(288),   kishp),
     &   (kmap(289),    kiso),(kmap(290),   kison),(kmap(291),    klim),
     &   (kmap(292),   kmahp),(kmap(293),    kmct),(kmap(294),    kmoh),
     &   (kmap(295),   knald),(kmap(296),    koli),(kmap(297),    kolt),
     &   (kmap(298),   konit),(kmap(299),    kop1),(kmap(300),    kop2),
     &   (kmap(301),   kora1),(kmap(302),   kora2),(kmap(303),    kpaa),
     &   (kmap(304),   kphen),(kmap(305),    kppn),(kmap(306),   kuald),
     &   (kmap(307),    kxym),(kmap(308),    kxyo),(kmap(309),    kxyp),
     &   (kmap(310),   kbzc3),(kmap(311),    kbzo),(kmap(312),   kmac3),
     &   (kmap(313),   kmco3),(kmap(314),    ko3p),(kmap(315),   krco3),
     &   (kmap(316),   kro2c),(kmap(317),   kro2x),(kmap(318),   ktbuo),
     &   (kmap(319),   kxace),(kmap(320),   kxacr),(kmap(321),   kxaf1),
     &   (kmap(322),   kxaf2),(kmap(323),   kxaf3),(kmap(324),   kxbac),
     &   (kmap(325),   kxbal),(kmap(326),   kxcch),(kmap(327),    kxco),
     &   (kmap(328),   kxgld),(kmap(329),   kxgly),(kmap(330),   kxhch),
     &   (kmap(331),   kxho2),(kmap(332),   kxipr),(kmap(333),   kxma3),
     &   (kmap(334),   kxmac),(kmap(335),   kxmc3),(kmap(336),   kxmek),
     &   (kmap(337),   kxmeo),(kmap(338),   kxmgl),(kmap(339),   kxmvk),
     &   (kmap(340),   kxno2),(kmap(341),    kxoh),(kmap(342),   kxpd2),
     &   (kmap(343),   kxrc3),(kmap(344),   kxrch),(kmap(345),   kxrn3),
     &   (kmap(346),   kxtbu),(kmap(347),   ky6px),(kmap(348),   kyapx),
     &   (kmap(349),   kyrpx),(kmap(350),   kzrn3),(kmap(351),   kacro),
     &   (kmap(352),   kacye),(kmap(353),   kafg1),(kmap(354),   kafg2),
     &   (kmap(355),   kafg3),(kmap(356),   kalk1),(kmap(357),   kalk2),
     &   (kmap(358),   kalk3),(kmap(359),   kalk4),(kmap(360),   kalk5),
     &   (kmap(361),   kapin),(kmap(362),   karo1),(kmap(363),   karo2),
     &   (kmap(364),   kb124),(kmap(365),   kbacl),(kmap(366),   kbald),
     &   (kmap(367),   kbd13),(kmap(368),   kccho),(kmap(369),   kco3h),
     &   (kmap(370),   kcooh),(kmap(371),   kethe),(kmap(372),   khcho),
     &   (kmap(373),   kiprd),(kmap(374),   kmacr),(kmap(375),    kmek),
     &   (kmap(376),   kmpan),(kmap(377),    kmvk),(kmap(378),   kmxyl),
     &   (kmap(379),   knphe),(kmap(380),   kole1),(kmap(381),   kole2),
     &   (kmap(382),   koxyl),(kmap(383),   kpan2),(kmap(384),   kpbzn),
     &   (kmap(385),   kprd2),(kmap(386),   kprpe),(kmap(387),   kpxyl),
     &   (kmap(388),   kr6px),(kmap(389),   krapx),(kmap(390),   krcho),
     &   (kmap(391),   krno3),(kmap(392),   kro3h),(kmap(393),   ksesq),
     &   (kmap(394),   ktolu),(kmap(395),     kxn),(kmap(396),   ktr01),
     &   (kmap(397),   ktr02),(kmap(398),   ktr03),(kmap(399),   ktr04),
     &   (kmap(400),   ktr05),(kmap(401),   ktr06),(kmap(402),   ktr07),
     &   (kmap(403),   ktr08),(kmap(404),   ktr09)
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henso20  -- Henry's Law constant at STP for SO2 (molar/atm)
c     tfactso2 -- Temperature dependence of SO2 Henry's Law constant (1/K)
c     nbin     -- Number of aerosol size bins
c     cwmin    -- Minimum cloud water threshold (g/m3)
c     tamin    -- Cloud water freezing threshold (K)
c
c-----------------------------------------------------------------------
c
      real henso20
      real tfactso2
      real cwmin
      real tamin
c
      common /depchm/ henso20, tfactso2, cwmin, tamin
c
      integer nbin
c
      common /aerochm/ nbin
c
c-----------------------------------------------------------------------
c     Pointers to photolysis reactions that receive T,P adjustments
c
      integer jno2rxn,jo3rxn,jhcho1rxn,jhcho2rxn,jch3chorxn
      common /phtrxn/ jno2rxn,jo3rxn,jhcho1rxn,jhcho2rxn,jch3chorxn
c
c-----------------------------------------------------------------------
c     Pointers to heterogeneous reactions to be adjusted in KHETERO
c                    N2O5 + H2O -> 2 HNO3       (1)
c                    N2O5 + HCL -> CLN2 + HNO3  (2)
c                    NTR2 + H2O -> HNO3         (3)
c                    INTR + H2O -> HNO3         (4)
c                    SO2        -> SULF         (5)
c                    GLY        -> CGLY         (6)
c                   MGLY        -> CGLY         (7)
c
      integer, parameter :: NHETRXN = 7
      integer, dimension(NHETRXN) :: ihetrxn
      common /hetrxn/ ihetrxn
c
c-----------------------------------------------------------------------
c     Pointer used to lookup pig chemistry rate constant
c
c     ipigrxn  -- pointer to NO self-reaction
c                    NO + NO -> 2 NO2
c
      integer ipigrxn
c
      common /pigrxn/ ipigrxn

      integer, parameter :: l_O    = 1,
     &                      l_O1D  = 2,
     &                      l_OH   = 3,
     &                      l_HO2  = 4,
     &                      l_NO3  = 5,
     &                      l_N2O5 = 6,
     &                      l_NO   = 7,
     &                      l_NO2  = 8,
     &                      l_O3   = 9,
     &                      l_HONO = 10,
     &                      l_HNO3 = 11,
     &                      l_CO   = 12,
     &                      l_FORM = 13,
     &                      l_SO2  = 14,
     &                      l_SULF = 15
      integer, parameter :: NRPIGCHEM = 23
      integer, parameter :: NSPIGCHEM = 15
      integer, dimension(NSPIGCHEM) :: lsmap
      integer, dimension(NRPIGCHEM) :: lrmap
      common /pigchemmap/ lsmap, lrmap
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines
c
c     dtaero     -- user input time interval (min) between calls to aerosol 
c                   routines
c     dt_aero    -- adjusted (for I/O frequency) time interval (min) between
c                   calls to aerosol routines
c
      real dtaero
      real dt_aero
c
      common /aero_t/ dtaero, dt_aero
