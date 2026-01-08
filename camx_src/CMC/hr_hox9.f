      subroutine hr_hox9(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_HOX9 solves the HOx family using Hertel's equations
c
c     Copyright 1996 - 2025
c     Ramboll
c     Created by the CMC version 5.2.7
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        EBISOLV
c
c --- Argument definitions:
c        y0  - initial Conc for this step (ppm)
c        yh  - current Conc iteration (ppm)
c        y1  - next Conc iteration (ppm)
c        H2O - water vapor Conc (ppm)
c        M   - total gas Conc (ppm)
c        O2  - oxygen Conc (ppm)
c        CH4 - methane Conc (ppm)
c        H2  - hydrogen Conc (ppm)
c        ny  - dimension of y0, y1 and yh
c        rk  - rate constants (ppm-n hr-1)
c        r   - reaction rates (hr-1)
c        nr  - dimension of rk and r
c        dt  - time step (hr)
c
c --- Includes:
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
c --- Arguments:
      integer ny, nr
      real y0(ny+1), y1(ny+1), yh(ny+1)
      real rk(nr), r(nr)
      real H2O, M, O2, CH4, H2, dt
c
c --- Local variables:
      real N2
      real t1, t2, t3, A, B, C, Q
      real newOH, newHO2, newHONO, newPNA
      real lsOH, lsHO2, lsHONO, lsPNA, lsO1D, yldOH, self
      real rOH_HO2, rHO2_OH, rOH_HONO, rHONO_OH, rHO2_PNA, rPNA_HO2
c
c --- Entry Point
c
      N2 = M - O2
c
c --- O1D loss
c
      lsO1D = 
     &       + rk( 40)*O2
     &       + rk( 41)*N2
     &       + rk( 42)*H2O
c
c --- Yield of OH from O1D
c
      yldOH  = ( 
     &         ( 2.000)*rk( 42)*H2O ) / lsO1D
c
c --- Conversion of HONO to OH
c
      rHONO_OH  = dt*(
     &       +          rk(  7) )
c
c --- Conversion of HO2 to OH
c
      rHO2_OH  = dt*(
     &       +          rk( 35)*yh(lO3)
     &       +          rk( 50)*yh(lNO)
     &       + ( 0.700)*rk( 59)*yh(lNO3)
     &       + ( 0.440)*rk(231)*yh(lACO3)
     &       + ( 0.440)*rk(232)*yh(lRCO3)
     &       + ( 0.150)*rk(233)*yh(lACTP) )
c
c --- Other OH production terms
c
      newOH = y0(lOH) + dt*(
     &       + ( 2.000)*r(  4)
     &       +          r(  8)
     &       + ( 0.200)*r(  9)
     &       + ( 0.340)*r( 19)
     &       +          r( 28)
     &       +          r( 29)
     &       +          r( 30)
     &       + ( 0.080)*r(126)
     &       + ( 0.220)*r(127)
     &       + ( 0.460)*r(128)
     &       + ( 0.280)*r(129)
     &       + ( 0.250)*r(130)
     &       + ( 0.850)*r(131)
     &       + ( 0.850)*r(132)
     &       + ( 0.190)*r(133)
     &       + ( 0.160)*r(134)
     &       + ( 0.100)*r(135)
     &       + ( 0.050)*r(136)
     &       + ( 0.050)*r(137)
     &       + ( 0.050)*r(138)
     &       + ( 0.050)*r(139)
     &       + ( 0.500)*r(157)
     &       + ( 0.280)*r(160)
     &       + ( 0.490)*r(161)
     &       + ( 0.158)*r(162)
     &       + ( 0.390)*r(163)
     &       + ( 0.158)*r(164)
     &       + ( 0.390)*r(165)
     &       +          r(  3)*yldOH )
c
c --- Conversion of PNA to HO2
c
      rPNA_HO2  = dt*(
     &       + ( 0.800)*rk(  9)
     &       +          rk( 67) )
c
c --- Conversion of OH to HO2
c
      rOH_HO2  = dt*(
     &       +          rk( 34)*yh(lO3)
     &       +          rk( 43)*H2
     &       +          rk( 47)*yh(lH2O2)
     &       +          rk( 58)*yh(lNO3)
     &       +          rk( 69)*yh(lSO2)
     &       +          rk( 70)*yh(lCO)
     &       + ( 0.049)*rk( 75)*yh(lHC8)
     &       + ( 0.350)*rk( 80)*yh(lACE)
     &       + ( 0.648)*rk( 81)*yh(lBEN)
     &       + ( 0.177)*rk( 82)*yh(lTOL)
     &       + ( 0.177)*rk( 83)*yh(lXYM)
     &       + ( 0.177)*rk( 84)*yh(lXYP)
     &       + ( 0.177)*rk( 85)*yh(lXYO)
     &       +          rk( 89)*yh(lHCHO)
     &       +          rk( 95)*yh(lHKET)
     &       +          rk( 99)*yh(lGLY)
     &       + ( 0.520)*rk(101)*yh(lDCB1)
     &       + ( 0.520)*rk(102)*yh(lDCB2)
     &       + ( 0.560)*rk(103)*yh(lDCB3)
     &       + ( 0.730)*rk(105)*yh(lPHEN)
     &       + ( 0.730)*rk(106)*yh(lCSL)
     &       +          rk(107)*yh(lEPX)
     &       +          rk(109)*yh(lMOH)
     &       +          rk(110)*yh(lEOH)
     &       +          rk(111)*yh(lROH)
     &       +          rk(112)*yh(lETEG)
     &       +          rk(117)*yh(lORA1) )
c
c --- Other HO2 production terms
c
      newHO2 = y0(lHO2) + dt*(
     &       + ( 2.000)*r( 11)
     &       +          r( 12)
     &       +          r( 13)
     &       + ( 1.220)*r( 15)
     &       +          r( 18)
     &       + ( 0.660)*r( 19)
     &       + ( 2.000)*r( 23)
     &       +          r( 24)
     &       + ( 1.500)*r( 25)
     &       + ( 1.500)*r( 26)
     &       +          r( 27)
     &       +          r( 28)
     &       +          r( 29)
     &       +          r( 31)
     &       + ( 0.150)*r(126)
     &       + ( 0.320)*r(127)
     &       + ( 0.070)*r(128)
     &       + ( 0.300)*r(129)
     &       + ( 0.250)*r(130)
     &       + ( 0.100)*r(131)
     &       + ( 0.100)*r(132)
     &       + ( 0.140)*r(133)
     &       + ( 0.110)*r(134)
     &       + ( 0.072)*r(135)
     &       +          r(136)
     &       +          r(137)
     &       +          r(138)
     &       + ( 1.500)*r(139)
     &       +          r(148)
     &       +          r(152)
     &       +          r(153)
     &       + ( 1.500)*r(157)
     &       + ( 0.290)*r(160)
     &       + ( 0.010)*r(161)
     &       + ( 0.308)*r(162)
     &       + ( 0.010)*r(163)
     &       + ( 0.308)*r(164)
     &       + ( 0.010)*r(165)
     &       +          r(172)
     &       +          r(173)
     &       + ( 0.660)*r(174)
     &       + ( 0.200)*r(175)
     &       + ( 0.606)*r(176)
     &       +          r(177)
     &       + ( 0.780)*r(178)
     &       + ( 0.830)*r(179)
     &       + ( 0.918)*r(180)
     &       + ( 0.950)*r(182)
     &       + ( 0.500)*r(183)
     &       + ( 0.950)*r(185)
     &       + ( 0.950)*r(186)
     &       + ( 0.950)*r(187)
     &       + ( 0.880)*r(188)
     &       + ( 0.820)*r(189)
     &       +          r(190)
     &       + ( 0.670)*r(194)
     &       + ( 0.770)*r(195)
     &       + ( 0.500)*r(197)
     &       + ( 0.300)*r(198)
     &       +          r(199)
     &       +          r(202)
     &       +          r(204)
     &       +          r(205)
     &       + ( 0.740)*r(248)
     &       +          r(249)
     &       + ( 0.894)*r(250)
     &       + ( 0.842)*r(251)
     &       + ( 0.910)*r(252)
     &       +          r(253)
     &       +          r(254)
     &       +          r(255)
     &       + ( 1.600)*r(256)
     &       +          r(257)
     &       + ( 2.000)*r(258)
     &       + ( 2.000)*r(259)
     &       +          r(260)
     &       + ( 2.000)*r(261)
     &       + ( 2.000)*r(262)
     &       + ( 2.000)*r(263)
     &       +          r(264)
     &       +          r(265)
     &       +          r(266)
     &       + ( 0.900)*r(267)
     &       + ( 0.900)*r(268)
     &       + ( 0.500)*r(269)
     &       + ( 0.834)*r(270)
     &       +          r(271)
     &       + ( 0.500)*r(272)
     &       +          r(273)
     &       +          r(274)
     &       +          r(275)
     &       +          r(276)
     &       +          r(277)
     &       + ( 2.000)*r(278)
     &       +          r(279)
     &       +          r(280)
     &       + ( 2.000)*r(281)
     &       + ( 0.500)*r(282)
     &       +          r(283)
     &       +          r(284)
     &       + ( 0.500)*r(285)
     &       + ( 0.394)*r(286)
     &       + ( 0.342)*r(287)
     &       + ( 0.303)*r(288)
     &       + ( 0.500)*r(289)
     &       + ( 0.500)*r(290)
     &       + ( 0.500)*r(291)
     &       + ( 0.600)*r(292)
     &       +          r(294)
     &       +          r(295)
     &       +          r(297)
     &       +          r(298)
     &       +          r(299)
     &       + ( 0.500)*r(300)
     &       + ( 0.500)*r(301)
     &       + ( 0.500)*r(302)
     &       + ( 0.330)*r(306)
     &       + ( 0.500)*r(307)
     &       + ( 0.500)*r(309)
     &       + ( 0.500)*r(310)
     &       + ( 0.500)*r(311)
     &       + ( 2.000)*r(314)
     &       +          r(315)
     &       +          r(317)
     &       +          r(319)
     &       +          r(322)
     &       +          r(323)
     &       + ( 0.254)*r(324)
     &       + ( 0.488)*r(325)
     &       + ( 0.820)*r(326)
     &       +          r(327)
     &       + ( 0.790)*r(328)
     &       + ( 0.860)*r(329)
     &       +          r(330)
     &       +          r(332)
     &       + ( 0.500)*r(333)
     &       +          r(335)
     &       +          r(336)
     &       +          r(337)
     &       +          r(338)
     &       +          r(339)
     &       +          r(340)
     &       + ( 0.670)*r(344)
     &       +          r(345)
     &       +          r(347)
     &       + ( 0.300)*r(348)
     &       +          r(349)
     &       +          r(352)
     &       +          r(354)
     &       +          r(355)
     &       +          r(358)
     &       + ( 0.500)*r(359) )
c
c
c --- Conversion of OH to HONO
c
      rOH_HONO  = dt*(
     &       +          rk( 49)*yh(lNO) )
c
c --- Other HONO production terms
c
      newHONO = y0(lHONO) 
c
c --- Conversion of HO2 to PNA
c
      rHO2_PNA  = dt*(
     &       +          rk( 66)*yh(lNO2) )
c
c --- Other PNA production terms
c
      newPNA = y0(lHNO4) 
c
c --- Net loss of OH
c
      lsOH = 1.0 + dt*(
     &       +          rk( 34)*yh(lO3)
     &       +          rk( 43)*H2
     &       +          rk( 44)*yh(lHO2)
     &       +          rk( 47)*yh(lH2O2)
     &       +          rk( 49)*yh(lNO)
     &       +          rk( 53)*yh(lHONO)
     &       +          rk( 56)*yh(lNO2)
     &       +          rk( 57)*yh(lHNO3)
     &       +          rk( 58)*yh(lNO3)
     &       +          rk( 68)*yh(lHNO4)
     &       +          rk( 69)*yh(lSO2)
     &       +          rk( 70)*yh(lCO)
     &       +          rk( 71)*CH4
     &       +          rk( 72)*yh(lETH)
     &       +          rk( 73)*yh(lHC3)
     &       +          rk( 74)*yh(lHC5)
     &       +          rk( 75)*yh(lHC8)
     &       +          rk( 76)*yh(lETE)
     &       +          rk( 77)*yh(lOLT)
     &       +          rk( 78)*yh(lOLI)
     &       +          rk( 79)*yh(lDIEN)
     &       + ( 0.350)*rk( 80)*yh(lACE)
     &       +          rk( 81)*yh(lBEN)
     &       +          rk( 82)*yh(lTOL)
     &       +          rk( 83)*yh(lXYM)
     &       +          rk( 84)*yh(lXYP)
     &       +          rk( 85)*yh(lXYO)
     &       +          rk( 86)*yh(lISO)
     &       +          rk( 87)*yh(lAPI)
     &       +          rk( 88)*yh(lLIM)
     &       +          rk( 89)*yh(lHCHO)
     &       +          rk( 90)*yh(lACD)
     &       +          rk( 91)*yh(lALD)
     &       +          rk( 92)*yh(lACT)
     &       +          rk( 93)*yh(lMEK)
     &       +          rk( 94)*yh(lKET)
     &       +          rk( 95)*yh(lHKET)
     &       +          rk( 96)*yh(lMACR)
     &       +          rk( 97)*yh(lMVK)
     &       +          rk( 98)*yh(lUALD)
     &       +          rk( 99)*yh(lGLY)
     &       +          rk(100)*yh(lMGLY)
     &       +          rk(101)*yh(lDCB1)
     &       +          rk(102)*yh(lDCB2)
     &       +          rk(103)*yh(lDCB3)
     &       +          rk(104)*yh(lBALD)
     &       +          rk(105)*yh(lPHEN)
     &       +          rk(106)*yh(lCSL)
     &       +          rk(107)*yh(lEPX)
     &       +          rk(108)*yh(lMCT)
     &       +          rk(109)*yh(lMOH)
     &       +          rk(110)*yh(lEOH)
     &       +          rk(111)*yh(lROH)
     &       +          rk(112)*yh(lETEG)
     &       + ( 0.650)*rk(113)*yh(lOP1)
     &       + ( 0.990)*rk(114)*yh(lOP2)
     &       +          rk(116)*yh(lMAHP)
     &       +          rk(117)*yh(lORA1)
     &       +          rk(118)*yh(lORA2)
     &       + ( 0.650)*rk(119)*yh(lPAA)
     &       +          rk(120)*yh(lPAN)
     &       +          rk(121)*yh(lPPN)
     &       +          rk(122)*yh(lMPAN)
     &       +          rk(123)*yh(lONIT)
     &       +          rk(124)*yh(lNALD)
     &       +          rk(125)*yh(lISON)
     &       +          rk(364)*yh(lACT)
     &       +          rk(365)*yh(lECH4)
     &       +          rk(366)*yh(lDMS)
     &       +          rk(367)*yh(lDMS)*O2 )
c
c --- Loss of HO2, excluding self-reaction
c     (net either HO2 or OH produced)
c
      lsHO2 = 1.0 + rHO2_OH + rHO2_PNA + dt*(
     &       +          rk( 44)*yh(lOH)
     &       +          rk( 51)*yh(lNO)
     &       + ( 0.300)*rk( 59)*yh(lNO3)
     &       +          rk(212)*yh(lMO2)
     &       +          rk(213)*yh(lETHP)
     &       +          rk(214)*yh(lHC3P)
     &       +          rk(215)*yh(lHC5P)
     &       +          rk(216)*yh(lHC8P)
     &       +          rk(217)*yh(lETEP)
     &       +          rk(218)*yh(lOLTP)
     &       +          rk(219)*yh(lOLIP)
     &       +          rk(220)*yh(lBENP)
     &       +          rk(221)*yh(lTLP1)
     &       +          rk(222)*yh(lTOLP)
     &       +          rk(223)*yh(lPER1)
     &       +          rk(224)*yh(lXYL1)
     &       +          rk(225)*yh(lXYLP)
     &       +          rk(226)*yh(lPER2)
     &       +          rk(227)*yh(lXYOP)
     &       +          rk(228)*yh(lISOP)
     &       +          rk(229)*yh(lAPIP)
     &       +          rk(230)*yh(lLIMP)
     &       + ( 0.560)*rk(231)*yh(lACO3)
     &       + ( 0.560)*rk(232)*yh(lRCO3)
     &       + ( 0.850)*rk(233)*yh(lACTP)
     &       +          rk(234)*yh(lMEKP)
     &       +          rk(235)*yh(lKETP)
     &       +          rk(236)*yh(lMACP)
     &       +          rk(237)*yh(lMCP)
     &       +          rk(238)*yh(lMVKP)
     &       +          rk(239)*yh(lUALP)
     &       +          rk(240)*yh(lADDC)
     &       +          rk(241)*yh(lCHO)
     &       +          rk(242)*yh(lMCTP)
     &       +          rk(243)*yh(lORAP)
     &       +          rk(244)*yh(lOLNN)
     &       +          rk(245)*yh(lOLND)
     &       +          rk(246)*yh(lADCN)
     &       +          rk(247)*yh(lXO2) )
c
c --- HO2 self-reaction
c
      self = dt*2.0*( 
     &       +          rk( 45)
     &       +          rk( 46)*H2O )
c
c --- Loss of HONO
c
      lsHONO = 1.0 + dt*(
     &       +          rk(  7)
     &       +          rk( 53)*yh(lOH) )
c
c --- Loss of PNA
c
      lsPNA = 1.0 + dt*(
     &       +          rk(  9)
     &       +          rk( 67)
     &       +          rk( 68)*yh(lOH) )
c
c --- Collect common terms
c
      t1 = 1.0 / ( lsOH*lsHONO - rHONO_OH*rOH_HONO )
      t2 = rOH_HO2*t1
      t3 = rPNA_HO2 / lsPNA
c
c --- Solve for HO2
c
      A = self
      B = lsHO2 - t3*rHO2_PNA - t2*rHO2_OH*lsHONO
      C = newHO2 + t3 * newPNA + t2*( newOH*lsHONO + newHONO*rHONO_OH )
      Q = -0.5 * (B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C))
c
c --- Update Concentrations
c
      y1(lHO2)  = MAX(1.0E-18, MAX(Q/A ,-C/Q) )
      y1(lOH)   = MAX(1.0E-18, ( ( newOH + rHO2_OH*y1(lHO2) )*lsHONO + 
     &                                        rHONO_OH*newHONO ) * t1 )
      y1(lHNO4)  = MAX(1.0E-15, ( newPNA + rHO2_PNA*y1(lHO2) ) / lsPNA )
      y1(lHONO) = MAX(1.0E-15, ( newHONO + rOH_HONO*y1(lOH) ) / lsHONO )
c
      return
      end

