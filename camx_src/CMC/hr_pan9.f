      subroutine hr_pan9(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_PAN9 solves PAN and RCO3 using Hertel's equations
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
      real fwd, bck, self, A, B, C, Q
      real newRCO3, lsPAN, lsRCO3
c
c --- Entry Point
c
      N2 = M - O2
c
c
c --- New RCO3 excluding from PAN
c
      newRCO3 = dt*(0.0
     &       +          r( 14)
     &       + ( 0.784)*r( 15)
     &       + ( 0.900)*r( 16)
     &       + ( 0.500)*r( 17)
     &       +          r( 18)
     &       + ( 0.670)*r( 19)
     &       +          r( 24)
     &       + ( 0.250)*r( 25)
     &       + ( 0.250)*r( 26)
     &       +          r( 90)
     &       + ( 0.313)*r( 98)
     &       +          r(100)
     &       + ( 0.650)*r(119)
     &       + ( 0.090)*r(128)
     &       + ( 0.150)*r(129)
     &       + ( 0.100)*r(130)
     &       + ( 0.100)*r(133)
     &       + ( 0.280)*r(134)
     &       + ( 0.002)*r(135)
     &       +          r(149)
     &       +          r(154)
     &       +          r(193)
     &       + ( 0.230)*r(195)
     &       + ( 0.350)*r(196)
     &       + ( 0.700)*r(198)
     &       + ( 0.150)*r(233)
     &       + ( 0.500)*r(269)
     &       + ( 0.269)*r(272)
     &       + ( 1.160)*r(274)
     &       + ( 0.160)*r(310)
     &       +          r(343)
     &       + ( 0.538)*r(346)
     &       + ( 0.700)*r(348) )
c
c --- Loss of RCO3 excluding self-reaction
c
      lsRCO3 = 1.0 + dt*(0.0
     &       +          rk(166)*yh(lNO2)
     &       +          rk(191)*yh(lNO)
     &       +          rk(231)*yh(lHO2)
     &       +          rk(267)*yh(lMO2)
     &       +          rk(285)*yh(lETHP)
     &       +          rk(286)*yh(lHC3P)
     &       +          rk(287)*yh(lHC5P)
     &       +          rk(288)*yh(lHC8P)
     &       +          rk(289)*yh(lETEP)
     &       +          rk(290)*yh(lOLTP)
     &       +          rk(291)*yh(lOLIP)
     &       +          rk(292)*yh(lBENP)
     &       +          rk(293)*yh(lTLP1)
     &       +          rk(294)*yh(lTOLP)
     &       +          rk(295)*yh(lPER1)
     &       +          rk(296)*yh(lXYL1)
     &       +          rk(297)*yh(lXYLP)
     &       +          rk(298)*yh(lPER2)
     &       +          rk(299)*yh(lXYOP)
     &       +          rk(300)*yh(lISOP)
     &       +          rk(301)*yh(lAPIP)
     &       +          rk(302)*yh(lLIMP)
     &       +          rk(304)*yh(lRCO3)
     &       + ( 0.500)*rk(305)*yh(lACTP)
     &       +          rk(306)*yh(lMEKP)
     &       +          rk(307)*yh(lKETP)
     &       + ( 0.731)*rk(308)*yh(lMACP)
     &       +          rk(309)*yh(lMCP)
     &       +          rk(311)*yh(lUALP)
     &       +          rk(312)*yh(lBALP)
     &       +          rk(313)*yh(lBAL1)
     &       +          rk(314)*yh(lADDC)
     &       +          rk(315)*yh(lMCTP)
     &       +          rk(316)*yh(lORAP)
     &       +          rk(317)*yh(lOLNN)
     &       +          rk(318)*yh(lOLND)
     &       +          rk(319)*yh(lADCN)
     &       +          rk(320)*yh(lXO2)
     &       +          rk(341)*yh(lNO3) )
c
c --- Loss of PAN
c
      lsPAN = 1.0 + dt*(0.0
     &       + rk( 32)
     &       + rk( 33)
     &       + rk(120)*yh(lOH)
     &       + rk(167) )
c
c
c --- Forward reaction of RCO3 to PAN
c
      fwd  = dt*(0.0
     &       +          rk(166)*yh(lNO2) )
c
c --- Backward reactions of PAN to RCO3
c
      bck  = dt*(0.0
     &       +          rk( 32)
     &       +          rk(167) )
c
c --- RCO3 self-reaction
c
      self = dt*(0.0
     &       + ( 2.000)*rk(303) )
c
c --- Solve for RCO3
c
      A = self*lsPAN
      B = ( lsPAN*lsRCO3) - (fwd*bck)
      C = lsPAN*( y0(lACO3) + newRCO3 ) + bck*y0(lPAN)
      Q = -0.5 * ( B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C) )
c
c --- Update Concentrations
c
      y1(lACO3) = MAX(1.0E-18, MAX(Q/A, -C/Q) )
      y1(lPAN)  = MAX(1.0E-15, (y0(lPAN) + fwd*y1(lACO3)) / lsPAN )
c
      return
      end

