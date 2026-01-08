      subroutine hr_nox9(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_NOX9 solves the NOx family using Hertel's equations
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
      real jO3_O1D, rNO2_O3P, rO3P_O3, rO3wNO
      real totO1D, totNO, totNO2, totO3P, totO3
      real lsO1D, lsNO, lsNO2, lsO3P, lsO3
      real yldO3P, rNO_NO2, rNO2_NO
      real t1, t2, t3, t4, t5, f1, f2, f3, f4, s1, s2, A, B, C, Q
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
c --- Yield of O3P from O1D
c
      yldO3P  = ( rk( 40)*O2 + rk( 41)*N2 ) / lsO1D
c
c --- Conversion of NO2 to NO
c     (rNO2_NO = rate of NO2 to NO)
c
      rNO2_NO  = dt*(
     &       +          rk(  1)
     &       +          rk( 54)*yh(lO) )
c
c --- NO production terms except rNO2_NO
c     (totNO = initial NO + new production)
c
      totNO = y0(lNO) + dt*(
     &       +          r(  5)
     &       +          r(  7)
     &       +          r( 61) )
c
c --- Net loss of NO
c     (net either NO or NO2 produced)
c
      lsNO = 1.0 + dt*(
     &       +          rk( 49)*yh(lOH)
     &       +          rk( 51)*yh(lHO2)
     &       + ( 0.065)*rk(174)*yh(lHC3P)
     &       + ( 0.136)*rk(175)*yh(lHC5P)
     &       + ( 0.261)*rk(176)*yh(lHC8P)
     &       + ( 0.030)*rk(178)*yh(lOLTP)
     &       + ( 0.050)*rk(179)*yh(lOLIP)
     &       + ( 0.082)*rk(180)*yh(lBENP)
     &       + ( 0.050)*rk(182)*yh(lTOLP)
     &       + ( 0.050)*rk(183)*yh(lPER1)
     &       + ( 0.050)*rk(185)*yh(lXYLP)
     &       + ( 0.050)*rk(186)*yh(lPER2)
     &       + ( 0.050)*rk(187)*yh(lXYOP)
     &       + ( 0.120)*rk(188)*yh(lISOP)
     &       + ( 0.180)*rk(189)*yh(lAPIP) )
c
c --- Conversion of NO to NO2 except O3+NO (rO3wNO)
c
      rNO_NO2 = dt*(
     &       +          rk( 48)*yh(lO)
     &       +          rk( 50)*yh(lHO2)
     &       + ( 2.000)*rk( 52)*yh(lNO)*O2
     &       +          rk( 60)*yh(lNO3)
     &       +          rk(172)*yh(lMO2)
     &       +          rk(173)*yh(lETHP)
     &       + ( 0.935)*rk(174)*yh(lHC3P)
     &       + ( 0.864)*rk(175)*yh(lHC5P)
     &       + ( 0.739)*rk(176)*yh(lHC8P)
     &       +          rk(177)*yh(lETEP)
     &       + ( 0.970)*rk(178)*yh(lOLTP)
     &       + ( 0.950)*rk(179)*yh(lOLIP)
     &       + ( 0.918)*rk(180)*yh(lBENP)
     &       +          rk(181)*yh(lTLP1)
     &       + ( 0.950)*rk(182)*yh(lTOLP)
     &       + ( 0.950)*rk(183)*yh(lPER1)
     &       +          rk(184)*yh(lXYL1)
     &       + ( 0.950)*rk(185)*yh(lXYLP)
     &       + ( 0.950)*rk(186)*yh(lPER2)
     &       + ( 0.950)*rk(187)*yh(lXYOP)
     &       + ( 0.880)*rk(188)*yh(lISOP)
     &       + ( 0.820)*rk(189)*yh(lAPIP)
     &       +          rk(190)*yh(lLIMP)
     &       +          rk(191)*yh(lACO3)
     &       +          rk(192)*yh(lRCO3)
     &       +          rk(193)*yh(lACTP)
     &       +          rk(194)*yh(lMEKP)
     &       +          rk(195)*yh(lKETP)
     &       +          rk(196)*yh(lMACP)
     &       +          rk(197)*yh(lMCP)
     &       +          rk(198)*yh(lMVKP)
     &       +          rk(199)*yh(lUALP)
     &       +          rk(200)*yh(lBALP)
     &       +          rk(201)*yh(lBAL1)
     &       +          rk(202)*yh(lADDC)
     &       +          rk(203)*yh(lMCTP)
     &       +          rk(204)*yh(lORAP)
     &       +          rk(205)*yh(lOLNN)
     &       +          rk(206)*yh(lOLND)
     &       +          rk(207)*yh(lADCN)
     &       +          rk(208)*yh(lXO2) )
c
c --- Remaining NO2 production
c
      totNO2 = y0(lNO2) + dt*(
     &       +          r(  6)
     &       +          r(  8)
     &       + ( 0.800)*r(  9)
     &       +          r( 31)
     &       +          r( 32)
     &       +          r( 53)
     &       +          r( 58)
     &       + ( 0.700)*r( 59)
     &       +          r( 60)
     &       + ( 2.000)*r( 62)
     &       +          r( 64)
     &       +          r( 67)
     &       +          r( 68)
     &       +          r(122)
     &       +          r(123)
     &       +          r(124)
     &       + ( 0.680)*r(151)
     &       + ( 0.500)*r(157)
     &       +          r(159)
     &       +          r(167)
     &       +          r(169)
     &       +          r(171)
     &       +          r(206)
     &       +          r(207)
     &       +          r(273)
     &       + ( 0.500)*r(282)
     &       + ( 0.700)*r(283)
     &       +          r(309)
     &       +          r(318)
     &       + ( 0.700)*r(319)
     &       +          r(322)
     &       +          r(323)
     &       +          r(324)
     &       +          r(325)
     &       +          r(326)
     &       +          r(327)
     &       +          r(328)
     &       +          r(329)
     &       +          r(330)
     &       +          r(331)
     &       +          r(332)
     &       +          r(333)
     &       +          r(334)
     &       +          r(335)
     &       +          r(336)
     &       +          r(337)
     &       +          r(338)
     &       +          r(339)
     &       +          r(340)
     &       +          r(341)
     &       +          r(342)
     &       +          r(343)
     &       +          r(344)
     &       +          r(345)
     &       +          r(346)
     &       +          r(347)
     &       +          r(348)
     &       +          r(349)
     &       +          r(350)
     &       +          r(351)
     &       +          r(352)
     &       +          r(353)
     &       +          r(354)
     &       +          r(355)
     &       + ( 2.000)*r(356)
     &       + ( 2.000)*r(357)
     &       + ( 0.500)*r(359)
     &       +          r(360)
     &       +          r(361) )
c
c --- Net loss of NO2
c     (net either NO or NO2 produced)
c
      lsNO2 = 1.0 + dt*(
     &       +          rk( 37)*yh(lO3)
     &       +          rk( 55)*yh(lO)
     &       +          rk( 56)*yh(lOH)
     &       +          rk( 63)*yh(lNO3)
     &       +          rk( 66)*yh(lHO2)
     &       +          rk(166)*yh(lACO3)
     &       +          rk(168)*yh(lRCO3)
     &       +          rk(170)*yh(lMACP)
     &       +          rk(209)*yh(lBAL2)
     &       +          rk(210)*yh(lCHO)
     &       +          rk(211)*yh(lMCTO) )
c
c --- Production of O3P except NO2+hv
c
      totO3P = y0(lO) + dt*(
     &       +          r(  2)
     &       +          r(  3)*yldO3P
     &       +          r(  6)
     &       + ( 0.090)*r(129) )
c
c --- Net loss of O3P
c
      lsO3P = 1.0 + dt*(
     &       +          rk( 38)*O2*M
     &       +          rk( 39)*yh(lO3)
     &       +          rk( 48)*yh(lNO)
     &       +          rk( 54)*yh(lNO2)
     &       +          rk( 55)*yh(lNO2) )
c
c --- Production of O3 except O3P+O2
c
      totO3 = y0(lO3)
c
c --- Net loss of O3 except O3+NO (rO3wNO)
c
      lsO3 = 1.0 + dt*(
     &       +          rk(  2)
     &       +          rk(  3)
     &       +          rk( 34)*yh(lOH)
     &       +          rk( 35)*yh(lHO2)
     &       +          rk( 37)*yh(lNO2)
     &       +          rk( 39)*yh(lO)
     &       +          rk(126)*yh(lETE)
     &       +          rk(127)*yh(lOLT)
     &       +          rk(128)*yh(lOLI)
     &       +          rk(129)*yh(lDIEN)
     &       +          rk(130)*yh(lISO)
     &       +          rk(131)*yh(lAPI)
     &       +          rk(132)*yh(lLIM)
     &       +          rk(133)*yh(lMACR)
     &       +          rk(134)*yh(lMVK)
     &       +          rk(135)*yh(lUALD)
     &       +          rk(136)*yh(lDCB1)
     &       +          rk(137)*yh(lDCB2)
     &       +          rk(138)*yh(lDCB3)
     &       +          rk(139)*yh(lEPX)
     &       +          rk(140)*yh(lMCTO) )
c
c --- Specific reactions
c
      jO3_O1D  = rk(  3)
      rNO2_O3P = dt*rk(  1)
      rO3P_O3  = dt*rk( 38)*O2*M
      rO3wNO   = dt*rk( 36)
c
c --- Collect common terms
c
      t1 = rNO2_O3P / lsNO2
      t2 = rNO2_NO / lsNO2
      t3 = rNO_NO2 / lsNO
      t4 = rO3P_O3  / lsO3P
      t5 = t3*totNO - t2*totNO2
      f1 = 1.0 + t2 + t3
      f2 = t1*t4
      f3 = lsO3*lsNO + rO3wNO*totNO
      f4 = totO3 + totO3P * t4
c
c --- Solve for change in NO and NO2 
c
      A = rO3wNO * (f1 - f2)
      B = f1*f3 + rO3wNO*(f2*(totNO2 - totNO) + f4 + t5)
      C = rO3wNO*totNO*(f4 + totNO2 * f2) + f3*t5
      Q = -0.5 * (B + SIGN(1.0,B)*SQRT(B*B - 4.0*A*C))
      Q = MAX(Q/A ,C/Q)
c
c --- Update Concentrations
c
      y1(lNO)  = MAX(1.0E-15, (totNO + Q) / lsNO )
      y1(lNO2) = MAX(1.0E-15, (totNO2 - Q) / lsNO2 )
      s1 = totO3P + rNO2_O3P*y1(lNO2)
      s2 = t4*s1
      y1(lO3)  = MAX(1.0E-15, (totO3 + s2) / 
     &                             ( lsO3 + rO3wNO * y1(lNO) ) )
      y1(lO)   = MAX(1.0E-20, s1 / lsO3P )
      totO1D = jO3_O1D*y1(lO3)
      y1(lO1D) = MAX(1.0E-25, totO1D / lsO1D )
c
      return
      end

