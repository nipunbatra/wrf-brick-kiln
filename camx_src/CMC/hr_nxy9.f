      subroutine hr_nxy9(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_NXY9 solves NO3 and N2O5 using Hertel's equations
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
      real self, fwd, bck
      real newNO3, lsNO3, lsN2O5, A, B, C, Q
c
c --- Entry Point
c
      N2 = M - O2
c
c
c --- New NO3 excluding from N2O5
c
      newNO3 = dt*(0.0
     &       + ( 0.200)*r(  9)
     &       +          r( 33)
     &       +          r( 37)
     &       +          r( 55)
     &       +          r( 57)
     &       +          r(120)
     &       +          r(121) )
c
c --- Loss of NO3 excluding self-reaction(s)
c
      lsNO3 = 1.0 + dt*(0.0
     &       +          rk(  5)
     &       +          rk(  6)
     &       +          rk( 58)*yh(lOH)
     &       +          rk( 59)*yh(lHO2)
     &       +          rk( 60)*yh(lNO)
     &       +          rk( 61)*yh(lNO2)
     &       +          rk( 63)*yh(lNO2)
     &       +          rk(141)*yh(lETE)
     &       +          rk(142)*yh(lOLT)
     &       +          rk(143)*yh(lOLI)
     &       +          rk(144)*yh(lDIEN)
     &       +          rk(145)*yh(lISO)
     &       +          rk(146)*yh(lAPI)
     &       +          rk(147)*yh(lLIM)
     &       +          rk(148)*yh(lHCHO)
     &       +          rk(149)*yh(lACD)
     &       +          rk(150)*yh(lALD)
     &       +          rk(151)*yh(lMACR)
     &       +          rk(152)*yh(lUALD)
     &       +          rk(153)*yh(lGLY)
     &       +          rk(154)*yh(lMGLY)
     &       +          rk(155)*yh(lPHEN)
     &       +          rk(156)*yh(lCSL)
     &       +          rk(157)*yh(lEPX)
     &       +          rk(158)*yh(lMCT)
     &       +          rk(159)*yh(lMPAN)
     &       +          rk(322)*yh(lMO2)
     &       +          rk(323)*yh(lETHP)
     &       +          rk(324)*yh(lHC3P)
     &       +          rk(325)*yh(lHC5P)
     &       +          rk(326)*yh(lHC8P)
     &       +          rk(327)*yh(lETEP)
     &       +          rk(328)*yh(lOLTP)
     &       +          rk(329)*yh(lOLIP)
     &       +          rk(330)*yh(lBENP)
     &       +          rk(331)*yh(lTLP1)
     &       +          rk(332)*yh(lTOLP)
     &       +          rk(333)*yh(lPER1)
     &       +          rk(334)*yh(lXYL1)
     &       +          rk(335)*yh(lXYLP)
     &       +          rk(336)*yh(lPER2)
     &       +          rk(337)*yh(lXYOP)
     &       +          rk(338)*yh(lISOP)
     &       +          rk(339)*yh(lAPIP)
     &       +          rk(340)*yh(lLIMP)
     &       +          rk(341)*yh(lACO3)
     &       +          rk(342)*yh(lRCO3)
     &       +          rk(343)*yh(lACTP)
     &       +          rk(344)*yh(lMEKP)
     &       +          rk(345)*yh(lKETP)
     &       +          rk(346)*yh(lMACP)
     &       +          rk(347)*yh(lMCP)
     &       +          rk(348)*yh(lMVKP)
     &       +          rk(349)*yh(lUALP)
     &       +          rk(350)*yh(lBALP)
     &       +          rk(351)*yh(lBAL1)
     &       +          rk(352)*yh(lADDC)
     &       +          rk(353)*yh(lMCTP)
     &       +          rk(354)*yh(lORAP)
     &       +          rk(355)*yh(lOLNN)
     &       +          rk(356)*yh(lOLND)
     &       +          rk(357)*yh(lADCN)
     &       +          rk(361)*yh(lXO2)
     &       +          rk(368)*yh(lDMS) )
c
c --- Loss of N2O5
c
      lsN2O5  = 1.0 + dt*(0.0
     &       +          rk( 64)
     &       +          rk( 65)*H2O )
c
c --- Forward reactions of NO3 to N2O5
c
      fwd  = dt*(0.0
     &       +          rk( 63)*yh(lNO2) )
c
c --- Backward reactions of N2O5 to NO3
c
      bck  = dt*(0.0
     &       +          rk( 64) )
c
c --- NO3 self-reaction excluding N2O5 production
c
      self = dt*(0.0
     &       + ( 2.000)*rk( 62) )
c
c --- Solve for NO3
c
      A = self*lsN2O5
      B = (lsN2O5*lsNO3) - (fwd*bck)
      C = lsN2O5*( y0(lNO3) + newNO3 ) + bck*y0(lN2O5)
      Q = -0.5 * ( B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C) )
c
c --- Update Concentrations
c
      y1(lNO3)  = MAX(1.0E-18, MAX(Q/A, -C/Q) )
      y1(lN2O5) = MAX(1.0E-15, (y0(lN2O5) + fwd*y1(lNO3)) / lsN2O5 )
c
      return
      end

