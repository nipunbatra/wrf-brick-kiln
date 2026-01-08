      subroutine hr_xo7(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_XO7 solves halogen oxides by EBI iterations
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
      integer iter, nit
      real TINY
      parameter (TINY = 1.0E-20)
      real gainI, lossI
      real gainIO, lossIO
      real gainOIO, lossOIO
      real r215, r216, r217, r218, r219, r220, r221, r222, r223, r224
      real r225, r226, r227, r229
c
c --- Entry Point
c
      N2 = M - O2
      r215 = r(215)
      r216 = r(216)
      r227 = r(227)
      r229 = r(229)
c
      nit = 2
c
      do iter = 1,nit
c
c --- Update reaction rates that change with iteration of y1
c
      r217 = rk(217)*y1(lI)*yh(lO3)
      r218 = rk(218)*y1(lIO)
      r219 = rk(219)*y1(lIO)*y1(lIO)
      r220 = rk(220)*y1(lIO)*yh(lHO2)
      r221 = rk(221)*y1(lIO)*yh(lNO)
      r222 = rk(222)*y1(lIO)*yh(lNO2)
      r223 = rk(223)*y1(lOIO)
      r224 = rk(224)*y1(lOIO)*yh(lOH)
      r225 = rk(225)*y1(lOIO)*y1(lIO)
      r226 = rk(226)*y1(lOIO)*yh(lNO)
c
c --- Gain and loss terms for the XO species
c
      gainI = 0.0
     &       + ( 2.000)*r215
     &       +          r216
     &       +          r218
     &       + ( 0.400)*r219
     &       +          r221
     &       +          r223
     &       +          r227
     &       +          r229
c
      lossI = 0.0
     &       +          r217
c
      gainIO = 0.0
     &       +          r217
     &       +          r226
c
      lossIO = 0.0
     &       +          r218
     &       + ( 2.000)*r219
     &       +          r220
     &       +          r221
     &       +          r222
     &       +          r225
c
      gainOIO = 0.0
     &       + ( 0.400)*r219
     &       +          r227
c
      lossOIO = 0.0
     &       +          r223
     &       +          r224
     &       +          r225
     &       +          r226
c
c --- EBI solution for the XO species
c
      y1(lI) = MIN(1.0, MAX(TINY,(y0(lI) + gainI*dt)
     &                         / (1.0 + lossI*dt/y1(lI))))
c
      y1(lIO) = MIN(1.0, MAX(TINY,(y0(lIO) + gainIO*dt)
     &                         / (1.0 + lossIO*dt/y1(lIO))))
c
      y1(lOIO) = MIN(1.0, MAX(TINY,(y0(lOIO) + gainOIO*dt)
     &                         / (1.0 + lossOIO*dt/y1(lOIO))))
c
      enddo
c
      return
      end

