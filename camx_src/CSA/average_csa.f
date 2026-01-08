      subroutine average_csa(icsa_grid,m1,m2,m3,igrid_mod,dt,
     &                   ncol,nrow,nlay,nspc,
     &                   tempk,press,conc,avcnc,ipa_cel)
      use procan
      implicit none
c
c----CAMx v7.32 250801
c
c     AVERAGE computes time-averaged concentrations. This version is for CSA arrays.
c
c     Copyright 1996 - 2025
c     Ramboll
c
c     Modifications: 
c
c     Input arguments:
c        icsa_grid          grid number for CSA domain
c        m1                 number of columns in met arrays
c        m2                 number of rows in met arrays
c        m3                 number of layers in met arrays
c        igrid_mod          index of modeling grid
c        dt                 time step for present grid concentration (s)
c        ncol               number of columns in conc arrays 
c        nrow               number of rows in conc arrays 
c        nlay               number of layers in conc arrays 
c        nspc               number of species in conc array
c        tempk              temperature field (K)
c        press              pressure field (mb)
c        conc               instant species concentration (umol/m3)
c        avcnc              average species concentration (gas=ppm)
c        ipa_cel            index of grid cell into CSA cells
c     Output arguments:
c        avcnc              average species concentration (gas=ppm)
c
c     Routines Called:
c        none
c
c     Called by:
c        CSA_AVRG
c
      include "camx.prm"
      include "camx.inc"
 
      integer icsa_grid
      integer m1
      integer m2
      integer m3
      integer igrid_mod
      integer ncol
      integer nrow
      integer nlay
      real    dt
      integer nspc
      real    tempk(m1,m2,m3),press(m1,m2,m3)
      real    avcnc(ncol,nrow,nlay,nspc)
      real    conc(ncol,nrow,nlay,nspc)
      integer ipa_cel(m1,m2,m3)
c
c-----Local variables ---
c
      integer lsp, i, j, k, icell, jcell, kcell, idxpa
      real    dtfact, convfac, tmp
c
c-----Entry point
c
c-----Increment running average
c
      dtfact = dt/(dtout*60.)
      do lsp = 1,nspc
        convfac = 1.
        do j = 2,m2-1
          do i = 2,m1-1
            do k=1,m3
                idxpa = ipa_cel(i,j,k)
                if( idxpa .GT. 0 ) then
                   if( ipanst(idxpa) .EQ. igrid_mod .AND. 
     &                     ipadom(idxpa) .EQ. icsa_grid ) then
                      icell = ipax(idxpa) - i_sw(ipadom(idxpa)) + 1
                      jcell = ipay(idxpa) - j_sw(ipadom(idxpa)) + 1
                      kcell = ipaz(idxpa) - b_lay(ipadom(idxpa)) + 1
                      tmp = 273./tempk(i,j,k)*press(i,j,k)/1013.
                      convfac = 1./(densfac*tmp)
                      avcnc(icell,jcell,kcell,lsp) = 
     &                      convfac*conc(icell,jcell,kcell,lsp)*dtfact +
     &                                      avcnc(icell,jcell,kcell,lsp)
                   endif
                endif
            enddo
          enddo
        enddo
      enddo
c
      return
      end
