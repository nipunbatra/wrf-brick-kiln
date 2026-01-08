      subroutine csa_avrg(igrid_mod)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use camxcom
      use grid
      use camxfld
      use tracer
      use procan
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c         Calls all of the routines that perform the averaging for 
c         the various CSA concentation fields. This will do all of
c         CSA dubdomains.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c       igrid_mode  grid number of modeling grid
c     Output:  
c
c    Called by:
c       CAMX
c    Subroutines called:
c       AVERAGE
c
c     Copyright 1996 - 2025
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      implicit none
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument delarations:
c-----------------------------------------------------------------------
c
      integer :: igrid_mod
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: icsa_grid
      integer :: ncol_csa
      integer :: nrow_csa
      integer :: nlay_csa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      do icsa_grid = 1,npadom
          if( ipagrd(icsa_grid) .EQ. igrid_mod ) then
            ncol_csa = i_ne(icsa_grid) - i_sw(icsa_grid) + 1
            nrow_csa = j_ne(icsa_grid) - j_sw(icsa_grid) + 1
            nlay_csa = t_lay(icsa_grid) - b_lay(icsa_grid) + 1
            call average_csa(icsa_grid,ncol(igrid_mod),nrow(igrid_mod),nlay(igrid_mod),
     &        igrid_mod,deltat(igrid_mod)/2.0,ncol_csa,nrow_csa,nlay_csa,
     &            ntotsp,tempk(iptr3d(igrid_mod)),press(iptr3d(igrid_mod)),
     &               csaconc(iptrcsa(icsa_grid)),csaavrg(iptrcsa(icsa_grid)),
     &                                       ipacl_3d(iptr3d_full(igrid_mod)) )
          endif
      enddo
      return
      end
