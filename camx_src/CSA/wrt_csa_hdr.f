      subroutine wrt_csa_hdr(endtim,enddate)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use camxcom
      use grid
      use procan
      use tracer
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c         Writes the metadata and the other header variables to the
c         NetCDF for the CSA grids.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        endtim  - model end time (HHMM)
c        enddate - model end date (YYJJJ)
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
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument delarations:
c-----------------------------------------------------------------------
c
      integer :: enddate
      real    :: endtim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 :: action
c
      character*20, dimension(MXTRSP) :: spec_units
      character*20, dimension(MXTRSP) :: spec_long_name
      character*60, dimension(MXTRSP) :: spec_desc
      character*60, dimension(MXTRSP) :: spec_coords
c
      integer :: icsa_grid           
      integer :: nocl_grid           
      integer :: nrow_grid           
      integer :: nlay_grid           
      integer :: igrid_mod
      integer :: ncol_csa
      integer :: nrow_csa
      integer :: nlay_csa
c
      real :: dx_csa
      real :: dy_csa
      real :: dx_grid
      real :: dy_grid
      real :: orgx_csa
      real :: orgy_csa
      real :: orgx_grid
      real :: orgy_grid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- loop over grids ---
c
      do icsa_grid = 1,npadom
         igrid_mod = ipadom(icsa_grid)
         dx_csa   = 1000.*(delx/float(meshold(igrid_mod)))
         dy_csa   = 1000.*(dely/float(meshold(igrid_mod)))
         orgx_grid = xorg + delx*(inst1(igrid_mod)-1) - dx_grid
         orgy_grid = yorg + dely*(jnst1(igrid_mod)-1) - dy_grid
         orgx_csa = 1000.*(orgx_grid + (i_sw(icsa_grid) - 1)*dx_grid)
         orgy_csa = 1000.*(orgy_grid + (j_sw(icsa_grid) - 1)*dy_grid)
         ncol_csa = i_ne(icsa_grid) - i_sw(icsa_grid) + 1
         nrow_csa = j_ne(icsa_grid) - j_sw(icsa_grid) + 1
         nlay_csa = t_lay(icsa_grid) - b_lay(icsa_grid) + 1
         if( llatlon ) then
            dx_csa   = dx_csa/1000.
            dy_csa   = dy_csa/1000.
            orgx_csa = orgx_csa/1000.
            orgy_csa = orgy_csa/1000.
         endif
c
         call ncf_set_vars_base()
         call ncf_set_tstep(begdate,begtim,enddate,endtim)
         write(action,'(A,I2)') 'Writing output file for CSA domain: ',
     &                                                         icsa_grid
         call ncf_set_specatt_ddm(spec_units,spec_long_name,spec_desc,
     &                                                        spec_coords)
         call ncf_set_global_csa('CSA       ',icsa_grid,begdate,begtim,
     &                                                  enddate,endtim)
         call ncf_wrt_dim(action,iow_csa(icsa_grid),icsa_grid,ncol_csa,
     &                                          nrow_csa,nlay_csa,ntotsp)
         call ncf_wrt_global(action,iow_csa(icsa_grid),ntotsp,ptname,.FALSE.)
         call ncf_wrt_vars_base(action,iow_csa(icsa_grid))
         call ncf_wrt_vars_species(action,iow_csa(icsa_grid),ncol_csa,
     &                     nrow_csa,nlay_csa,ntotsp,ptname,spec_units,
     &                           spec_long_name,spec_desc,spec_coords,4)
         call ncf_enddef_file(action,iow_csa(icsa_grid))
         call ncf_wrt_data_gridcsa(action,iow_csa(icsa_grid),icsa_grid,
     &             ncol_csa,nrow_csa,orgx_csa,orgy_csa,dx_csa,dy_csa,nlay_csa)
      enddo
c
      return
      end
