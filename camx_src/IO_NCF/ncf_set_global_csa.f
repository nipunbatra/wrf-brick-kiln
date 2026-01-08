c**** NCF_SET_GLOBAL_CSA
c
      subroutine ncf_set_global_csa(inname,this_csagrd,begin_date,
     &                        begin_time,ending_date,ending_time)
      use camxcom
      use grid
      use procan
      use tracer
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the global file attributes for the NetCDF file
c
c     Copyright 1996 - 2025
c     Ramboll
c      Argument description:
c       Inputs:
c         inname      C filetype of this kind of file
c         this_csagrd I grid index
c         begin_date  I model begin date (YYJJJ)
c         begin_time  R model begin time
c         ending_date I model end date (YYJJJ)
c         ending_time R model end time
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'ncf_iodat.inc'
      include 'filunit.inc'
      include 'flags.inc'
      include 'chmdat.inc'
      include 'namelist.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*10 inname
      integer      this_csagrd
      integer      begin_date
      real         begin_time
      integer      ending_date
      real         ending_time
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: igrid_mod
      real    :: dx_grid
      real    :: dy_grid
      real    :: orgx_grid
      real    :: orgy_grid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- get current date/time ---
c
      call getime(ncf_cdate,ncf_ctime)
      ncf_wdate = ncf_cdate
      ncf_wtime = ncf_ctime
c
c --- domain definition attributes ---
c
      igrid_mod = ipagrd(this_csagrd)
      ncf_iutm = iuzon
      ncf_istag = 0
      ncf_cproj = 2
      dx_grid   = delx/float(meshold(igrid_mod))
      dy_grid   = dely/float(meshold(igrid_mod))
      orgx_grid = xorg + delx*(inst1(igrid_mod)-1) - dx_grid
      orgy_grid = yorg + dely*(jnst1(igrid_mod)-1) - dy_grid
      ncf_xorig = 1000.*(orgx_grid + (i_sw(this_csagrd) - 1)*dx_grid)
      ncf_yorig = 1000.*(orgy_grid + (j_sw(this_csagrd) - 1)*dy_grid)
      ncf_xcell = 1000.*dx_grid
      ncf_ycell = 1000.*dy_grid
      if( llatlon ) then
         ncf_cproj = 0
         ncf_gdtyp = 1
         ncf_xorig = ncf_xorig/1000.
         ncf_yorig = ncf_yorig/1000.
         ncf_xcell = ncf_xcell/1000.
         ncf_ycell = ncf_ycell/1000.
      else if( lutm ) then
         ncf_cproj = 1
         ncf_gdtyp = 5
      else if( lambrt ) then
         ncf_cproj = 2
         ncf_gdtyp = 2
      else if( lrpolar ) then
         ncf_cproj = 3
         ncf_gdtyp = 4
      else if( lpolar ) then
         ncf_cproj = 4
         ncf_gdtyp = 6
      else if( lmerc ) then
         ncf_cproj = 5
         ncf_gdtyp = 7
      endif
      ncf_xcent = polelon
      ncf_ycent = polelat
      ncf_p_alp = tlat1
      ncf_p_bet = tlat2
      ncf_p_gam = polelon
      ncf_ncols = i_ne(this_csagrd) - i_sw(this_csagrd) + 1
      ncf_nrows = j_ne(this_csagrd) - j_sw(this_csagrd) + 1
      ncf_nlays = t_lay(this_csagrd) - b_lay(this_csagrd) + 1

      ncf_nthik = 1
      ncf_nvars = nspec + NCF_BASE_VARS
c
c --- file description attributes ---
c
      ncf_name = inname
      ncf_note = runmsg
      ncf_itzon = itzon
      ncf_ftype = 1
      ncf_vgtyp = 6
      ncf_vgtop = 10000.
      ncf_vglvls = 0.
      ncf_gdnam = "CAMx v7.32"
      ncf_upnam = "CAMx v7.32"
      ncf_filedesc = inname
      ncf_sdate = Start_Date_Hour(1)*1000 + MOD(begin_date,1000)
      ncf_stime = begin_time*100
      if(dtout .LT. 60. ) then
        ncf_tstep = INT(dtout)*100
      else
        ncf_tstep = INT(dtout/60.)*10000
      endif
      ncf_grid_id = ipagrd(this_csagrd)
      ncf_i_grid_start = i_sw(this_csagrd)
      ncf_i_grid_end = i_ne(this_csagrd)
      ncf_j_grid_start = j_sw(this_csagrd)
      ncf_j_grid_end = j_ne(this_csagrd)
      ncf_grid_mesh_factor = meshold(ipagrd(this_csagrd))
c
c --- simulation description attributes ---
c
      ncf_flexi_nest = 0
      if( lflexi ) ncf_flexi_nest = 1
      ncf_advection = Advection_Solver
      ncf_chem_solver = Chemistry_Solver
      ncf_pig = PiG_Submodel
      ncf_probing_tool = Probing_Tool
      ncf_total_species = nspec
      ncf_radical_species = nrad
      ncf_gas_species = ngas
      ncf_pm_species = naero
      ncf_reactions = nreact
      ncf_drydep = Drydep_Model
      ncf_wetdep = 0
      if( lwet ) ncf_wetdep = 1
      ncf_acm2 = 0
      if( lacm2 ) ncf_acm2 = 1
      ncf_cig_model = 0
      if( Subgrid_Convection ) ncf_cig_model = 1
      ncf_surface_model = 0
      if( lsrfmod ) ncf_surface_model = 1
      ncf_inline_ix_emiss = 0
      if( lixemis ) ncf_inline_ix_emiss = 1
      ncf_bidi_nh3_drydep = 0
      if( lbidinh3 ) ncf_bidi_nh3_drydep = 1
      ncf_super_stepping = 0
      if( lsuper ) ncf_super_stepping = 1
      ncf_gridded_emiss = 0
      if( larsrc ) ncf_gridded_emiss = 1
      ncf_point_emiss = 0
      if( lptsrc ) ncf_point_emiss = 1
      ncf_ignore_emiss_dates = 0
      if( le1day ) ncf_ignore_emiss_dates = 1
      ncf_output_3d = 0
      if( l3davg(1) ) ncf_output_3d = 1
      ncf_pig_sample_grid = 0
      ncf_conventions = "CF-1.6"
      ncf_history = 'Generated by '//PRMVERSION
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
