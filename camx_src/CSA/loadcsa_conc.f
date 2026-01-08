c*** LOADCSA_CONC
c
      subroutine loadcsa_conc(ncol_csa,nrow_csa,nlay_csa,ncsa,gridcsa,
     &         icl_csa,jcl_csa,kcl_csa,icl_mod,jcl_mod,kcl_mod,
     &         numspc_csa,nspc,ngas,nrad,nhet,ncol_mod,nrow_mod,nlay_mod,gridmod)
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine loads the regular model concentrations into the
c     appropriate location for the CSA grids. This is done each hour.
c
c     Copyright 1996 - 2025
c     Ramboll
c
c   Argument descriptions:
c       ncol_csa   I  number of cells in X direction in CSA grid
c       nrow_csa   I  number of cells in Y direction in CSA grid
c       nlay_csa   I  number of layers in CSA grid
c       ncsa       I  number of total CSA species
c       gridcsa    R  4-D array of CSA values
c       icl_csa    I  the X grid location of current cell of CSA grid
c       jcl_csa    I  the Y grid location of current cell of CSA grid
c       kcl_csa    I  the vertical grid location of current layer of CSA grid
c       icl_mod    I  the X grid location of current cell of model grid
c       jcl_mod    I  the Y grid location of current cell of model grid
c       kcl_mod    I  the vertical grid location of current layer of model grid
c       numspc_csa I  number of CSA families
c       nspc       I  number of modeled species
c       ngas       I  number of gaseous species
c       nrad       I  number of radical species
c       nhet       I  number of heterogeneous rxns
c       ncol_mod   I  number of cells in X direction in modelling grid
c       nrow_mod   I  number of cells in Y direction in modelling grid
c       nlay_mod   I  number of layers in modelling grid
c       gridmod    R  gridded array of regular model  concentrations
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer ncol_csa
      integer nrow_csa
      integer nlay_csa
      integer ncsa
      real    gridcsa(ncol_csa,nrow_csa,nlay_csa,ncsa)
      integer icl_csa
      integer jcl_csa
      integer kcl_csa
      integer icl_mod
      integer jcl_mod
      integer kcl_mod
      integer numspc_csa
      integer nspc
      integer ngas
      integer nrad
      integer nhet
      integer ncol_mod
      integer nrow_mod
      integer nlay_mod
      real    gridmod(ncol_mod,nrow_mod,nlay_mod,nspc)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   icsa, ispc, idxcsa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- loop over the modeled species ---
c
      do ispc=nrad+1,ngas
c
c  --- loop over the number of CSA families ---
c
         do icsa=1,numspc_csa
c
c  --- calculate the index into the CSA species list ---
c
            idxcsa = iptddm(ispc)+icsa-1
            gridcsa(icl_csa,jcl_csa,kcl_csa,idxcsa) = 
     &                  gridmod(icl_mod,jcl_mod,kcl_mod,ispc)
c
c  --- next family ---
c
         enddo
c
c  --- next species ---
c
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
