c*** LOADCSA_SENS
c
      subroutine loadcsa_sens(filflg,numcol,numrow,numlay,ncsa,gridval,
     &                  icl_csa,jcl_csa,kcl_csa,numspc_csa,nspc,
     &                  ngas,nrad,nhet,sens,convfac)
      use tracer
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine loads the sensitivies which are stored in a
c     4-D array from/into a 2-D array in conclusion/preperation 
c     of the CSA chemistry routine.  The 2-D array contains the 
c     family of sensitivities for one cell and all modeled species.
c     The flag "filflg" determines wether the values in the 4-D
c     array are loaded into the 2-D array or vice-versa.  Sensitivities
c     are converted from umol/m3 to ppm for chemistry and back.
c
c     Copyright 1996 - 2025
c     Ramboll
c
c   Argument descriptions:
c       filflg     L  flag for determining which direction to fill
c                     .TRUE.  = put 2-D values into 4-D gridded array
c                     .FALSE. = put 4-D values into 2-D array
c       numcol     I  number of cells in X direction
c       numrow     I  number of cells in Y direction
c       numlay     I  number of layers 
c       ncsa       I  number of total CSA species
c       gridval    R  4-D array of CSA values
c       icl_csa    I  the X grid location of current cell
c       jcl_csa    I  the Y grid location of current cell
c       kcl_csa    I  the vertical grid location of current layer
c       numspc_csa I  number of CSA families
c       nspc       I  number of modeled species
c       ngas       I  number of gaseous species
c       nrad       I  number of radical species
c       nhet       I  number of heterogeneous rxns
c       sens       R  2-D array for this cell and species
c       convfac    R  conversion from ppm to umol/m3
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
      logical filflg
      integer numcol
      integer numrow
      integer numlay
      integer ncsa
      real    gridval(numcol,numrow,numlay,ncsa)
      integer icl_csa
      integer jcl_csa
      integer kcl_csa
      integer numspc_csa
      integer nspc
      integer nrad
      real    sens(numspc_csa,nspc+nhet)
      real    convfac
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
c  --- task 1: load 4-D array
c
      if( filflg ) then
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
            gridval(icl_csa,jcl_csa,kcl_csa,idxcsa) = sens(icsa,ispc)*convfac
c
c  --- next family ---
c
          enddo
c
c  --- next species ---
c
        enddo
c
c  --- radical species ---
c
        do ispc=1,nrad
          do icsa=1,numspc_csa
            idxcsa = iptddm(ispc)+icsa-1
            gridval(icl_csa,jcl_csa,kcl_csa,idxcsa) = sens(icsa,ispc)
          enddo
        enddo
c
c  --- particulate species & k_hetero -> no need to load these sensitivities back
c
c
c  --- task 2: load 2-D array
c
      else
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
            sens(icsa,ispc) = gridval(icl_csa,jcl_csa,kcl_csa,idxcsa)/convfac
c
c  --- next family ---
c
          enddo
c
c  --- next species ---
c
        enddo
c
c  --- radical species ---
c
        do ispc=1,nrad
          do icsa=1,numspc_csa
            idxcsa = iptddm(ispc)+icsa-1
            sens(icsa,ispc) = gridval(icl_csa,jcl_csa,kcl_csa,idxcsa)
          enddo
        enddo
c
c  --- particulate species ---
c
        do ispc=ngas+1,nspc
          do icsa=1,numspc_csa
            idxcsa = iptddm(ispc)+icsa-1
            sens(icsa,ispc) = gridval(icl_csa,jcl_csa,kcl_csa,idxcsa)
          enddo
        enddo
c
c  --- k_hetero ---
c
        do ispc=nspc+1,nspc+nhet
          do icsa=1,numspc_csa
            sens(icsa,ispc) = 0.0 ! reset to zero
          enddo
        enddo
c
c  --- done ---
c
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
