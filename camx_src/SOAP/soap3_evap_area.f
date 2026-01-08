C**** SOAP3_EVAP_AREA
c
      subroutine soap3_evap_area( igrd, ncol_grid, nrow_grid,
     &                    nlay_grid, nlay_ems, evap_factors,
     &                                  temperature, area_emis, numspcs)
      use chmstry
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     This routine applies factors to evaporate POA to SVOA and
c     updates the emissions array so emissions are injected into
c     the appropriate model species.
c
c      Copyright 1996 - 2025
c     Ramboll
c
c     Input arguments:
c        igrd             grid index
c        ncol_grid        number of columns
c        nrow_grid        number of rows
c        nlay_grid        number of layers in grid
c        nlay_ems         number of layers in emissions foles
c        numspcs          number of model species
c        temperature      temperature in aloft layers
c     Output arguments:
c        area_emis        area emissions rate (grams/s)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/24   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'soap3.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer ncol_grid
      integer nrow_grid
      integer nlay_grid
      integer nlay_ems
      integer numspcs
      real    evap_factors(3,NUM_SOAP_EMISS)
      real    temperature(ncol_grid,nrow_grid,nlay_grid)
      real    area_emis(ncol_grid,nrow_grid,nlay_ems,numspcs)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer imod, iems, idx_soap_spec(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      integer icl, jcl, klay, isoa, ipt, idx_svoa, idx_hoa
      real    evap_fact, evap_mass, convfac
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c ---- loop over model species ---
c
       idx_svoa = -9
       idx_soap_spec = -9
       do iems=1,nemspc
c
c ---- find the emissions species that are used in SOAP3 evap ---
c
          if( emspcname(iems)(1:5) .EQ. 'SVOA ' ) idx_svoa = iems
          if( emspcname(iems)(1:4) .EQ. 'POA ' .OR. 
     &                    emspcname(iems)(1:4) .EQ. 'HOA ') idx_hoa = iems
          do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
            if( emspcname(iems) .EQ. soap_emiss_names(isoa) ) 
     &                                         idx_soap_spec(isoa) = iems
          enddo
       enddo
c
c  ---- don't change anything if SVOA or HOA not found ---
c
       if( idx_svoa .LE. 0 ) goto 9999
       if( idx_hoa .LE. 0 ) goto 9999
c
c  --- loop over grid and update all affected species ----
c
       do klay=1,nlay_ems
         do jcl=1,nrow_grid
           do icl=1,ncol_grid
             area_emis(icl,jcl,klay,idx_hoa) = 0. 
c
c  --- see if POA species was found ---
c
             do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
                if( idx_soap_spec(isoa) .LE. 0 ) cycle
                if( area_emis(icl,jcl,klay,idx_soap_spec(isoa)) .LE. 0 ) cycle
c
c ---- evaporate POA into SVOA ---
c
                 evap_fact = MIN(1.0,
     &                       (evap_factors(1,isoa) * temperature(icl,jcl,klay)**2 + 
     &                        evap_factors(2,isoa) * temperature(icl,jcl,klay) + 
     &                        evap_factors(3,isoa) ))
                 evap_mass = area_emis(icl,jcl,klay,idx_soap_spec(isoa)) * 
     &                       (1. - evap_fact)
                 if( ksvoa .NE. 0 .and. ksvoa .NE. numspcs+1 )
     &                area_emis(icl,jcl,klay,idx_svoa) = 
     &                     MAX(0.,area_emis(icl,jcl,klay,idx_svoa) +
     &                                 evap_mass/mole_weight(ksvoa) )
c
c ---- back out the evap emissions from each POA species ---
c
                 area_emis(icl,jcl,klay,idx_hoa) = area_emis(icl,jcl,klay,idx_hoa) + 
     &              MAX( 0.,area_emis(icl,jcl,klay,idx_soap_spec(isoa)) - evap_mass)
             enddo
           enddo
         enddo
       enddo
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
