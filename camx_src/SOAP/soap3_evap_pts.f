C**** SOAP3_EVAP_PTS
c
      subroutine soap3_evap_pts( idx_pntid,idx_layer_beg,idx_layer_end,
     &                      idx_svoa,idx_hoa,evap_factors,temperature,
     &                          pntemis_soap3,pntemis_svoa,pntemis_hoa)
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
c        idx_pntid        idxex of this point source in global arrays
c        idx_layer_beg    index of the layer to inject the emissions into
c        idx_layer_end    index of the layer to inject the emissions into
c        temperature      temperature in aloft layers
c        pntemis_soap3    point emissioons rate for POA_XX (grams/s)
c     Output arguments:
c        idx_hoa          index of HOA in species array
c        idx_svoa          index of HOA in species array
c        pntemis_svoa     point emissions rate for SVOA (moles/s)
c        pntemis_hoa      point emissions rate for HOA (grams/s)
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
      integer idx_pntid
      integer idx_layer_beg
      integer idx_layer_end
      integer idx_hoa
      integer idx_svoa
      real    evap_factors(3,NUM_SOAP_EMISS)
      real    temperature(*)
      real    pntemis_soap3(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB,MXLAYER)
      real    pntemis_svoa(MXLAYER)
      real    pntemis_hoa(MXLAYER)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer imod, iems, klay, isoa, ipt
      integer idx_soap_spec(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      real    evap_fact, evap_mass, convfac, poa_emiss
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c ---- loop over model species ---
c
      idx_svoa = -9
      idx_hoa = -9
      idx_soap_spec = -9
      do iems=1,nemspc
c
c ---- find the emissions species that are used in SOAP3 evap ---
c
          if( emspcname(iems)(1:5) .EQ. 'SVOA ' ) idx_svoa = iems
          if( emspcname(iems)(1:4) .EQ. 'POA ' ) idx_hoa = iems
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
c  ---- loop over layers ---
c
      pntemis_hoa = 0.
      pntemis_svoa = 0.
      do klay=idx_layer_beg,idx_layer_end
c
c  --- see if POA species was found ---
c
           do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
              if( idx_soap_spec(isoa) .LE. 0 ) cycle
              if( pntemis_soap3(isoa,klay) .LE. 0. ) cycle
c
c ---- evaporate POA into SVOA ---
c
              evap_fact = MIN(1.0,
     &                    (evap_factors(1,isoa) * temperature(klay)**2 + 
     &                     evap_factors(2,isoa) * temperature(klay) + 
     &                     evap_factors(3,isoa) ))
              evap_mass = pntemis_soap3(isoa,klay) * (1. - evap_fact)
              if( ksvoa .GT. 0 .and. ksvoa .NE. nspec+1 )
     &              pntemis_svoa(klay) = pntemis_svoa(klay) + 
     &                 MAX(0., pntemis_svoa(klay) + evap_mass/mole_weight(ksvoa))
c
c ---- back out the vcap emissions from each POA species ---
c
              pntemis_hoa(klay) =  pntemis_hoa(klay) +
     &             MAX( 0., pntemis_soap3(isoa,klay) - evap_mass)
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
