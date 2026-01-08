C**** SA_SOAP3_EVAP_PTS
c
      subroutine sa_soap3_evap_pts( idx_pntid,idx_layer_beg,idx_layer_end,
     &                     numpts,num_trac_emiss,istart_hoa,iend_hoa,
     &                   istart_svoa,iend_svoa,evap_factors,temperature,
     &                               sa_pnt_emiss,pntemis_hoa,pntemis_svoa)
      use chmstry
      use tracer
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
c        numpts           number of points in arrays
c        num_trac_emiss   number of tracer emissions species
c        istart_hoa       starting index in array for HOA emissions
c        iend_hoa         ending index in array for HOA emissions
c        istart_svoa      starting index in array for SVOA emissions
c        iend_svoa        ending index in array for SVOA emissions
c        temperature      temperature in aloft layers
c        sa_pnt_emiss     point emissions rate for SA
c     Output arguments:
c        pntemis_hoa      point emissions by layer HOA tracer species
c        pntemis_svoa     point emissions by layer HOA tracer species
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
      integer numpts
      integer istart_hoa
      integer iend_hoa
      integer istart_svoa
      integer iend_svoa
      integer num_trac_emiss
      real    evap_factors(3,NUM_SOAP_EMISS)
      real    temperature(*)
      real    sa_pnt_emiss(numpts,num_trac_emiss)
      real    pntemis_svoa(istart_svoa:iend_svoa,MXLAYER)
      real    pntemis_hoa(istart_hoa:iend_hoa,MXLAYER)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer klay, isoa, ispc_hoa, ispc_svoa, ispc_poa
      integer idx_soap_spec(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      real    evap_fact, evap_mass, convfac
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
       pntemis_hoa = 0.
       pntemis_svoa = 0.
c
c
c ---- loop over model species ---
c
      ispc_svoa = istart_svoa - 1
      ispc_poa = -1
      do ispc_hoa=istart_hoa,iend_hoa
         ispc_svoa = ispc_svoa + 1
         ispc_poa = ispc_poa + 1
c
c  ---- loop over layers ---
c
          do klay=idx_layer_beg,idx_layer_end
c
c  --- see if POA species was found ---
c
              do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
                 if( sa_pnt_emiss(idx_pntid,soap3_iemcls(isoa)+ispc_poa) .LE. 0. ) cycle
c
c ---- evaporate POA into SVOA ---
c
                 evap_fact = MIN(1.0,
     &                       (evap_factors(1,isoa) * temperature(klay)**2 + 
     &                        evap_factors(2,isoa) * temperature(klay) + 
     &                        evap_factors(3,isoa) ))
                 evap_mass = sa_pnt_emiss(idx_pntid,soap3_iemcls(isoa)+ispc_poa) * 
     &                       (1. - evap_fact)
                 pntemis_svoa(ispc_svoa,klay) = 
     &                pntemis_svoa(ispc_svoa,klay) + 
     &                 MAX(0., pntemis_svoa(ispc_svoa,klay) + 
     &                                  evap_mass/mole_weight(ksvoa))
c
c ---- back out the vcap emissions from each POA species ---
c
                 pntemis_hoa(ispc_hoa,klay) = pntemis_hoa(ispc_hoa,klay) +
     &             MAX( 0., sa_pnt_emiss(idx_pntid,soap3_iemcls(isoa)+ispc_poa) - evap_mass)
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
