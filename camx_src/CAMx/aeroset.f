      subroutine aeroset(dsec_i,ierr)
      use camxcom
      use grid
      use chmstry
      use filunit
c
c----CAMx v7.32 250801
c
c     AEROSET prepares for the AERO chemistry routines in CAMx
c
c      Copyright 1996 - 2025
c     Ramboll
c  
c     Modifications:
c        12/29/06   --bkoo--     Added species mapping for the updated SOA scheme
c        01/10/07   --bkoo--     Included camx_aero.inc
c                                Set pointers for aqueous PM
c        03/18/14   --bkoo--     Added pointers for benzene SOA
c        12/07/14   --bkoo--     Set pointers for VBS species
c        01/08/16   --bkoo--     Updated for SAPRC99->SAPRC07TC
c                                            Revised VBS
c        08/25/16   --bkoo--     Updated for new SOAP
c        09/16/16   --pk--       Added pointers for GLY, MGLY & GLYD
c        09/16/16   --bkoo--     Set WTMOL here (moved from CHMDAT.F)
c        10/16/17   --bkoo--     Added sfacJSOA; removed LSOAP2
c        11/22/17   --bkoo--     Updated WTMOL and separated Ca & Mg from CO3
c        02/21/20   --rlb--      Added CB6r5 (idmech=1)
c        03/24/21   --gy--       Added CB7 (idmech=7)
c        02/15/22   --gy--       Added RACM2s21 (idmech=9)
c        11/03/23   --cemery--   Removed CB05 (Mech 6)
c        05/17/24   --cemery--   Removed CMU PM option
c
c     Input arguments: 
c        dsec_i              User supplied section cutpoints 
c 
c     Output arguments: 
c        ierr                Error flag
c            
c     Routines called: 
c        UPTIME
c            
c     Called by: 
c        READCHM
c 
      include 'camx.prm'
      include 'soap.inc'
      include 'vbs.inc'
      include 'section.inc'
      include 'section_aq.inc'
      include 'camx_aero.inc'
c
      real*8  dsec_i(nsecp1)
      integer ierr
c
      real    dtio
      integer nsteps, i, n
c     
c-----Entry point 
c
      ierr = 0
c
c-----Reset dtaero to hit output interval exactly
c
      dtio = amin1(60.,dtinp,dtems,dtout)
      nsteps = INT( 0.999*dtio/dtaero ) + 1
      dt_aero = dtio/nsteps
      do n = 1,ngrid
        time_aero(n) = time
        date_aero(n) = date
        call uptime(time_aero(n),date_aero(n),60.*dt_aero)
        aero_dt(n)  = 0.0
      enddo
c
c-----Set pointers for aqueous PM
c
      if (idmech.EQ.5) then ! SAPRC07TC
        khpo_c = kh2o2
        kfoa_c = kfacd
        kmhp_c = kcooh
        kpaa_c = kco3h
        kohp_c = krooh
        kopa_c = kro3h
        kgly_c = kgly
        kmgly_c = kmgly
        kglyd_c = kglyd
      elseif (idmech.EQ.9) then ! RACM2s21
        khpo_c = kh2o2
        kfoa_c = kora1
        kmhp_c = kop1
        kpaa_c = kpaa
        kohp_c = kop2
        kopa_c = nspec+1
        kgly_c = kgly
        kmgly_c = kmgly
        kglyd_c = nspec+1
      elseif (idmech.EQ.1 .or. idmech.EQ.3 
     &         .or. idmech.EQ.4 .or. idmech.EQ.7 ) then ! CB6 or CB7
!bkoo - to be updated: Add ISPX to hydroperoxides for CB6?
        khpo_c = kh2o2
        kfoa_c = kfacd
        kmhp_c = kmepx
        kpaa_c = kpacd
        kohp_c = krooh
        kopa_c = nspec+1
        kgly_c = kgly
        kmgly_c = kmgly
        kglyd_c = kglyd
      endif
c
      if (lvbs .and. LVBSPREPROC) then
        kpap_c(0) = kpap0
        kpap_c(1) = kpap1
        kpap_c(2) = kpap2
        kpap_c(3) = kpap3
        kpap_c(4) = kpap4
        kpcp_c(0) = kpcp0
        kpcp_c(1) = kpcp1
        kpcp_c(2) = kpcp2
        kpcp_c(3) = kpcp3
        kpcp_c(4) = kpcp4
        kpfp_c(0) = kpfp0
        kpfp_c(1) = kpfp1
        kpfp_c(2) = kpfp2
        kpfp_c(3) = kpfp3
        kpfp_c(4) = kpfp4
      endif
c
c-----SOAP parameters based on new chamber data (Hodzic et al., 2016, ACP)
c
      if( luse_soap3_simple ) then
         y_h      = yh_3simple
         y_l      = yl_3simple
         y_n      = yn_3simple
         csatS    = csatS_3simple
         cstempS  = cstempS_3simple
         deltahS  = deltahS_3simple
         kpolya   = kpolya_3simple
         kpolyb   = kpolyb_3simple
         sfacJSOA = sfacJSOA_3simple
      else if( luse_soap3_complx ) then
         y_h      = yh_3complx
         y_l      = yl_3complx
         y_n      = yn_3complx
         csatS    = csatS_3complx
         cstempS  = cstempS_3complx
         deltahS  = deltahS_3complx
         kpolya   = kpolya_3complx
         kpolyb   = kpolyb_3complx
         sfacJSOA = sfacJSOA_3complx
      endif
c
c-----Molecular weights for PM chemistry modules;
c     Used in AEROCHEM_CF, AEROCHEM_CMU, and DDM routines
c
      wtmol( 1) =  96.1 ! SO4
      wtmol( 2) =  18.  ! NH4
      wtmol( 3) =  62.  ! NO3
      wtmol( 4) =  35.5 ! Cl
      wtmol( 5) =  23.  ! Na
      wtmol( 6) =  39.1 ! K
      wtmol( 7) =  60.  ! CO3
      wtmol( 8) =  40.1 ! Ca
      wtmol( 9) =  24.3 ! Mg
      wtmol(10) =  55.8 ! Fe
      wtmol(11) =  54.9 ! Mn
      if (lvbs) then
        ksoac_c = kPBS0
        wtmol(12) = mwPBS0
      else
        ksoac_c = kboa0
        wtmol(12) = mwboa0
      endif
c
      return
      end
