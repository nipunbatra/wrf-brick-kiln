c*** PSATSA
c
      subroutine psatsa(numcol,numrow,numlay,igrid,icell,jcell,kcell,
     &                                         nspc,delcon,nprec,yfac)
      use grid
      use chmstry
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the "chemistry" adjustments to the tracer 
c     species for PSAT.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     the regular model chemistry.  This is essentially an adjustment
c     for production or decay.
c
c      Copyright 1996 - 2025
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c       numcol  I number of columns in this slice
c       numrow  I number of rows in this slice
c       numlay  I number of layers in this slice
c       igrid     I  grid number
c       icell     I  the X grid location of current cell
c       jcell     I  the X grid location of current cell
c       kcell     I  the vertical grid location of current layer
c       nspc      I  number of model species
c       delcon    R  array of change in concentrations total
c       nprec     I  number of SOA precursors
c       yfac      R  array of weighting factors for SOA
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     09/28/03   --gwilson--  Original development
c     12/29/06   --bkoo--     Revised for the updated SOA scheme
c     05/04/07   --gwilson--  Changed call to cyctpnsa so that conversion
c                             factor is set for lower bound values -
c                             the conversion is passed here from chemdriv
c     03/18/14   --bkoo--     Revised for benzene SOA
c     08/25/16   --bkoo--     Updated for new SOAP
c     10/16/17   --bkoo--     Fixed DELARO that was missing IVOA changes
c                             Updated for SOA photolysis
c     01/12/18   --bkoo--     Removed BNZA/TOLA/XYLA/ISP/TRP
c     01/09/19   --cemery--   Added DMS
c     06/18/19   --bkoo--     Updated delARO/delTRP/delSQT for SAPRC07T
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
      integer   numcol
      integer   numrow
      integer   numlay
      integer   igrid
      integer   icell
      integer   jcell
      integer   kcell
      integer   nspc
      real      delcon(7,nspc+1)
      integer   nprec
      real      yfac(3,nprec)
c
c --- The nprec dimension of yfac is
c     1 - BNZA; 2 - TOLA; 3 - XYLA; 4 - IVOA
c     5 - SVOA; 6 - ISP ; 7 - TRP ; 8 - APN; 9 - SQT
c     except yfac(1,1) and (2,1) are a special case and contain
c     weighting factors for computing NOx-dependent yields from ARO
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   iacg2, iacg1, iaoa0, iaro, iisp, itrp, iapin, isqt
      integer   iso2, idms, ihg0, i
      integer*8 idxcel, idx, jdx, kdx, ldx, mdx, idx0, idx2
      real      delARO, delISP, delTRP, delSQT, delso2, delps4, deldms, delhg2
      real      delACG2, delACG1, delNV1, delIVOA, delSVOA, delAPIN
      real      sumaro, sumisp, sumtrp, sumsqt, sumso2, sumps4, sumdms
      real      sumhg0, sumhg2, sumivoa, sumsvoa, sumapin

      real, parameter :: wtboa0 = 220.0 ! must be consistent with MWBOA0 in SOAP.INC
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  DBLE(ipsa3d(igrid)-1) + DBLE(icell) + 
     &            DBLE(numcol) * DBLE(jcell-1) + 
     &                 DBLE(numcol) * DBLE(numrow) * DBLE(kcell-1)
c
c----------------------------------------------------------------------
c  ---  SOA tracers ----
c----------------------------------------------------------------------
c
      if( lsoa ) then
c
c---- assign these pointers to those used for RACM
c
         delARO = delcon(1,kbenz) + delcon(1,ktol ) + delcon(1,kxyl )   !CB6 & CB7 & SAPRC(BENZ) & RACM2(TOL)
     &          + delcon(1,ktolu) + delcon(1,karo1)                     !SAPRC
     &          + delcon(1,kmxyl) + delcon(1,koxyl) + delcon(1,kpxyl)   !SAPRC
     &          + delcon(1,kb124) + delcon(1,karo2)                     !SAPRC
     &          + delcon(1,kben )                                       !RACM2
     &          + delcon(1,kxym ) + delcon(1,kxyo ) + delcon(1,kxyp )   !RACM2

         delISP = delcon(1,kisop)                                       !CB6 & CB7 & SAPRC
     &          + delcon(1,kiso )                                       !RACM2

         delTRP = delcon(1,kterp)                                       !CB6 & CB7
     &          + delcon(1,klim)                                        !RACM2

         delAPIN = delcon(1,kapin)                                      !CB7
     &           + delcon(1,kapi)                                       !RACM2

         delSQT = delcon(1,ksqt )                                       !CB6 & CB7 & RACM2
     &          + delcon(1,ksesq)                                       !SAPRC

         delACG2 = delcon(1,kacg2)

         delACG1 = delcon(1,kacg1)

         delNV1 = delcon(3,kaoa0) - delcon(5,kaoa0) ! adjust non-vol Anthro-SOA delta [ug/m3] to exclude
     &                            - delcon(6,kaoa0) ! changes due to Anthro-SOA polymerization/photolysis

         delIVOA = delcon(1,kivoa )                                     !CB6 & CB7 & RACM2 & SAPRC

         delSVOA = delcon(1,ksvoa )                                     !CB6 & CB7 & RACM2 & SAPRC
c
c  --- get the sum of the tracer species for ARO/AG2 yield ---
c
         iacg2 = iptcls(idxipt(ITRAG2)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            iacg2 = iacg2 + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(iacg2)*yfac(1,1)
     &                      + ylrates(iacg2)*yfac(2,1) )
         enddo
c
c  --- make adjustment to AG2 tracers based on change in AG2 and
c      yield rates for ARO ---
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRAG2)),nptcls(idxipt(ITRAG2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delACG2 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO/ACG1 yield ---
c
         iacg1 = iptcls(idxipt(ITRAG1)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            iacg1 = iacg1 + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(iacg1)*yfac(1,1)
     &                      + ylrates(iacg1)*yfac(2,1) )
         enddo
c
c  --- make adjustment to ACG1 tracers based on change in ACG1 and
c      yield rates for ARO ---
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRAG1)),nptcls(idxipt(ITRAG1))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delACG1 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO/AOA0 yield ---
c      SPECIAL CASE: non-volatile CG; ARO -> AOA0 directly (skip CG)
c
         iaoa0 = iptcls(idxipt(ITRPA0)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iaoa0 = iaoa0 + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(iaoa0)*yfac(1,1)
     &                      + ylrates(iaoa0)*yfac(2,1) )
         enddo
c
c  --- make adjustment to AOA0 tracers based on change in AOA0 and
c      yield rates for ARO ---
c      SPECIAL CASE: non-volatile CG; ARO -> AOA0 directly (skip CG)
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRPA0)),nptcls(idxipt(ITRPA0))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delNV1 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO ---
c
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            sumaro = sumaro + ptconc(idx)*wtkoh(i)
         enddo
c
c  --- make adjustment to ARO based on change in aromatic precursors and kOH ---
c
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delARO * 
     &                                  ptconc(idx)*wtkoh(i) / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ISP/TRP/API/SQT ---
c
         sumisp = 0.
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            sumisp = sumisp + ptconc(idx)
         enddo

         sumtrp = 0.
         do i=iptcls(idxipt(ITRTRP)),nptcls(idxipt(ITRTRP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            sumtrp = sumtrp + ptconc(idx)
         enddo

         sumapin = 0.
         do i=iptcls(idxipt(ITRAPN)),nptcls(idxipt(ITRAPN))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            sumapin = sumapin + ptconc(idx)
         enddo

         sumsqt = 0.
         do i=iptcls(idxipt(ITRSQT)),nptcls(idxipt(ITRSQT))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            sumsqt = sumsqt + ptconc(idx)
         enddo
c
c  --- make adjustment to BG2 tracers based on change in ISP/TRP/API/SQT and
c      distribution of ISP/TRP/API/SQT ----
c
         iisp  = iptcls(idxipt(ITRISP)) - 1
         itrp  = iptcls(idxipt(ITRTRP)) - 1
         iapin = iptcls(idxipt(ITRAPN)) - 1
         isqt  = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRBG2)),nptcls(idxipt(ITRBG2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(isqt-1)
            iapin = iapin + 1
            mdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iapin-1)
            ptconc(idx) = ptconc(idx) - delISP * yfac(1,6) *
     &                                           ptconc(jdx) / sumisp
     &                                - delTRP * yfac(1,7) *
     &                                           ptconc(kdx) / sumtrp
     &                                - delAPIN * yfac(1,8) *
     &                                           ptconc(mdx) / sumapin
     &                                - delSQT * yfac(1,9) *
     &                                           ptconc(ldx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to BG1 tracers based on change in ISP/TRP/API/SQT and
c      distribution of ISP/TRP/API/SQT ----
c
         iisp = iptcls(idxipt(ITRISP)) - 1
         itrp = iptcls(idxipt(ITRTRP)) - 1
         iapin = iptcls(idxipt(ITRAPN)) - 1
         isqt = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRBG1)),nptcls(idxipt(ITRBG1))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(isqt-1)
            iapin = iapin + 1
            mdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(mdx-1)
            ptconc(idx) = ptconc(idx) - delISP * yfac(2,6) *
     &                                           ptconc(jdx) / sumisp
     &                                - delTRP * yfac(2,7) *
     &                                           ptconc(kdx) / sumtrp
     &                                - delAPIN * yfac(2,8) *
     &                                           ptconc(ldx) / sumapin
     &                                - delSQT * yfac(2,9) *
     &                                           ptconc(ldx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to BOA0 tracers based on change in ISP/TRP/API/SQT and
c      distribution of ISP/TRP/API/SQT ----
c      SPECIAL CASE: non-volatile CG; ISP/TRP/API/SQT -> BOA0 directly (skip CG)
c
         iisp = iptcls(idxipt(ITRISP)) - 1
         itrp = iptcls(idxipt(ITRTRP)) - 1
         iapin = iptcls(idxipt(ITRAPN)) - 1
         isqt = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRPB0)),nptcls(idxipt(ITRPB0))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(isqt-1)
            iapin = iapin + 1
            mdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iapin-1)
            ptconc(idx) = ptconc(idx) - wtboa0 * 
     &                                ( delISP * yfac(3,6) *
     &                                           ptconc(jdx) / sumisp
     &                                + delTRP * yfac(3,7) *
     &                                           ptconc(kdx) / sumtrp
     &                                + delAPIN * yfac(3,8) *
     &                                           ptconc(ldx) / sumapin
     &                                + delSQT * yfac(3,9) *
     &                                           ptconc(ldx) / sumsqt )
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to ISP based on change in ISP ---
c
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delISP * ptconc(idx) / sumisp
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to TRP based on change in TRP ---
c
         do i=iptcls(idxipt(ITRTRP)),nptcls(idxipt(ITRTRP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delTRP * ptconc(idx) / sumtrp
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to APN based on change in APN ---
c
         do i=iptcls(idxipt(ITRAPN)),nptcls(idxipt(ITRAPN))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delAPIN * ptconc(idx) / sumAPIN
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to SQT based on change in SQT ---
c
         do i=iptcls(idxipt(ITRSQT)),nptcls(idxipt(ITRSQT))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delSQT * ptconc(idx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to IVA based on change in IVA ---
c
         sumIVOA = 0.
         do i=iptcls(idxipt(ITRIVA)),nptcls(idxipt(ITRIVA))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumIVOA = sumIVOA + ptconc(idx)
         enddo
         do i=iptcls(idxipt(ITRIVA)),nptcls(idxipt(ITRIVA))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delIVOA * ptconc(idx) / sumIVOA
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to SVA based on change in SVA ---
c
         sumSVOA = 0.
         do i=iptcls(idxipt(ITRSVA)),nptcls(idxipt(ITRSVA))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumSVOA = sumSVOA + ptconc(idx)
         enddo
         do i=iptcls(idxipt(ITRSVA)),nptcls(idxipt(ITRSVA))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delSVOA * ptconc(idx) / sumSVOA
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
      endif
c
c----------------------------------------------------------------------
c  ---  Sulfate tracers ---
c----------------------------------------------------------------------
c
      if( lsulfate ) then
         delso2 = delcon(3,kso2)
         delps4 = delcon(3,kpso4)
c
c  --- get the sum of the tracer species for SO2 ---
c
         sumso2 = 0.
         do i=iptcls(idxipt(ITRSO2)),nptcls(idxipt(ITRSO2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumso2 = sumso2 + ptconc(idx)
         enddo
c
c  --- get the sum of the tracer species for DMS ---
c
         if( ldmschm ) then
            idms = iptcls(idxipt(ITRDMS)) - 1
            deldms = delcon(1,kdms)
            sumdms = 0.
            do i=iptcls(idxipt(ITRDMS)),nptcls(idxipt(ITRDMS))
               idx = idxcel + numcol * numrow * numlay * (i-1)
               sumdms = sumdms + ptconc(idx)
            enddo
         endif
c
c  --- if PSO4 change is positive, make adjustment to PS4 
c      tracers based on distribtion of SO2 ---
c
         if( delps4 .GT. 0. ) then
            iso2 = iptcls(idxipt(ITRSO2)) - 1
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
               iso2 = iso2 + 1
               jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iso2-1)
               ptconc(idx) = ptconc(idx) + delps4 * ptconc(jdx) / sumso2
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
            enddo
c
c  --- if PSO4 change is negative, make adjustment to PS4 
c      tracers based on distribtion of PS4 ---
c
         else
c
c  --- get the sum of the tracer species for PS4 ---
c
            sumps4 = 0.
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
               sumps4 = sumps4 + ptconc(idx)
            enddo
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
               ptconc(idx) = ptconc(idx) + delps4 * ptconc(idx) / sumps4
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
            enddo
         endif
c
c  --- make adjustment to SO2 tracers based on change in SO2 (and DMS) ---
c
         do i=iptcls(idxipt(ITRSO2)),nptcls(idxipt(ITRSO2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delso2 * ptconc(idx) / sumso2
            if( ldmschm ) then
               idms = idms + 1
               jdx = idxcel + numcol * numrow * numlay * (idms-1)
               ptconc(idx) = ptconc(idx) - deldms * ptconc(jdx) / sumdms
            endif
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to DMS tracers based on change in DMS ---
c
         if( ldmschm ) then
           do i=iptcls(idxipt(ITRDMS)),nptcls(idxipt(ITRDMS))
              idx = idxcel + numcol * numrow * numlay * (i-1)
              ptconc(idx) = ptconc(idx) + deldms * ptconc(idx) / sumdms
              ptconc(idx) = MAX(BNDLPT,ptconc(idx))
           enddo
         endif
      endif
c
c----------------------------------------------------------------------
c  ---  Mercury tracers ----
c----------------------------------------------------------------------
c
      if( lmercury ) then
         delhg2 = delcon(3,khg2)

c  --- get the sum of the tracer species ---
c
         sumhg0 = 0.
         do i=iptcls(idxipt(ITRHG0)),nptcls(idxipt(ITRHG0))
            idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumhg0 = sumhg0 + ptconc(idx0)
         enddo
         sumhg2 = 0.
         do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
            idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumhg2 = sumhg2 + ptconc(idx2)
         enddo
c
c  --- HG0 being oxidized to HG2 ----
c
         if( delhg2 .GT. 0 ) then
            ihg0 = iptcls(idxipt(ITRHG0)) - 1
            do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
               idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
               ihg0 = ihg0 + 1
               idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(ihg0-1)
               ptconc(idx2) = ptconc(idx2) + 
     &                                   delhg2 * ptconc(idx0) / sumhg0
               ptconc(idx2) = MAX(BNDLPT,ptconc(idx2))
               ptconc(idx0) = ptconc(idx0) - 
     &                                   delhg2 * ptconc(idx2) / sumhg2
               ptconc(idx0) = MAX(BNDLPT,ptconc(idx0))
            enddo
c
c  --- HG2 beign reduced to HG0 ----
c
         else
            ihg0 = iptcls(idxipt(ITRHG0)) - 1
            do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
               idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                      DBLE(numlay) * DBLE(i-1)
               ihg0 = ihg0 + 1
               idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(ihg0-1)
               ptconc(idx2) = ptconc(idx2) + 
     &                                  delhg2 * ptconc(idx2) / sumhg2
               ptconc(idx2) = MAX(BNDLPT,ptconc(idx2))
               ptconc(idx0) = ptconc(idx0) - 
     &                                  delhg2 * ptconc(idx0) / sumhg0
               ptconc(idx0) = MAX(BNDLPT,ptconc(idx0))
            enddo
         endif
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
