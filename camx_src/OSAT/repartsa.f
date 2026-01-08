c*** REPARTSA.F
c
      subroutine repartsa(numcol,numrow,numlay,igrid,
     &                                  icell,jcell,kcell,modcon,delcon)
      use grid
      use chmstry
      use tracer
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine repartitions the particulate and gaseous portions
c     nitric acid and ammonia tracers based on the ratio found in the
c     regular model.
c
c      Copyright 1996 - 2025
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c       numcol  I number of columns in this slice
c       numrow  I number of rows in this slice
c       numlay  I number of layers in this slice
c       igrid   I  grid number
c       icell   I  the X grid location of current cell
c       jcell   I  the X grid location of current cell
c       kcell   I  the vertical grid location of current layer
c       modcon  R  array of model species concentrations
c       delcon  R  array of change in model species concentrations
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     09/28/03   --gwilson--  Original development
c     12/29/06   --bkoo--     Revised for the updated SOA scheme
c     05/05/07   --gwilson--  Added code to make an adjustment to
c                             tracer species if there is stray from 
c                             the regular - added to handle mass errors
c                             in regular model chemistry routines
c     08/25/16   --bkoo--     Updated for new SOAP
c     09/16/16   --bkoo--     Updated for in-cloud SOA formation
c     10/16/17   --bkoo--     Updated for SOA photolysis
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
      integer numcol
      integer numrow
      integer numlay
      integer igrid
      integer icell
      integer jcell
      integer kcell
      real    modcon(*)
      real    delcon(7,MXSPEC+1)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c    FUZZ - fuzz value for comparing against lower bound
c
      real FUZZ
c
      parameter( FUZZ = 10.0 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer*8 idxcel, idxhno3, idxpno3, idxnh3, idxpnh4
      integer*8 idxpo1, idxpo2, idxpo3, idxpo4
      integer*8 idxacg2, idxacg1, idxbcg2, idxbcg1
      integer*8 idxppa, idxppb, idxisp
      integer   i, ipno3, ipnh4, ipo1, ipo2, ipo3, ipo4
      integer   ippa, ippb
      real      wthno3, wtnh3, conhno3, conpno3, connh3, conpnh4
      real      conacg2, conacg1, conbcg2, conbcg1
      real      conpo1, conpo2, conpo3, conpo4, conppa, conppb
      real      wtacg2, wtacg1, wtbcg2, wtbcg1
      real      ratio, sumtrac, contrac, conadjust
      real      pratio1, pratio2
      real      sumisp, delppb
c
c-----------------------------------------------------------------------
c    Data statements:
c     Note: MWs for CG species should be consistent with those defined
c           in SOAP.INC
c-----------------------------------------------------------------------
c
      data wtnh3  /18.0/
      data wthno3 /62.0/
      data wtacg2 /150.0/
      data wtacg1 /150.0/
      data wtbcg2 /180.0/
      data wtbcg1 /180.0/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  DBLE(ipsa3d(igrid)-1) + DBLE(icell) + 
     &    DBLE(numcol) * DBLE(jcell-1) + DBLE(numcol * numrow *(kcell-1))
c
      if( lnitrate ) then
c
c  --- calculate the ratio of gas to  particluate nitric acid ---
c
         conhno3 = modcon(khno3) * wthno3
         conpno3 = modcon(kpno3)
         ratio = 0.
         if( conhno3+conpno3 .GT. 0. )
     &                    ratio = conhno3 / (conpno3 + conhno3)
c
c  --- get the sum of all tracer species ---
c
         contrac = 0.0
         ipno3 = iptcls(idxipt(ITRPN3)) - 1
         do i=iptcls(idxipt(ITRHN3)),nptcls(idxipt(ITRHN3))
            idxhno3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ipno3 = ipno3 + 1
            idxpno3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(ipno3-1)
            sumtrac = ptconc(idxhno3)*wthno3 + ptconc(idxpno3)
            if( sumtrac .GT. FUZZ*BNDLPT ) contrac = contrac + sumtrac
         enddo
         if( contrac .GT. 0. ) then
             conadjust = (conhno3+conpno3) / contrac
         else
             conadjust = 1.0
         endif
c
c  --- adjust the nitric acid tracers ----
c
         ipno3 = iptcls(idxipt(ITRPN3)) - 1
         do i=iptcls(idxipt(ITRHN3)),nptcls(idxipt(ITRHN3))
            idxhno3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ipno3 = ipno3 + 1
            idxpno3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(ipno3-1)
c
c  --- adjust for mass imbalance, if needed ---
c
            sumtrac = ptconc(idxhno3)*wthno3 + ptconc(idxpno3)
            if( sumtrac .GT. FUZZ*BNDLPT ) sumtrac = conadjust * sumtrac
c
c  --- redistribute the gas and particulates ---
c
            ptconc(idxhno3) = (ratio * sumtrac ) / wthno3
            ptconc(idxhno3) = MAX(ptconc(idxhno3),BNDLPT)
            ptconc(idxpno3) = sumtrac - ptconc(idxhno3) * wthno3
            ptconc(idxpno3) = MAX(ptconc(idxpno3),BNDLPT)
         enddo
c
c  --- calculate the ratio of gas to particluate ammonia ---
c
         connh3 = modcon(knh3) * wtnh3
         conpnh4 = modcon(kpnh4)
         ratio = 0.
         if( connh3+conpnh4 .GT. 0. )
     &                    ratio = connh3 / (conpnh4 + connh3)
c
c  --- get the sum of all tracer species ---
c
         contrac = 0.0
         ipnh4 = iptcls(idxipt(ITRPN4)) - 1
         do i=iptcls(idxipt(ITRNH3)),nptcls(idxipt(ITRNH3))
            idxnh3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            ipnh4 = ipnh4 + 1
            idxpnh4 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(ipnh4-1)
            sumtrac = ptconc(idxnh3)*wtnh3 + ptconc(idxpnh4)
            if( sumtrac .GT. FUZZ*BNDLPT ) contrac = contrac + sumtrac
         enddo
         if( contrac .GT. 0. ) then
             conadjust = (connh3+conpnh4) / contrac
         else
             conadjust = 1.0
         endif
c
c  --- adjust the ammonia tracers ----
c
         ipnh4 = iptcls(idxipt(ITRPN4)) - 1
         do i=iptcls(idxipt(ITRNH3)),nptcls(idxipt(ITRNH3))
            idxnh3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
            ipnh4 = ipnh4 + 1
            idxpnh4 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(ipnh4-1)	
c        
c  --- adjust for mass imbalance, if needed ---
c  
            sumtrac = ptconc(idxnh3)*wtnh3 + ptconc(idxpnh4)
            if( sumtrac .GT. FUZZ*BNDLPT ) sumtrac = conadjust * sumtrac
c
c  --- redistribute the gas and particulates ---
c
            ptconc(idxnh3) = ratio * sumtrac / wtnh3
            ptconc(idxnh3) = MAX(ptconc(idxnh3),BNDLPT)
            ptconc(idxpnh4) = sumtrac - ptconc(idxnh3) * wtnh3
            ptconc(idxpnh4) = MAX(ptconc(idxpnh4),BNDLPT)
         enddo
      endif
c
c  --- secondary aerosols ---
c
      if( lsoa ) then
c
c  --- apportion photolyzed AOA0 & BOA0 first
c
         conppa = modcon(kaoa0) - delcon(5,kaoa0) ! get equil con by reverting the changes
     &                          - delcon(6,kaoa0) ! due to SOA polymerization & photolysis
         pratio1 = 0.
         if ( conppa .GT. 0. ) pratio1 = -delcon(6,kaoa0) / conppa

         do i=iptcls(idxipt(ITRPA0)),nptcls(idxipt(ITRPA0))
           idxppa = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
           ptconc(idxppa) = ptconc(idxppa) - pratio1*ptconc(idxppa)
         enddo

         conppb = modcon(kboa0) - delcon(5,kboa0) ! get equil con by reverting the changes
     &                          - delcon(6,kboa0) ! due to SOA polymerization, photolysis,
     &                          - delcon(7,kboa0) ! and in-cloud SOA formation
         pratio1 = 0.
         if ( conppb .GT. 0. ) pratio1 = -delcon(6,kboa0) / conppb

         do i=iptcls(idxipt(ITRPB0)),nptcls(idxipt(ITRPB0))
           idxppb = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
           ptconc(idxppb) = ptconc(idxppb) - pratio1*ptconc(idxppb)
         enddo
c
c  --- calculate the ratio of AG2 to PA2 ---
c
         conacg2 = modcon(kacg2) * wtacg2
         conpo1 = modcon(kaoa2) - delcon(5,kaoa2) ! get equil con by reverting the changes
     &                          - delcon(6,kaoa2) ! due to SOA polymerization & photolysis
         ratio = 0.
         if( conacg2+conpo1 .GT. 0. )
     &                    ratio = conacg2 / (conpo1 + conacg2)
         pratio1 = 0.
         if ( conpo1 .GT. 0. ) pratio1 = -delcon(6,kaoa2) / conpo1
         conpo1 = conpo1 + delcon(6,kaoa2)
         pratio2 = 0.
         if ( conpo1 .GT. 0. ) pratio2 = -delcon(5,kaoa2) / conpo1
c
c  --- adjust the AG2/PA2 ratio ---
c
         ipo1 = iptcls(idxipt(ITRPA2)) - 1
         ippa = iptcls(idxipt(ITRPA0)) - 1
         do i=iptcls(idxipt(ITRAG2)),nptcls(idxipt(ITRAG2))
            idxacg2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
            ipo1 = ipo1 + 1
            ippa = ippa + 1
            idxpo1 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(ipo1-1)
            idxppa = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(ippa-1)
            sumtrac = ptconc(idxacg2)*wtacg2 + ptconc(idxpo1)
            ptconc(idxacg2) = (ratio * sumtrac ) / wtacg2              ! by partitioning
            ptconc(idxpo1) = (1.0-ratio) * sumtrac                   ! by partitioning
            ptconc(idxpo1) = (1.0-pratio1) * ptconc(idxpo1)          ! by photolysis
            ptconc(idxppa) = ptconc(idxppa) + pratio2*ptconc(idxpo1) ! by polymerization
            ptconc(idxpo1) = (1.0-pratio2) * ptconc(idxpo1)          ! by polymerization
            ptconc(idxacg2) = MAX(ptconc(idxacg2),BNDLPT)
            ptconc(idxpo1) = MAX(ptconc(idxpo1),BNDLPT)
            ptconc(idxppa) = MAX(ptconc(idxppa),BNDLPT)
         enddo
c
c  --- calculate the ratio of ACG1 to AOA1 ---
c
         conacg1 = modcon(kacg1) * wtacg1
         conpo2 = modcon(kaoa1) - delcon(5,kaoa1) ! get equil con by reverting the changes
     &                          - delcon(6,kaoa1) ! due to SOA polymerization & photolysis
         ratio = 0.
         if( conacg1+conpo2 .GT. 0. )
     &                    ratio = conacg1 / (conpo2 + conacg1)
         pratio1 = 0.
         if ( conpo2 .GT. 0. ) pratio1 = -delcon(6,kaoa1) / conpo2
         conpo2 = conpo2 + delcon(6,kaoa1)
         pratio2 = 0.
         if ( conpo2 .GT. 0. ) pratio2 = -delcon(5,kaoa1) / conpo2
c
c  --- adjust the ACG1/AOA1 ratio ---
c
         ipo2 = iptcls(idxipt(ITRPA1)) - 1
         ippa = iptcls(idxipt(ITRPA0)) - 1
         do i=iptcls(idxipt(ITRAG1)),nptcls(idxipt(ITRAG1))
            idxacg1 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            ipo2 = ipo2 + 1
            ippa = ippa + 1
            idxpo2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(ipo2-1)
            idxppa = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(ippa-1)
            sumtrac = ptconc(idxacg2)*wtacg1 + ptconc(idxpo2)
            ptconc(idxacg1) = (ratio * sumtrac ) / wtacg1              ! by partitioning
            ptconc(idxpo2) = (1.0-ratio) * sumtrac                   ! by partitioning
            ptconc(idxpo2) = (1.0-pratio1) * ptconc(idxpo2)          ! by photolysis
            ptconc(idxppa) = ptconc(idxppa) + pratio2*ptconc(idxpo2) ! by polymerization
            ptconc(idxpo2) = (1.0-pratio2) * ptconc(idxpo2)          ! by polymerization
            ptconc(idxacg1) = MAX(ptconc(idxacg1),BNDLPT)
            ptconc(idxpo2) = MAX(ptconc(idxpo2),BNDLPT)
            ptconc(idxppa) = MAX(ptconc(idxppa),BNDLPT)
         enddo
c
c  --- calculate the ratio of BG2 to PB2 ---
c
         conbcg2 = modcon(kbcg2) * wtbcg2
         conpo3 = modcon(kboa2) - delcon(5,kboa2) ! get equil con by reverting the changes
     &                          - delcon(6,kboa2) ! due to SOA polymerization & photolysis
         ratio = 0.
         if( conbcg2+conpo3 .GT. 0. )
     &                    ratio = conbcg2 / (conpo3 + conbcg2)
         pratio1 = 0.
         if ( conpo3 .GT. 0. ) pratio1 = -delcon(6,kboa2) / conpo3
         conpo3 = conpo3 + delcon(6,kboa2)
         pratio2 = 0.
         if ( conpo3 .GT. 0. ) pratio2 = -delcon(5,kboa2) / conpo3
c
c  --- adjust the BG2/PB2 ratio ---
c
         ipo3 = iptcls(idxipt(ITRPB2)) - 1
         ippb = iptcls(idxipt(ITRPB0)) - 1
         do i=iptcls(idxipt(ITRBG2)),nptcls(idxipt(ITRBG2))
            idxBCG2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            ipo3 = ipo3 + 1
            ippb = ippb + 1
            idxpo3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(ipo3-1)
            idxppb = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(ippb-1)
            sumtrac = ptconc(idxbcg2)*wtbcg2 + ptconc(idxpo3)
            ptconc(idxbcg2) = (ratio * sumtrac ) / wtbcg2              ! by partitioning
            ptconc(idxpo3) = (1.0-ratio) * sumtrac                   ! by partitioning
            ptconc(idxpo3) = (1.0-pratio1) * ptconc(idxpo3)          ! by photolysis
            ptconc(idxppb) = ptconc(idxppb) + pratio2*ptconc(idxpo3) ! by polymerization
            ptconc(idxpo3) = (1.0-pratio2) * ptconc(idxpo3)          ! by polymerization
            ptconc(idxbcg2) = MAX(ptconc(idxbcg2),BNDLPT)
            ptconc(idxpo3) = MAX(ptconc(idxpo3),BNDLPT)
            ptconc(idxppb) = MAX(ptconc(idxppb),BNDLPT)
         enddo
c
c  --- calculate the ratio of BG1 to PB1 ---
c
         conbcg1 = modcon(kbcg1) * wtbcg1
         conpo4 = modcon(kboa1) - delcon(5,kboa1) ! get equil con by reverting the changes
     &                          - delcon(6,kboa1) ! due to SOA polymerization & photolysis
         ratio = 0.
         if( conbcg1+conpo4 .GT. 0. )
     &                    ratio = conbcg1 / (conpo4 + conbcg1)
         pratio1 = 0.
         if ( conpo4 .GT. 0. ) pratio1 = -delcon(6,kboa1) / conpo4
         conpo4 = conpo4 + delcon(6,kboa1)
         pratio2 = 0.
         if ( conpo4 .GT. 0. ) pratio2 = -delcon(5,kboa1) / conpo4
c
c  --- adjust the BG1/PB1 ratio ---
c
         ipo4 = iptcls(idxipt(ITRPB1)) - 1
         ippb = iptcls(idxipt(ITRPB0)) - 1
         do i=iptcls(idxipt(ITRBG1)),nptcls(idxipt(ITRBG1))
            idxbcg1 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            ipo4 = ipo4 + 1
            ippb = ippb + 1
            idxpo4 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(ipo4-1)
            idxppb = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(ippb-1)
            sumtrac = ptconc(idxbcg1)*wtbcg1 + ptconc(idxpo4)
            ptconc(idxbcg1) = (ratio * sumtrac ) / wtbcg1              ! by partitioning
            ptconc(idxpo4) = (1.0-ratio) * sumtrac                   ! by partitioning
            ptconc(idxpo4) = (1.0-pratio1) * ptconc(idxpo4)          ! by photolysis
            ptconc(idxppb) = ptconc(idxppb) + pratio2*ptconc(idxpo4) ! by polymerization
            ptconc(idxpo4) = (1.0-pratio2) * ptconc(idxpo4)          ! by polymerization
            ptconc(idxbcg1) = MAX(ptconc(idxbcg1),BNDLPT)
            ptconc(idxpo4) = MAX(ptconc(idxpo4),BNDLPT)
            ptconc(idxppb) = MAX(ptconc(idxppb),BNDLPT)
         enddo
c
c  --- get the sum of the tracer species for ISP ---
c
         sumisp = 0.
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idxisp = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            sumisp = sumisp + ptconc(idxisp)
         enddo
c
c  --- adjust PB0 for in-cloud SOA formation based on distribution of ISP ---
c
         delppb = delcon(7,kboa0)
         ippb = iptcls(idxipt(ITRPB0)) - 1
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idxisp = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            ippb = ippb + 1
            idxppb = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(ippb-1)
            ptconc(idxppb) = ptconc(idxppb) + delppb * ptconc(idxisp) / sumisp
            ptconc(idxppb) = MAX(ptconc(idxppb),BNDLPT)
         enddo
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
