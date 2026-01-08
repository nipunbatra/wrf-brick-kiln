c**** STABSA
c
      subroutine stabsa(idate,begtim,jdate,endtim)
      use filunit
      use chmstry
      use tracer
      use rtracchm
      use grid
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the flags that determine if each species in the
c   species list should be added to the VOC or NOx tracers.  The carbon
c   number for the VOC species is also put into the appropriate place
c   in the array.
c      Inputs:
c        idate  I   date of the beginning of the simulation (YYJJJ)
c        begtim R   hour of the begining of simulation
c        jdate  I   date of the ending of the simulation (YYJJJ)
c        endtim R   hour of the endng of simulation
c
c      Copyright 1996 - 2025
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/27/96   --gwilson--    Original development
c     09/20/03   --gwilson--    Completely changed the way tracer classes
c                               are identified to deal with PSAT
c     08/20/06   --bkoo--       Added ETOH, MTBE & MBUT for updated SAPRC99
c                               Fixed rkohsprc & rmirsprc
c     09/01/06   --bkoo--       Updated rkohcbiv for ISOP, MEOH & ETOH
c     12/29/06   --bkoo--       Revised for the updated SOA scheme
c     01/08/07   --bkoo--       Added Mechanism 6 (CB05)
c     07/16/07   --bkoo--       Added HRVOC
c     10/25/07   --gyarwood-    Make MIR and kOH values per carbon
c     10/25/09   --gwilson--    Fixed coded related to timing tracers
c     08/23/13   --bkoo--       Set LVOCSP flags of SOA precursors for DDM-PM
c     03/18/14   --bkoo--       Added BNZA
c     06/24/16   --bkoo--       Updated Mechanism 4 (CB6r4)
c     08/25/16   --bkoo--       Updated for new SOAP
c     11/27/17   --bkoo--       Updated for CB6r2h
c     01/12/18   --bkoo--       Removed BNZA/TOLA/XYLA/ISP/TRP
c     01/16/18   --bkoo--       Updated for SAPRC07T
c     06/18/19   --bkoo--       Expanded spcsoap & rkohsoap for SAPRC07T
c     02/21/20   --rlb--        Added Mechanism 1 (CB6r5)
c     03/25/21   --gy--         Added Mechanism 7 (CB7)
c     09/02/21   --gy--         Added Mechanism 9 (RACM2s21)
c     11/03/23   --cemery--     Removed CB05 (Mech 6)
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'soap.inc'
      include 'soap3.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      integer   jdate
      real      begtim
      real      endtim
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   NSASPRC    I    the number of VOC species names to be checked (SAPRC)
c   NSARC2     I    the number of VOC species names to be checked (RACM2)
c   NSACB6     I    the number of VOC species names to be checked (CB6)
c   NSACB7     I    the number of VOC species names to be checked (CB7)
c   NSASOAP    I    the number of SOA-producing VOC species for PSAT
c   MXSANAM    I    maximum size of the arrays
c   NHRVSPRC   I    the number of HRVOC species names (SAPRC)
c   NHRVRC2    I    the number of HRVOC species names (RACM2)
c   NHRVCB6    I    the number of HRVOC species names (CB6)
c   NHRVCB7    I    the number of HRVOC species names (CB7)
c   MXHRVNAM   I    maximum size of the HRVOC arrays
c
      integer   NSASPRC, NSACB6, NSACB7, NSARC2, MXSANAM
      integer   NSASOAP
      integer   NHRVSPRC, NHRVCB6, NHRVCB7, NHRVRC2, MXHRVNAM
c
      parameter( NSASPRC = 30 )
      parameter( NSARC2  = 27 )
      parameter( NSACB6  = 19 )
      parameter( NSACB7  = 21 )
      parameter( NSASOAP = 12 ) ! Only the precursors used for calculating
                                ! reactivity/yield weighting factors
      parameter( MXSANAM = MAX(NSASPRC,
     &                         NSACB6, NSACB7, NSARC2) )
      parameter( NHRVSPRC = 5 )
      parameter( NHRVRC2  = 4 )
      parameter( NHRVCB6  = 3 )
      parameter( NHRVCB7  = 3 )
      parameter( MXHRVNAM = MAX(NHRVSPRC,
     &                          NHRVCB6, NHRVCB7, NHRVRC2) )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 spcnam(MXSANAM), hrvnam(MXHRVNAM)
      character*10 spcsprc(NSASPRC), spcrc2(NSARC2),
     &             spccb6(NSACB6), spccb7(NSACB7)
      character*10 hrvsprc(NHRVSPRC), hrvrc2(NHRVRC2),
     &             hrvcb6(NHRVCB6), hrvcb7(NHRVCB7)
      character*10 spcsoap(NSASOAP)
      character*10 nameyld(MXALCLS,MXSPEC)
      character*10 nonam, no2nam, hononam, o3name, namecls(MXALCLS,MXSPEC)
      integer      numcls(MXALCLS), ispc, nnames, i
      integer      idxo3n, idxo3v, idxoon, idxoov
      integer      numyld(MXALCLS), nhrvoc, mvecedge
      real         carbon(MXSANAM), rkoh(MXSANAM), rmir(MXSANAM)
      real         carbsprc(NSASPRC),rkohsprc(NSASPRC),rmirsprc(NSASPRC)
      real         carbrc2(NSARC2),  rkohrc2(NSARC2),  rmirrc2(NSARC2)
      real         carbcb6(NSACB6),  rkohcb6(NSACB6),  rmircb6(NSACB6)
      real         carbcb7(NSACB7),  rkohcb7(NSACB7),  rmircb7(NSACB7)
      real         rkohsoap(NSASOAP)
      real         coefcon(MXALCLS,MXSPEC), coeflux(MXALCLS,MXSPEC)
      real         yieldH(MXALCLS,MXSPEC),yieldL(MXALCLS,MXSPEC)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data spcrc2  /'HCHO      ', 'ACD       ', 'ALD       ',
     &              'ACT       ', 'MEK       ', 'KET       ',
     &              'MOH       ', 'EOH       ', 'ROH       ',
     &              'ETEG      ', 'ISO       ', 'API       ',
     &              'LIM       ', 'BEN       ', 'TOL       ',
     &              'XYM       ', 'XYP       ', 'XYO       ',
     &              'ETH       ', 'HC3       ', 'HC5       ',
     &              'HC8       ', 'ETE       ', 'OLT       ',
     &              'OLI       ', 'ACE       ', 'DIEN      ' /
c
      data carbrc2  /1.0,          2.0,          3.0,
     &               3.0,          4.0,          5.0,
     &               1.0,          2.0,          3.0,
     &               2.0,          5.0,          10.0,
     &               10.0,         6.0,          7.1,
     &               8.9,          9.0,          8.9, 
     &               2.0,          3.6,          5.6,
     &               7.9,          2.0,          3.8,
     &               5.0,          2.0,          4.0 /
c
      data rkohrc2  /1.24E+04,     2.21E+04,     2.82E+04,
     &               2.05E+02,     1.64E+03,     4.28E+03,
     &               1.32E+03,     4.74E+03,     7.51E+03,
     &               2.17E+04,     1.48E+05,     7.82E+04,
     &               2.38E+05,     1.80E+03,     8.77E+03,
     &               3.41E+04,     2.11E+04,     2.01E+04,
     &               3.56E+02,     3.28E+03,     6.56E+03,
     &               1.67E+04,     1.21E+04,     4.52E+04,
     &               1.05E+05,     1.10E+03,     9.83E+04 /
c
      data rmirrc2  /3.41,         4.46,         6.67,
     &               0.47,         2.61,         4.97,
     &               0.36,         1.34,         2.18,
     &               4.62,         12.0,         12.7,
     &               14.9,         1.12,         5.23,
     &               21.6,         18.0,         17.7,
     &               0.18,         1.01,         1.05,
     &               3.46,         4.08,         5.36,
     &               13.4,         0.28,         8.67 /
c
      data hrvrc2  /'ETE       ','OLT       ','OLI       ',
     &              'DIEN      ' /
c
      data spcsprc /'HCHO      ','CCHO      ','RCHO      ',
     &              'ACET      ','MEK       ','MEOH      ',
     &              'ETHE      ','PRPE      ','BD13      ',
     &              'ISOP      ','APIN      ','ACYE      ',
     &              'BENZ      ','TOLU      ','MXYL      ',
     &              'OXYL      ','PXYL      ','B124      ',
     &              'ETOH      ','ALK1      ','ALK2      ',
     &              'ALK3      ','ALK4      ','ALK5      ',
     &              'OLE1      ','OLE2      ','ARO1      ',
     &              'ARO2      ','TERP      ','SESQ      '/
c
      data carbsprc /1.0,         2.0,         3.0,
     &               3.0,         4.0,         1.0,
     &               2.0,         3.0,         4.0,
     &               5.0,        10.0,         2.0,
     &               6.0,         7.0,         8.0,
     &               8.0,         8.0,         9.0,
     &               2.0,         2.0,         2.5,
     &               4.0,         5.4,         8.3,
     &               5.2,         5.4,         7.2,
     &               8.8,        10.0,        15.0/
c
      data rkohsprc  /1.25E4,     2.21E4,     2.93E4,
     &                2.77E2,     1.74E3,     1.32E3,
     &                1.21E4,     3.89E4,     9.83E4,
     &                1.49E5,     7.72E4,     1.12E3,
     &                1.80E3,     8.31E3,     3.41E4,
     &                2.01E4,     2.11E4,     4.80E4,
     &                4.74E3,     3.66E2,     1.62E3,
     &                3.40E3,     6.42E3,     1.40E4,
     &                5.33E4,     9.57E4,     1.16E4,
     &                4.56E4,     1.44E5,     1.44E5/
c
      data rmirsprc  /5.92,       6.00,       8.57,
     &                0.43,       2.23,       0.45,
     &                5.26,      10.23,      14.21,
     &               15.05,      12.79,       0.52,
     &                1.17,       7.69,      21.57,
     &               16.90,      12.92,      22.22,
     &                1.46,       0.18,       0.45,
     &                1.44,       4.77,       3.01,
     &               10.53,      15.30,       5.97,
     &               20.88,      11.48,       5.74/
c
      data hrvsprc /'ETHE      ','PRPE      ','BD13      ',
     &              'OLE1      ','OLE2      '/
c
      data spccb6  /'PAR       ','ETHA      ','MEOH      ',
     &              'ETOH      ','ETH       ','OLE       ',
     &              'IOLE      ','ISOP      ','TERP      ',
     &              'FORM      ','ALD2      ','ALDX      ',
     &              'TOL       ','XYL       ','PRPA      ',
     &              'BENZ      ','ETHY      ','ACET      ',
     &              'KET       '/
c
      data carbcb6  /1.,          2.,          1.,
     &               2.,          2.,          2.,
     &               4.,          5.,         10.,
     &               1.,          2.,          2.,
     &               7.,          8.,          3.,
     &               6.,          2.,          3.,
     &               1./
c
      data rkohcb6  /1.20E3,      3.56E2,      1.32E3,
     &               4.74E3,      1.16E4,      4.22E4,
     &               8.85E4,      1.48E5,      1.00E5,
     &               1.25E4,      2.21E4,      2.82E4,
     &               8.32E3,      2.73E4,      1.58E3,
     &               1.80E3,      1.11E3,      2.60E2,
     &               0.00E0/
c
      data rmircb6  /0.509,       0.135,       0.480,
     &               1.53,        4.95,        9.66,
     &               16.0,        12.7,        9.91,
     &               4.87,        5.80,        8.35,
     &               7.39,        20.5,        0.541,
     &               1.39,        0.487,       0.564,
     &               4.86/
c
      data hrvcb6  /'ETH       ','IOLE      ','OLE       '/
c
      data spccb7  /'PAR       ','ETHA      ','MEOH      ',
     &              'ETOH      ','ETH       ','OLE       ',
     &              'IOLE      ','ISOP      ','TERP      ',
     &              'FORM      ','ALD2      ','ALDX      ',
     &              'TOL       ','XYL       ','PRPA      ',
     &              'BENZ      ','ETHY      ','ACET      ',
     &              'KET       ','APIN      ','SQT       '/
c
      data carbcb7  /1.,          2.,          1.,
     &               2.,          2.,          2.,
     &               4.,          5.,         10.,
     &               1.,          2.,          2.,
     &               7.,          8.,          3.,
     &               6.,          2.,          3.,
     &               1.,         10.,         15./
c
      data rkohcb7  /1.23E3,      3.56E2,      1.32E3,
     &               4.74E3,      1.16E4,      4.22E4,
     &               8.85E4,      1.48E5,      1.95E5,
     &               1.25E4,      2.21E4,      2.82E4,
     &               8.32E3,      2.73E4,      1.58E3,
     &               1.80E3,      1.11E3,      2.60E2,
     &               1.48E3,      4.04E4,      2.96E5/
c
      data rmircb7  /0.36,        0.155,       0.353,
     &               1.33,        3.83,        7.61,
     &               12.9,        10.5,        11.0,
     &               3.86,        4.93,        5.65,
     &               6.52,        16.4,        0.38,
     &               1.02,        0.41,        0.51,
     &               0.94,        9.31,        5.87/
c
      data hrvcb7  /'ETH       ','IOLE      ','OLE       '/
c
      data spcsoap /'BENZ      ',
     &              'TOL       ','TOLU      ','ARO1      ',
     &              'XYL       ','OXYL      ','MXYL      ',
     &              'PXYL      ','B124      ','ARO2      ',
     &              'IVOA      ','SVOA      '/
c
      data rkohsoap / 1.80E3,
     &                8.32E3,     8.32E3,      8.32E3,
     &                2.73E4,     2.73E4,      2.73E4,
     &                2.73E4,     2.73E4,      2.73E4,
     &                1.98E4,     1.98E4/
c
      data nonam  /'NO        '/
c
      data no2nam /'NO2       '/
c
      data hononam /'HONO      '/
c
      data o3name /'O3        '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- initialize the arrays for reactivity ---
c
      do i=1,nspec
         crbnum(i) = 0.
         rkohrt(i) = 0.
         rmirrt(i) = 0.
         lusespc(i) = .FALSE.
      enddo
c
c  --- calculate the beginning of the various tracer types ---
c      there will be (ngroup+1) if there is an extra group for the
c      "leftover" group  ----
c      
      if( lsa_ioric ) then
          numic_spcs = ncls_ioric
      else
          numic_spcs = 1
      endif
      if( lsa_iorbc ) then
          numbc_spcs = ncls_iorbc
          numtc_spcs = 0
          if( .NOT. lsa_iortc ) numtc_spcs = 1
      else
          numtc_spcs = 0
          if( lbndry ) then
              numbc_spcs = 5
          else
              numbc_spcs = 1
          endif
      endif
      nbdic = numic_spcs + numbc_spcs + numtc_spcs 
      if( ngroup .EQ. 0 ) then
          ncount = nregin 
      else
          if( leftovr_area ) then
             ncount = (ngroup + 1) * nregin
          else
             ncount = ngroup * nregin
          endif
      endif
c
c  ---- load the local arrays based on the mechanism ---
c
      if( idmech .EQ. 1 .OR. idmech .EQ. 4) then
         nnames = NSACB6
         call spccb6r4sa(namecls,numcls,coefcon,coeflux,
     &                      nameyld,numyld,yieldH,yieldL,mwspec)
         do i=1,nnames
            spcnam(i) = spccb6(i)
            carbon(i) = carbcb6(i)
            rkoh(i) = rkohcb6(i)/carbon(i)
            rmir(i) = rmircb6(i)/carbon(i)
         enddo
         nhrvoc = NHRVCB6
         do i=1,nhrvoc
            hrvnam(i) = hrvcb6(i)
         enddo
      elseif( idmech .EQ. 3 ) then
         nnames = NSACB6
         call spccb6r5hsa(namecls,numcls,coefcon,coeflux,
     &                      nameyld,numyld,yieldH,yieldL,mwspec)
         do i=1,nnames
            spcnam(i) = spccb6(i)
            carbon(i) = carbcb6(i)
            rkoh(i) = rkohcb6(i)/carbon(i)
            rmir(i) = rmircb6(i)/carbon(i)
         enddo
         nhrvoc = NHRVCB6
         do i=1,nhrvoc
            hrvnam(i) = hrvcb6(i)
         enddo
      elseif( idmech .EQ. 5 ) then
         nnames = NSASPRC
         call spcsprcsa(namecls,numcls,coefcon,coeflux,
     &                  nameyld,numyld,yieldH,yieldL,mwspec)
         do i=1,nnames
            spcnam(i) = spcsprc(i)
            carbon(i) = carbsprc(i)
            rkoh(i) = rkohsprc(i)/carbon(i)
            rmir(i) = rmirsprc(i)/carbon(i)
         enddo
         nhrvoc = NHRVSPRC
         do i=1,nhrvoc
            hrvnam(i) = hrvsprc(i)
         enddo
      elseif( idmech .EQ. 7 ) then
         nnames = NSACB7
         call spccb7sa(namecls,numcls,coefcon,coeflux,
     &                 nameyld,numyld,yieldH,yieldL,mwspec)
         do i=1,nnames
            spcnam(i) = spccb7(i)
            carbon(i) = carbcb7(i)
            rkoh(i) = rkohcb7(i)/carbon(i)
            rmir(i) = rmircb7(i)/carbon(i)
         enddo
         nhrvoc = NHRVCB7
         do i=1,nhrvoc
            hrvnam(i) = hrvcb7(i)
         enddo
      elseif( idmech .EQ. 9 ) then        ! RACM2s21
         nnames = NSARC2
         call spcrc2sa(namecls,numcls,coefcon,coeflux,
     &                 nameyld,numyld,yieldH,yieldL,mwspec)
         do i=1,nnames
            spcnam(i) = spcrc2(i)
            carbon(i) = carbrc2(i)
            rkoh(i) = rkohrc2(i)/carbon(i)
            rmir(i) = rmirrc2(i)/carbon(i)
         enddo
         nhrvoc = NHRVRC2
         do i=1,nhrvoc
            hrvnam(i) = hrvrc2(i)
         enddo
      else
         write(iout,'(//,a)') 'ERROR in STABSA:'
         write(iout,'(/,1X,A,I10)')
     &              'Unknown chemical mechanism ID number: ',idmech
         write(iout,'(/,1X,2A)')
     &              'Ozone source apportionment is not available ',
     &              'for this chemical mechanism'
         call camxerr()
      endif
c
c  --- make sure there is enough space ---
c
      ntrcls = 0
      do i=1,MXALCLS
        if( numcls(i) .GT. 0 ) ntrcls = ntrcls + 1
      enddo
c
c  --- call routine to allocate arrays by class ---
c
      call alloc_tracer_class(nspec)
      call zeros(trspmap,nspec*ntrcls)
      call zeros(fluxmap,nspec*ntrcls)
      call zeros(yhratmap,nspec*ntrcls)
      call zeros(ylratmap,nspec*ntrcls)
c
c  --- initialize the pointers into the gridded conc array ---
c
      ntrcls = 0
      iptrbeg = 1
      idxo3n = 0
      idxo3v = 0
      idxoon = 0
      idxoov = 0
      do i=1,MXALCLS
        if( numcls(i) .GT. 0 ) then
           ntrcls = ntrcls + 1
           iptcls(ntrcls) = iptrbeg
           ipttrc(ntrcls) = iptrbeg
           iptrbeg = iptrbeg + ncount + nbdic
           nptcls(ntrcls) = iptrbeg - 1
           npttrc(ntrcls) = iptrbeg - 1
           iemcls(ntrcls) = iptcls(ntrcls) + nbdic
           idxcls(ntrcls) = i
           idxipt(i) = ntrcls
           if( i .EQ. ITRO3N ) idxo3n = ntrcls
           if( i .EQ. ITRO3V ) idxo3v = ntrcls
           if( i .EQ. ITROON ) idxoon = ntrcls
           if( i .EQ. ITROOV ) idxoov = ntrcls
        endif
      enddo
      soap3_ntrcls = ntrcls 
      if( luse_soap3 .AND. lsoa ) soap3_ntrcls = ntrcls + (IDX_SOAP_POA_BB-IDX_SOAP_POA_OP)+1
      ipttim = iptrbeg
      iemtim = ipttim
c
c   --- adjust the pointers for O3N and O3V (together make up O3)
c
      if( idxo3n .GT. 0 .AND. idxo3v .GT. 0 ) then
        ipttrc(idxo3v) = 0
        npttrc(idxo3n) = npttrc(idxo3v)
        npttrc(idxo3v) = 0
      endif
c
c   --- adjust the pointers for OON and OOV (together make up RGN)
c
      if( idxoon .GT. 0 .AND. idxoov .GT. 0 ) then
        ipttrc(idxoov) = 0
        npttrc(idxoon) = npttrc(idxoov)
        npttrc(idxoov) = 0
      endif
c
c   --- calculate the number of tracers at first time step
c       (1 timing release is updated later) ---
c
      nreles = 0
      ntotsp = ipttim-1
      notimespc = ipttim-1
      nsaspc = ipttim-1
c
c  --- calculate the number of timing tracers there will be and put
c      the names into the names array ---
c
      if( ntrtim .GT. 0 ) then
        ibegdt = idate
        btim = begtim/100.
        ienddt = jdate
        etim = endtim/100.
        if( etim .EQ. 0. ) then
            etim = 24.
            ienddt = ienddt - 1
        endif
        timnow = btim
        idtnow = ibegdt
        nhours = (ienddt-ibegdt)*24 + INT( etim - btim )
        do i=1,nhours
           if( MOD( INT(timnow), 24/ntrtim ) .EQ. 0 .OR. i .EQ. 1) then
              do j=1,nregin
                  ntotsp = ntotsp + 2
              enddo
           endif
           timnow = timnow + 1.0
           if( timnow .EQ. 24.0 ) then
               timnow = 0.
               idtnow = idtnow + 1
           endif
        enddo
        if( ntotsp .GT. MXTRSP ) goto 7000
      endif
c
c   --- if doing SOAP3 set the number species that includes the SOAP3
c       emissions species ---
c
      soap3_ntotsp = ntotsp
      soap3_nsaspc = nsaspc
      if( luse_soap3 .AND. lsoa ) then
         soap3_clsnam(IDX_SOAP_POA_OP) = 'POP'
         soap3_clsnam(IDX_SOAP_POA_GV) = 'PGV'
         soap3_clsnam(IDX_SOAP_POA_DV) = 'PDV'
         soap3_clsnam(IDX_SOAP_POA_MC) = 'PMC'
         soap3_clsnam(IDX_SOAP_POA_IC) = 'PIC'
         soap3_clsnam(IDX_SOAP_POA_AV) = 'PAV'
         soap3_clsnam(IDX_SOAP_POA_BB) = 'PBB'
         do i=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
            soap3_iptcls(i) = iptrbeg
            iptrbeg = iptrbeg + ncount
            soap3_nptcls(i) = iptrbeg - 1
            soap3_iemcls(i) = soap3_iptcls(i)
         enddo
         soap3_ntotsp = iptrbeg-1
         soap3_nsaspc = iptrbeg-1
      endif
c
c   --- call routine to allocate the arrays that depend on total
c       tracer species ---
c
      if( .NOT. lddm .AND. .NOT. lhddm ) then
         mvecedge = MAX(ncol(1),nrow(1))
         call alloc_tracer_specs(ngrid,ncol,nrow,nlay,nlayers_ems,
     &                                       tectyp,mvecedge,iout)
         call alloc_tracer_vdep(ngrid,ncol,nrow,ntotsp)
      endif
c
      if( ntrtim .EQ. 0 ) npttim = 0
c
c   --- loop over the species in the chemparam file and 
c       intialize the flags ---
c
      do ispc=1,nspec
         lvocsp(ispc) = .FALSE.
         lnoxsp(ispc) = .FALSE.
         lo3sp(ispc) = .FALSE.
         lvocsoa(ispc) = .FALSE.
         lhrvoc(ispc) = .FALSE.
c
c   --- loop over the tracer classes ---
c
         do i=1,ntrcls
              itr = idxcls(i)
c
c   ---- loop over the model species making up this class ---
c
               do imod=1,numcls(itr)
c
c   --- check the name, if it matches load the global array ---
c
                 if( spname(ispc) .EQ. namecls(itr,imod) ) then
                     trspmap(ispc,i) = coefcon(itr,imod)
                     fluxmap(ispc,i) = coeflux(itr,imod)
                     lusespc(ispc) = .TRUE.
                 endif
               enddo
c
c   ---- loop over the model species for yields ---
c
               do imod=1,numyld(itr)
c
c   --- check the name, if it matches load the global array ---
c
                 if( spname(ispc) .EQ. nameyld(itr,imod) ) then
                     yhratmap(ispc,i) = yieldH(itr,imod)
                     ylratmap(ispc,i) = yieldL(itr,imod)
                 endif
               enddo
c
c   ---- next tracer class ---
c
         enddo
c
c   --- check for O3 species ----
c
         if( spname(ispc) .EQ. o3name ) then
              lo3sp(ispc) = .TRUE.
         endif
c
c   --- check for NOx species ----
c
         if( spname(ispc) .EQ. nonam .OR. 
     &           spname(ispc) .EQ. no2nam .OR. spname(ispc) .EQ. hononam ) then
              lnoxsp(ispc) = .TRUE.
         endif
c
c   ---- check for reactive VOC species ---
c
         do i=1,nnames
            if( spname(ispc) .EQ. spcnam(i) ) then
               lvocsp(ispc) = .TRUE.
               crbnum(ispc) = carbon(i)
               rkohrt(ispc) = rkoh(i)
               rmirrt(ispc) = rmir(i)
            endif
         enddo
c
c   ---- check for SOA-producing VOC species for PSAT ---
c   NOTE: 1. Carbon numbers and MIR values are irrelevant to PSAT
c         2. Only needed for calculating reactivity/yield weighting factors
c
         if( lsoa ) then
           do i=1,NSASOAP
              if( spname(ispc) .EQ. spcsoap(i) ) then
                 lvocsoa(ispc) = .TRUE.
                 rkohrt(ispc) = rkohsoap(i)
              endif
           enddo
         endif
c
c   ---- check for additional VOC species that produce SOA ---
c   NOTE: only for DDM, otherwise they should NOT be added
c
         if ( lddm ) then
            if ( idmech .NE. 5 .OR. idmech .NE. 7 ) then
               if (spname(ispc).EQ.'SQT       ') lvocsp(ispc) = .TRUE.
            endif
            if( spname(ispc) .EQ. 'IVOA      ' ) lvocsp(ispc) = .TRUE.
            if( spname(ispc) .EQ. 'SVOA      ' ) lvocsp(ispc) = .TRUE.
         endif
c
c   ---- check for HRVOC species ---
c   NOTE: set for DDM
c
         do i=1,nhrvoc
            if( spname(ispc) .EQ. hrvnam(i) ) lhrvoc(ispc) = .TRUE.
         enddo
c
c  --- next simulation species ---
c
      enddo
c      
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error Messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in STABSA:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,'(2A)') ' You are using timing tracers but have not ',
     &                       'allocated space for the timing releases.'
      write(iout,*) 'Please change the value for parameter MXTRSP.'
      write(iout,*) 'It should be set to a value of at least: ',ntotsp
      call flush(iout)
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in STABSA:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,'(2A)') ' You are using SOAP3 but have not ',
     &                       'allocated space for the timing releases.'
      write(iout,*) 'Please change the value for parameter MXTRSP.'
      write(iout,*) 'It should be set to a value of at least: ',soap3_ntotsp
      call flush(iout)
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
