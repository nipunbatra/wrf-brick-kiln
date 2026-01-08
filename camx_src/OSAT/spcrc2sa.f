c**** SPCRC2SA
c
      subroutine spcrc2sa(namecls,numcls,coefcon,coeflux,
     &                     nameyld,numyld,yieldH,yieldL,molwt)
      use chmstry
      use tracer
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine initializes the species variables that will be used
c   to determine the species that will be included in each tracer
c   species class
c    Argument descriptions:
c     Input:  
c     Output:  
c       namecls  C  array of regular model species in each tracer class
c       numcls   I  the number of species contributing to this class
c       coefcon  R  2-D array of coefficients for making linear combo
c                   this one for concentrations and emissions
c       coeflux  R  2-D array of coefficients for making linear combo
c                   this one for fluxes
c       nameyld  C  array of regular model species for the yields of
c                   each tracer class
c       numyld   I  number of regular model species for the yields of
c                   each tracer class
c       yieldH   R  high-NOx yield rates for each species in each class
c       yieldL   R  low-NOx yield rates for each species in each class
c       molwt    R  molecular weight of model species
c
c     Copyright 1996 - 2025
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     09/02/21   --gyarwood-    Created RACM2s21 version (RC2 for short)
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'soap.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*10 namecls(MXALCLS,MXSPEC)
      integer      numcls(MXALCLS)
      real         coefcon(MXALCLS,MXSPEC)
      real         coeflux(MXALCLS,MXSPEC)
      character*10 nameyld(MXALCLS,MXSPEC)
      integer      numyld(MXALCLS)
      real         yieldH(MXALCLS,MXSPEC),yieldL(MXALCLS,MXSPEC)
      real         molwt(MXSPEC)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   NUMRC2    I   number of RACM2 species used by PSAT treatment
c
      integer NUMRC2
c
      parameter( NUMRC2 = 66 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 rc2nam(NUMRC2)
      integer*4    i, j
      real         wtrc2(NUMRC2)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c

      data rc2nam /
     &        'NO        ', 'NO2       ', 'HCHO      ', 'ACD       ',
     &        'ALD       ', 'ACT       ', 'MEK       ', 'KET       ',
     &        'MOH       ', 'EOH       ', 'ROH       ', 'ETEG      ',
     &        'ISO       ', 'API       ', 'LIM       ', 'BEN       ',
     &        'TOL       ', 'XYM       ', 'XYP       ', 'XYO       ',
     &        'ETH       ', 'HC3       ', 'HC5       ', 'HC8       ',
     &        'ETE       ', 'OLT       ', 'OLI       ', 'ACE       ',
     &        'SO2       ', 'HONO      ', 'NO3       ', 'N2O5      ',
     &        'PAN       ', 'PPN       ', 'MPAN      ', 'HNO4      ',
     &        'ONIT      ', 'ISON      ', 'NALD      ', 'DIEN      ',
     &        'HNO3      ', 'NH3       ', 'IVOA      ', 'SVOA      ',
     &        'ACG2      ', 'ACG1      ', 'BCG2      ', 'BCG1      ',
     &        'HG0       ', 'HG2       ', 'PSO4      ', 'PNO3      ',
     &        'PNH4      ', 'PEC       ', 'HOA       ', 'FCRS      ',
     &        'FPRM      ', 'CCRS      ', 'CPRM      ', 'AOA2      ',
     &        'AOA1      ', 'BOA2      ', 'BOA1      ', 'AOA0      ',
     &        'BOA0      ', 'HGP       '/ 

      data wtrc2 /
     &         46.0       ,  46.0       ,  30.0       ,  44.0       ,
     &         58.1       ,  58.1       ,  72.1       ,  86.1       ,
     &         32.0       ,  46.1       ,  60.1       ,  62.0       ,
     &         68.1       , 136.2       , 136.2       ,  78.1       ,
     &         92.1       , 119.9       , 121.2       , 119.9       ,
     &         30.1       ,  52.5       ,  80.5       , 112.8       ,
     &         28.0       ,  53.3       ,  70.1       ,  26.0       ,
     &         64.0       ,  46.0       ,  62.0       , 108.0       ,
     &        121.0       , 135.0       , 148.0       ,  79.0       ,
     &        119.1       , 147.1       , 105.0       ,  56.1       ,
     &         63.0       ,  17.0       , 212.0       , 212.0       ,
     &        150.0       , 150.0       , 180.0       , 180.0       ,
     &        200.59      , 200.59      ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- initialize all arrays ---
c
      do i=1,MXALCLS
        numcls(i) = 0
        numyld(i) = 0
        do j=1,MXSPEC
          namecls(i,j) = ' '
          nameyld(i,j) = ' '
          coefcon(i,j) = 0.0
          coeflux(i,j) = 0.0
          yieldH(i,j) = 0.0
          yieldL(i,j) = 0.0
        enddo
      enddo
c
c  --- set the molecular weight for each modeled species ---
c
      do i=1,nspec
        molwt(i) = 1.0
        do j=1,NUMRC2
           if( spname(i) .EQ. rc2nam(j) ) molwt(i) = wtrc2(j)
        enddo
      enddo
c
c   --- if doing the OZONE or NITRATE species ---
c
      if( lozone .OR. lnitrate ) then
c
c  ---- NIT species ---
c
         numcls(ITRNIT) = 2
         namecls(ITRNIT,1) = 'NO'
         coefcon(ITRNIT,1) = 1.0
         coeflux(ITRNIT,1) = 1.0
         namecls(ITRNIT,2) = 'HONO'
         coefcon(ITRNIT,2) = 1.0
         coeflux(ITRNIT,2) = 1.0
c
c  ---- RGN species ---
c
         numcls(ITRRGN) = 3
         namecls(ITRRGN,1) = 'NO2'
         coefcon(ITRRGN,1) = 1.0
         coeflux(ITRRGN,1) = 1.0
         namecls(ITRRGN,2) = 'NO3'
         coefcon(ITRRGN,2) = 1.0
         coeflux(ITRRGN,2) = 1.0
         namecls(ITRRGN,3) = 'N2O5'
         coefcon(ITRRGN,3) = 2.0
         coeflux(ITRRGN,3) = 2.0
c
c  ---- TPN species ---
c
         numcls(ITRTPN) = 4
         namecls(ITRTPN,1) = 'PAN'
         coefcon(ITRTPN,1) = 1.0
         coeflux(ITRTPN,1) = 1.0
         namecls(ITRTPN,2) = 'PPN'
         coefcon(ITRTPN,2) = 1.0
         coeflux(ITRTPN,2) = 1.0
         namecls(ITRTPN,3) = 'MPAN'
         coefcon(ITRTPN,3) = 1.0
         coeflux(ITRTPN,3) = 1.0
         namecls(ITRTPN,4) = 'HNO4'
         coefcon(ITRTPN,4) = 1.0
         coeflux(ITRTPN,4) = 1.0
c
c  ---- NTR species ---
c
         numcls(ITRNTR) = 3
         namecls(ITRNTR,1) = 'ONIT'
         coefcon(ITRNTR,1) = 1.0
         coeflux(ITRNTR,1) = 1.0
         namecls(ITRNTR,2) = 'ISON'
         coefcon(ITRNTR,2) = 1.0
         coeflux(ITRNTR,2) = 1.0
         namecls(ITRNTR,3) = 'NALD'
         coefcon(ITRNTR,3) = 1.0
         coeflux(ITRNTR,3) = 1.0
c
c  ---- HN3 species ---
c
         numcls(ITRHN3) = 1
         namecls(ITRHN3,1) = 'HNO3'
         coefcon(ITRHN3,1) = 1.0
         coeflux(ITRHN3,1) = 1.0
      endif
c
c   --- if doing the OZONE species ---
c
      if( lozone ) then
c
c  ---- VOC species ---
c
         numcls(ITRVOC) = 26
         namecls(ITRVOC,1) = 'HCHO'
         coefcon(ITRVOC,1) = 1.0
         coeflux(ITRVOC,1) = 1.0
         namecls(ITRVOC,2) = 'ACD'
         coefcon(ITRVOC,2) = 2.0
         coeflux(ITRVOC,2) = 2.0
         namecls(ITRVOC,3) = 'ALD'
         coefcon(ITRVOC,3) = 3.0
         coeflux(ITRVOC,3) = 3.0
         namecls(ITRVOC,4) = 'ACT'
         coefcon(ITRVOC,4) = 3.0
         coeflux(ITRVOC,4) = 3.0
         namecls(ITRVOC,5) = 'MEK'
         coefcon(ITRVOC,5) = 4.0
         coeflux(ITRVOC,5) = 4.0
         namecls(ITRVOC,6) = 'KET'
         coefcon(ITRVOC,6) = 5.0
         coeflux(ITRVOC,6) = 5.0
         namecls(ITRVOC,7) = 'ETE'
         coefcon(ITRVOC,7) = 2.0
         coeflux(ITRVOC,7) = 2.0
         namecls(ITRVOC,8) = 'OLT'
         coefcon(ITRVOC,8) = 3.8
         coeflux(ITRVOC,8) = 3.8
         namecls(ITRVOC,9) = 'OLI'
         coefcon(ITRVOC,9) = 5.0
         coeflux(ITRVOC,9) = 5.0
         namecls(ITRVOC,10) = 'ISO'
         coefcon(ITRVOC,10) = 5.0
         coeflux(ITRVOC,10) = 5.0
         namecls(ITRVOC,11) = 'API'
         coefcon(ITRVOC,11) = 10.0
         coeflux(ITRVOC,11) = 10.0
         namecls(ITRVOC,12) = 'LIM'
         coefcon(ITRVOC,12) = 10.0
         coeflux(ITRVOC,12) = 10.0
         namecls(ITRVOC,13) = 'BEN'
         coefcon(ITRVOC,13) = 6.0
         coeflux(ITRVOC,13) = 6.0
         namecls(ITRVOC,14) = 'TOL'
         coefcon(ITRVOC,14) = 7.1
         coeflux(ITRVOC,14) = 7.1
         namecls(ITRVOC,15) = 'XYM'
         coefcon(ITRVOC,15) = 8.9
         coeflux(ITRVOC,15) = 8.9
         namecls(ITRVOC,16) = 'XYP'
         coefcon(ITRVOC,16) = 9.0
         coeflux(ITRVOC,16) = 9.0
         namecls(ITRVOC,17) = 'XYO'
         coefcon(ITRVOC,17) = 8.9
         coeflux(ITRVOC,17) = 8.9
         namecls(ITRVOC,18) = 'MOH'
         coefcon(ITRVOC,18) = 1.0
         coeflux(ITRVOC,18) = 1.0
         namecls(ITRVOC,19) = 'EOH'
         coefcon(ITRVOC,19) = 2.0
         coeflux(ITRVOC,19) = 2.0
         namecls(ITRVOC,20) = 'ROH'
         coefcon(ITRVOC,20) = 3.0
         coeflux(ITRVOC,20) = 3.0
         namecls(ITRVOC,21) = 'ETH'
         coefcon(ITRVOC,21) = 2.0
         coeflux(ITRVOC,21) = 2.0
         namecls(ITRVOC,22) = 'HC3'
         coefcon(ITRVOC,22) = 3.6
         coeflux(ITRVOC,22) = 3.6
         namecls(ITRVOC,23) = 'HC5'
         coefcon(ITRVOC,23) = 5.6
         coeflux(ITRVOC,23) = 5.6
         namecls(ITRVOC,24) = 'HC8'
         coefcon(ITRVOC,24) = 7.9
         coeflux(ITRVOC,24) = 7.9
         namecls(ITRVOC,25) = 'ACE'
         coefcon(ITRVOC,25) = 2.0
         coeflux(ITRVOC,25) = 2.0
         namecls(ITRVOC,26) = 'ETEG'
         coefcon(ITRVOC,26) = 2.0
         coeflux(ITRVOC,26) = 2.0
         namecls(ITRVOC,27) = 'DIEN'
         coefcon(ITRVOC,27) = 4.0
         coeflux(ITRVOC,27) = 4.0
c
c  ---- O3-NOx species ---
c
         numcls(ITRO3N) = 1
         namecls(ITRO3N,1) = 'O3'
         coefcon(ITRO3N,1) = 0.5
         coeflux(ITRO3N,1) = 1.0
c
c  ---- O3-VOC species ---
c
         numcls(ITRO3V) = 1
         namecls(ITRO3V,1) = 'O3'
         coefcon(ITRO3V,1) = 0.5
         coeflux(ITRO3V,1) = 1.0
c
c  ---- Odd-oxygen in RGN from O3N ---
c
         numcls(ITROON) = 3
         namecls(ITROON,1) = 'NO2'
         coefcon(ITROON,1) = 0.5
         coeflux(ITROON,1) = 1.0
         namecls(ITROON,2) = 'NO3'
         coefcon(ITROON,2) = 0.5
         coeflux(ITROON,2) = 1.0
         namecls(ITROON,3) = 'N2O5'
         coefcon(ITROON,3) = 1.0
         coeflux(ITROON,3) = 2.0
c
c  ---- Odd-oxygen in RGN from O3V ---
c
         numcls(ITROOV) = 3
         namecls(ITROOV,1) = 'NO2'
         coefcon(ITROOV,1) = 0.5
         coeflux(ITROOV,1) = 1.0
         namecls(ITROOV,2) = 'NO3'
         coefcon(ITROOV,2) = 0.5
         coeflux(ITROOV,2) = 1.0
         namecls(ITROOV,3) = 'N2O5'
         coefcon(ITROOV,3) = 1.0
         coeflux(ITROOV,3) = 2.0
      endif
c
c   --- if doing the SULFATE species ---
c
      if( lsulfate ) then
c
c  ---- SO2 species ---
c
        numcls(ITRSO2) = 1
        namecls(ITRSO2,1) = 'SO2'
        coefcon(ITRSO2,1) = 1.0
        coeflux(ITRSO2,1) = 1.0
c
c  ---- PS4 species ---
c
        numcls(ITRPS4) = 1
        namecls(ITRPS4,1) = 'PSO4'
        coefcon(ITRPS4,1) = 1.0
        coeflux(ITRPS4,1) = 1.0
      endif
c
c   --- if doing the NITRATE species ---
c
      if( lnitrate ) then
c
c  ---- NIT species ---
c
c         defined above
c
c  ---- RGN species ---
c
c         defined above
c
c  ---- TPN species ---
c
c         defined above
c
c  ---- NTR species ---
c
c         defined above
c
c  ---- HN3 species ---
c
c         defined above
c
c  ---- PN3 species ---
c
          numcls(ITRPN3) = 1
          namecls(ITRPN3,1) = 'PNO3'
          coefcon(ITRPN3,1) = 1.0
          coeflux(ITRPN3,1) = 1.0
c
c  ---- NH3 species ---
c
          numcls(ITRNH3) = 1
          namecls(ITRNH3,1) = 'NH3'
          coefcon(ITRNH3,1) = 1.0
          coeflux(ITRNH3,1) = 1.0
c
c  ---- PN4 species ---
c
          numcls(ITRPN4) = 1
          namecls(ITRPN4,1) = 'PNH4'
          coefcon(ITRPN4,1) = 1.0
          coeflux(ITRPN4,1) = 1.0
      endif
c
c   --- if doing the SOA species ---
c
      if( lsoa ) then
c
c  ---- ARO species ---
c
          numcls(ITRARO) = 5
          namecls(ITRARO,1) = 'BEN'
          coefcon(ITRARO,1) = 1.0
          coeflux(ITRARO,1) = 1.0
          namecls(ITRARO,2) = 'TOL'
          coefcon(ITRARO,2) = 1.0
          coeflux(ITRARO,2) = 1.0
          namecls(ITRARO,3) = 'XYM'
          coefcon(ITRARO,3) = 1.0
          coeflux(ITRARO,3) = 1.0
          namecls(ITRARO,4) = 'XYP'
          coefcon(ITRARO,4) = 1.0
          coeflux(ITRARO,4) = 1.0
          namecls(ITRARO,5) = 'XYO'
          coefcon(ITRARO,5) = 1.0
          coeflux(ITRARO,5) = 1.0
c
c  ---- ISP species ---
c
          numcls(ITRISP) = 1
          namecls(ITRISP,1) = 'ISO '
          coefcon(ITRISP,1) = 1.0
          coeflux(ITRISP,1) = 1.0
c
c  ---- TRP species ---
c
          numcls(ITRTRP) = 1
          namecls(ITRTRP,1) = 'LIM'
          coefcon(ITRTRP,1) = 1.0
          coeflux(ITRTRP,1) = 1.0
c
c  ---- APN species ---
c
          numcls(ITRAPN) = 1
          namecls(ITRAPN,1) = 'API'
          coefcon(ITRAPN,1) = 1.0
          coeflux(ITRAPN,1) = 1.0
c
c  ---- SQT species ---
c
          numcls(ITRSQT) = 1
          namecls(ITRSQT,1) = 'SQT'
          coefcon(ITRSQT,1) = 1.0
          coeflux(ITRSQT,1) = 1.0
c
c  ---- AG2 species ---
c
          numcls(ITRAG2) = 1
          namecls(ITRAG2,1) = 'ACG2'
          coefcon(ITRAG2,1) = 1.0
          coeflux(ITRAG2,1) = 1.0
          numyld(ITRAG2) = 7
          nameyld(ITRAG2,1) = 'BEN'
          yieldH(ITRAG2,1) = MAX(EPSYLD, y_h(1,1))
          yieldL(ITRAG2,1) = MAX(EPSYLD, y_l(1,1))
          nameyld(ITRAG2,2) = 'TOL'
          yieldH(ITRAG2,2) = MAX(EPSYLD, y_h(1,2))
          yieldL(ITRAG2,2) = MAX(EPSYLD, y_l(1,2))
          nameyld(ITRAG2,3) = 'XYM'
          yieldH(ITRAG2,3) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRAG2,3) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRAG2,4) = 'XYP'
          yieldH(ITRAG2,4) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRAG2,4) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRAG2,5) = 'XYO'
          yieldH(ITRAG2,5) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRAG2,5) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRAG2,6) = 'IVOA'
          yieldH(ITRAG2,6) = MAX(EPSYLD, y_h(1,4))
          yieldL(ITRAG2,6) = MAX(EPSYLD, y_l(1,4))
          nameyld(ITRAG2,7) = 'SVOA'
          yieldH(ITRAG2,7) = MAX(EPSYLD, y_h(1,5))
          yieldL(ITRAG2,7) = MAX(EPSYLD, y_l(1,5))
c
c  ---- AG1 species ---
c
          numcls(ITRAG1) = 1
          namecls(ITRAG1,1) = 'ACG1'
          coefcon(ITRAG1,1) = 1.0
          coeflux(ITRAG1,1) = 1.0
          numyld(ITRAG1) = 7
          nameyld(ITRAG1,1) = 'BEN'
          yieldH(ITRAG1,1) = MAX(EPSYLD, y_h(2,1))
          yieldL(ITRAG1,1) = MAX(EPSYLD, y_l(2,1))
          nameyld(ITRAG1,2) = 'TOL'
          yieldH(ITRAG1,2) = MAX(EPSYLD, y_h(2,2))
          yieldL(ITRAG1,2) = MAX(EPSYLD, y_l(2,2))
          nameyld(ITRAG1,3) = 'XYM'
          yieldH(ITRAG1,3) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRAG1,3) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRAG1,4) = 'XYP'
          yieldH(ITRAG1,4) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRAG1,4) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRAG1,5) = 'XYO'
          yieldH(ITRAG1,5) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRAG1,5) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRAG1,6) = 'IVOA'
          yieldH(ITRAG1,6) = MAX(EPSYLD, y_h(2,4))
          yieldL(ITRAG1,6) = MAX(EPSYLD, y_l(2,4))
          nameyld(ITRAG1,7) = 'SVOA'
          yieldH(ITRAG1,7) = MAX(EPSYLD, y_h(2,5))
          yieldL(ITRAG1,7) = MAX(EPSYLD, y_l(2,5))
c
c  ---- BG2 species ---
c
          numcls(ITRBG2) = 1
          namecls(ITRBG2,1) = 'BCG2'
          coefcon(ITRBG2,1) = 1.0
          coeflux(ITRBG2,1) = 1.0
c
c  ---- BG1 species ---
c
          numcls(ITRBG1) = 1
          namecls(ITRBG1,1) = 'BCG1'
          coefcon(ITRBG1,1) = 1.0
          coeflux(ITRBG1,1) = 1.0
c
c  ---- PA2 species ---
c
          numcls(ITRPA2) = 1
          namecls(ITRPA2,1) = 'AOA2'
          coefcon(ITRPA2,1) = 1.0
          coeflux(ITRPA2,1) = 1.0
c
c  ---- PA1 species ---
c
          numcls(ITRPA1) = 1
          namecls(ITRPA1,1) = 'AOA1'
          coefcon(ITRPA1,1) = 1.0
          coeflux(ITRPA1,1) = 1.0
c
c  ---- PB2 species ---
c
          numcls(ITRPB2) = 1
          namecls(ITRPB2,1) = 'BOA2'
          coefcon(ITRPB2,1) = 1.0
          coeflux(ITRPB2,1) = 1.0
c
c  ---- PB1 species ---
c
          numcls(ITRPB1) = 1
          namecls(ITRPB1,1) = 'BOA1'
          coefcon(ITRPB1,1) = 1.0
          coeflux(ITRPB1,1) = 1.0
c
c  ---- PA0 species ---
c
          numcls(ITRPA0) = 1
          namecls(ITRPA0,1) = 'AOA0'
          coefcon(ITRPA0,1) = 1.0
          coeflux(ITRPA0,1) = 1.0
          ! SPECIAL CASE: non-volatile CG; anthro -> AOA0 directly (skip CG)
          numyld(ITRPA0) = 7
          nameyld(ITRPA0,1) = 'BEN'
          yieldH(ITRPA0,1) = MAX(EPSYLD, y_h(3,1))
          yieldL(ITRPA0,1) = MAX(EPSYLD, y_l(3,1))
          nameyld(ITRPA0,2) = 'TOL'
          yieldH(ITRPA0,2) = MAX(EPSYLD, y_h(3,2))
          yieldL(ITRPA0,2) = MAX(EPSYLD, y_l(3,2))
          nameyld(ITRPA0,3) = 'XYM'
          yieldH(ITRPA0,3) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPA0,3) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPA0,4) = 'XYP'
          yieldH(ITRPA0,4) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPA0,4) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPA0,5) = 'XYO'
          yieldH(ITRPA0,5) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPA0,5) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPA0,6) = 'IVOA'
          yieldH(ITRPA0,6) = MAX(EPSYLD, y_h(3,4))
          yieldL(ITRPA0,6) = MAX(EPSYLD, y_l(3,4))
          nameyld(ITRPA0,7) = 'SVOA'
          yieldH(ITRPA0,7) = MAX(EPSYLD, y_h(3,5))
          yieldL(ITRPA0,7) = MAX(EPSYLD, y_l(3,5))
c
c  ---- PB0 species ---
c
          numcls(ITRPB0) = 1
          namecls(ITRPB0,1) = 'BOA0'
          coefcon(ITRPB0,1) = 1.0
          coeflux(ITRPB0,1) = 1.0
c
c  ---- IVA species ---
c
          numcls(ITRIVA) = 1
          namecls(ITRIVA,1) = 'IVOA'
          coefcon(ITRIVA,1) = 1.0
          coeflux(ITRIVA,1) = 1.0
c
c  ---- SVA species ---
c
          numcls(ITRSVA) = 1
          namecls(ITRSVA,1) = 'SVOA'
          coefcon(ITRSVA,1) = 1.0
          coeflux(ITRSVA,1) = 1.0
c
c  ---- HOA species ---
c
          numcls(ITRHOA) = 1
          namecls(ITRHOA,1) = 'HOA'
          coefcon(ITRHOA,1) = 1.0
          coeflux(ITRHOA,1) = 1.0
      endif
c
c   --- if doing the PRIMARY species ---
c
      if( lprimary ) then
c
c  ---- PEC species ---
c
          numcls(ITRPEC) = 1
          namecls(ITRPEC,1) = 'PEC'
          coefcon(ITRPEC,1) = 1.0
          coeflux(ITRPEC,1) = 1.0
c
c  ---- PFC species ---
c
          numcls(ITRPFC) = 1
          namecls(ITRPFC,1) = 'FCRS'
          coefcon(ITRPFC,1) = 1.0
          coeflux(ITRPFC,1) = 1.0
c
c  ---- PFN species ---
c
          numcls(ITRPFN) = 1
          namecls(ITRPFN,1) = 'FPRM'
          coeflux(ITRPFN,1) = 1.0
          coefcon(ITRPFN,1) = 1.0
c
c  ---- PCC species ---
c
          numcls(ITRPCC) = 1
          namecls(ITRPCC,1) = 'CCRS'
          coefcon(ITRPCC,1) = 1.0
          coeflux(ITRPCC,1) = 1.0
c
c  ---- PCS species ---
c
          numcls(ITRPCS) = 1
          namecls(ITRPCS,1) = 'CPRM'
          coefcon(ITRPCS,1) = 1.0
          coeflux(ITRPCS,1) = 1.0

          if( lmineral ) then
c
c  ---- PFE species ---
c
              numcls(ITRPFE) = 1
              namecls(ITRPFE,1) = 'PFE'
              coefcon(ITRPFE,1) = 1.0
              coeflux(ITRPFE,1) = 1.0
c
c  ---- PMN species ---
c
              numcls(ITRPMN) = 1
              namecls(ITRPMN,1) = 'PMN'
              coefcon(ITRPMN,1) = 1.0
              coeflux(ITRPMN,1) = 1.0
c
c  ---- PMG species ---
c
              numcls(ITRPMG) = 1
              namecls(ITRPMG,1) = 'PMG'
              coefcon(ITRPMG,1) = 1.0
              coeflux(ITRPMG,1) = 1.0
c
c  ---- PK species ---
c
              numcls(ITRPK) = 1
              namecls(ITRPK,1) = 'PK'
              coefcon(ITRPK,1) = 1.0
              coeflux(ITRPK,1) = 1.0
c
c  ---- PCA species ---
c
              numcls(ITRPCA) = 1
              namecls(ITRPCA,1) = 'PCA'
              coefcon(ITRPCA,1) = 1.0
              coeflux(ITRPCA,1) = 1.0
c
c  ---- PAL species ---
c
              numcls(ITRPAL) = 1
              namecls(ITRPAL,1) = 'PAL'
              coefcon(ITRPAL,1) = 1.0
              coeflux(ITRPAL,1) = 1.0
c
c  ---- PSI species ---
c
              numcls(ITRPSI) = 1
              namecls(ITRPSI,1) = 'PSI'
              coefcon(ITRPSI,1) = 1.0
              coeflux(ITRPSI,1) = 1.0
c
c  ---- PTI species ---
c
              numcls(ITRPTI) = 1
              namecls(ITRPTI,1) = 'PTI'
              coefcon(ITRPTI,1) = 1.0
              coeflux(ITRPTI,1) = 1.0
          endif
      endif
c
c   --- if doing the MERCURY species ---
c
      if( lmercury ) then
c       
c  ---- HG0 species ---
c
          numcls(ITRHG0) = 1
          namecls(ITRHG0,1) = 'HG0'
          coefcon(ITRHG0,1) = 1.0
          coeflux(ITRHG0,1) = 1.0
c       
c  ---- HG2 species ---
c
          numcls(ITRHG2) = 1
          namecls(ITRHG2,1) = 'HG2'
          coefcon(ITRHG2,1) = 1.0
          coeflux(ITRHG2,1) = 1.0
c       
c  ---- PHG species ---
c
          numcls(ITRPHG) = 1
          namecls(ITRPHG,1) = 'HGP'
          coefcon(ITRPHG,1) = 1.0
          coeflux(ITRPHG,1) = 1.0
      endif
c       
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
