c----CAMx v7.32 250801
c
c     SOAP.INC contains parameters for OA modules
c
c      Copyright 1996 - 2025
c     Ramboll
c
c     Modifications:
c        08/25/16 -bkoo-   Revised to include parameters for multiple OA schemes
c        11/28/16 -bkoo-   Revised number of VBS precursors
c        10/16/17 -bkoo-   Added sfacJSOA; removed LSOAP2
c        04/22/19 -gy---   SOAP2.2: Update TRP SOA yields for RO2 autoxidation
c                                   JSOA reduced to 0.001 (Malecha et al. 2018)
c                                   kpolya deactivated because duplicative
c        04/05/24 -cemery- SOAP3: rename species, POA->HOA (unevaporated POA)
c                                 add SVOA (evaporated POA)
c                                 updated yields
c        04/05/24 -gy---   SOAP3: seperate NO3 reaction SOA yields 
c
c-----------------------------------------------------------------------
c     Parameters for SOAP
c-----------------------------------------------------------------------
c
      integer, parameter :: NPREC_S  = 9
      integer, parameter :: NSOAP    = 4   ! No. of CG/SOA pairs
      integer, parameter :: NSOAA    = 2   ! No. of anthropogenic CG/SOA pairs
      integer, parameter :: NNVOL    = 2   ! No. of non-volatile products
      integer, parameter :: NSOAPRXN = 15  ! No. of gas-phase reactions in SOACHEM

      real, parameter :: mwhoa  = 220.    ! MW of HOA
      real, parameter :: mwaoa0 = 220.    ! MW of AOA0
      real, parameter :: mwboa0 = 220.    ! MW of BOA0

      real, parameter :: mwsoap(NSOAP) =              ! MW of CG/SOA pairs (g/mol)
     &                                     (/ 150.,   ! ACG2
     &                                        150.,   ! ACG1
     &                                        180.,   ! BCG2
     &                                        180. /) ! BCG1

c
c     SOAP3 COMPLX parameters
c
      real, parameter :: yh_3complx(3,NPREC_S) = reshape(   ! high-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.187,   0.0,     0.004,          ! BNZA
     &              0.092,   0.012,   0.0,            ! TOLA
     &              0.064,   0.0,     0.005,          ! XYLA
     &              0.0,     0.0,     1.00,           ! IVOA
     &              1.10,    0.73,    0.18,           ! SVOA
     &              0.038,   0.008,   0.0,            ! ISP
     &              0.166,   0.038,   0.031,          ! TRP
     &              0.166,   0.038,   0.031,          ! APN
     &              1.417,   0.068,   0.0 /),         ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yl_3complx(3,NPREC_S) = reshape(   ! low-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.131,          ! BNZA
     &              0.0,     0.0,     0.126,          ! TOLA
     &              0.0,     0.0,     0.174,          ! XYLA
     &              0.0,     0.0,     1.00,           ! IVOA
     &              1.10,    0.73,    0.18,           ! SVOA
     &              0.038,   0.008,   0.0,            ! ISP
     &              0.166,   0.038,   0.031,          ! TRP
     &              0.166,   0.038,   0.031,          ! APN
     &              1.417,   0.068,   0.0 /),         ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yn_3complx(3,NPREC_S) = reshape(   ! NO3 reaction yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.00,           ! BNZA (not used)
     &              0.0,     0.0,     0.00,           ! TOLA (not used)
     &              0.0,     0.0,     0.00,           ! XYLA (not used)
     &              0.0,     0.0,     0.00,           ! IVOA (not used)
     &              0.0,     0.0,     0.00,           ! SVOA (not used)
     &              0.038,   0.008,   0.0,            ! ISP
     &              0.325,   0.083,   0.0,            ! TRP
     &              0.159,   0.026,   0.0,            ! APN
     &              1.417,   0.068,   0.0 /),         ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: csatS_3complx(NSOAP) =             ! Sat. conc. of CG/SOA pairs (ug/m3)
     &                                    (/ 14.,     ! ACG2
     &                                        0.31,   ! ACG1
     &                                       26.,     ! BCG2
     &                                        0.45 /) ! BCG1

      real, parameter :: cstempS_3complx(NSOAP) =           ! Ref. temp. for sat. conc. (K)
     &                                     (/ 300.,   ! ACG2
     &                                        300.,   ! ACG1
     &                                        300.,   ! BCG2
     &                                        300. /) ! BCG1

      real, parameter :: deltahS_3complx(NSOAP) =           ! Heat of vap. for CG/SOA pairs (J/mol)
     &                                  (/ 18000.,    ! ACG2
     &                                     18000.,    ! ACG1
     &                                     65000.,    ! BCG2
     &                                     85000.  /) ! BCG1
      real, parameter :: kpolya_3complx   = 0.0   ! polymerization rate for AOA0 [1/hr]; deactivated
      real, parameter :: kpolyb_3complx   = 0.0   ! polymerization rate for BOA0 [1/hr]; deactivated
      real, parameter :: sfacJSOA_3complx = 0.001 ! scaling factor for JSOA: JSOA = sfacJSOA * JNO2

c
c     SOAP3 SIMPLE parameters (non-volatile products only) 
c
      real, parameter :: yh_3simple(3,NPREC_S) = reshape(   ! high-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.10,           ! BNZA
     &              0.0,     0.0,     0.06,           ! TOLA
     &              0.0,     0.0,     0.04,           ! XYLA
     &              0.0,     0.0,     1.00,           ! IVOA
     &              0.0,     0.0,     0.98,           ! SVOA
     &              0.0,     0.0,     0.01,           ! ISP
     &              0.0,     0.0,     0.139,          ! TRP
     &              0.0,     0.0,     0.052,          ! APN
     &              0.0,     0.0,     0.32 /),        ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yl_3simple(3,NPREC_S) = reshape(   ! low-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.20,           ! BNZA
     &              0.0,     0.0,     0.12,           ! TOLA
     &              0.0,     0.0,     0.09,           ! XYLA
     &              0.0,     0.0,     1.00,           ! IVOA
     &              0.0,     0.0,     0.98,           ! SVOA
     &              0.0,     0.0,     0.01,           ! ISP
     &              0.0,     0.0,     0.139,          ! TRP
     &              0.0,     0.0,     0.052,          ! APN
     &              0.0,     0.0,     0.40 /),        ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yn_3simple(3,NPREC_S) = reshape(   ! NO3 reaction yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.00,           ! BNZA (not used)
     &              0.0,     0.0,     0.00,           ! TOLA (not used)
     &              0.0,     0.0,     0.00,           ! XYLA (not used)
     &              0.0,     0.0,     0.00,           ! IVOA (not used)
     &              0.0,     0.0,     0.00,           ! SVOA (not used)
     &              0.0,     0.0,     0.01,           ! ISP
     &              0.0,     0.0,     0.083,          ! TRP
     &              0.0,     0.0,     0.026,          ! APN
     &              0.0,     0.0,     0.40 /),        ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: csatS_3simple(NSOAP) =             ! Sat. conc. of CG/SOA pairs (ug/m3)
     &                                    (/ 14.,     ! ACG2
     &                                        0.31,   ! ACG1
     &                                       26.,     ! BCG2
     &                                        0.45 /) ! BCG1

      real, parameter :: cstempS_3simple(NSOAP) =           ! Ref. temp. for sat. conc. (K)
     &                                     (/ 300.,   ! ACG2
     &                                        300.,   ! ACG1
     &                                        300.,   ! BCG2
     &                                        300. /) ! BCG1

      real, parameter :: deltahS_3simple(NSOAP) =           ! Heat of vap. for CG/SOA pairs (J/mol)
     &                                  (/ 18000.,    ! ACG2
     &                                     18000.,    ! ACG1
     &                                     65000.,    ! BCG2
     &                                     85000.  /) ! BCG1
      real, parameter :: kpolya_3simple   = 0.0   ! polymerization rate for AOA0 [1/hr]; deactivated
      real, parameter :: kpolyb_3simple   = 0.0   ! polymerization rate for BOA0 [1/hr]; deactivated
      real, parameter :: sfacJSOA_3simple = 0.001 ! scaling factor for JSOA: JSOA = sfacJSOA * JNO2
c
c     SOAP3 1.5VDBS parameters: not implemented May 2024,review/update before using
c
      real, parameter :: yh_3dbs(3,NPREC_S) = reshape(   ! high-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.129,   0.204,   0.0,            ! BNZA
     &              0.187,   0.180,   0.0185,         ! TOLA
     &              0.0593,  0.0343,  0.0123,         ! XYLA
     &              1.38,    0.32,    0.0,            ! IVOA
     &              0.48,    0.62,    0.37,           ! SVOA
     &              0.0286,  0.00850, 0.0,            ! ISP
     &              0.0567,  0.034,   0.0433,         ! TRP
     &              0.0567,  0.034,   0.0433,         ! APN
     &              0.105,   0.372,   0.162 /),       ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yl_3dbs(3,NPREC_S) = reshape(   ! low-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.253,   0.0870,  0.0,            ! BNZA
     &              0.406,   0.211,   0.109,          ! TOLA
     &              0.206,   0.216,   0.142,          ! XYLA
     &              1.38,    0.32,    0.0,            ! IVOA
     &              0.48,    0.62,    0.37,           ! SVOA
     &              0.0588,  0.0108,  0.00335,        ! ISP
     &              0.1133,  0.068,   0.0433,         ! TRP
     &              0.1133,  0.068,   0.0433,         ! APN
     &              0.154,   0.453,   0.250 /),       ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yn_3dbs(3,NPREC_S) = reshape(   ! NO3 reaction yield [ppm/ppm]
c     placeholder: use yl_3dbs except APIN is 37% of TRP folowing SIMPLE
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.00,           ! BNZA (not used)
     &              0.0,     0.0,     0.00,           ! TOLA (not used)
     &              0.0,     0.0,     0.00,           ! XYLA (not used)
     &              0.0,     0.0,     0.00,           ! IVOA (not used)
     &              0.0,     0.0,     0.00,           ! SVOA (not used)
     &              0.0588,  0.0108,  0.00335,        ! ISP
     &              0.1133,  0.068,   0.0433,         ! TRP
     &              0.0421,  0.0253,  0.0161,         ! APN
     &              0.154,   0.453,   0.250 /),       ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: csatS_3dbs(NSOAP) =             ! Sat. conc. of CG/SOA pairs (ug/m3)
     &                                    (/ 14.,     ! ACG2
     &                                        0.31,   ! ACG1
     &                                       26.,     ! BCG2
     &                                        0.45 /) ! BCG1

      real, parameter :: cstempS_3dbs(NSOAP) =           ! Ref. temp. for sat. conc. (K)
     &                                     (/ 300.,   ! ACG2
     &                                        300.,   ! ACG1
     &                                        300.,   ! BCG2
     &                                        300. /) ! BCG1

      real, parameter :: deltahS_3dbs(NSOAP) =           ! Heat of vap. for CG/SOA pairs (J/mol)
     &                                  (/ 116000.,   ! ACG2
     &                                     147000.,   ! ACG1
     &                                     118000.,   ! BCG2
     &                                     123000. /) ! BCG1
      real, parameter :: kpolya_3dbs   = 0.0   ! polymerization rate for AOA0 [1/hr]; deactivated
      real, parameter :: kpolyb_3dbs   = 0.0   ! polymerization rate for BOA0 [1/hr]; deactivated
      real, parameter :: sfacJSOA_3dbs = 0.001 ! scaling factor for JSOA: JSOA = sfacJSOA * JNO2

c     pointers iksoap() to access rate constants rk()
c            1 - BNZA + OH; 2 - TOLA + OH; 3 - XYLA + OH
c            4 - ISP + OH;  5 - ISP + O3;  6 - ISP + NO3
c            7 - TRP + OH;  8 - TRP + O3;  9 - TRP + NO3
c           10 - SQT + OH; 11 - SQT + O3; 12 - SQT + NO3
c           13 - APN + OH; 14 - APN + O3; 15 - APN + NO3

      integer, dimension(NSOAPRXN) :: iksoap

      real :: kpolya, kpolyb, sfacJSOA
      real, dimension(NSOAP) :: csatS, cstempS, deltahS
      real, dimension(3,NPREC_S) :: y_h, y_l, y_n

      logical luse_soap3_complx
      logical luse_soap3_simple
      logical luse_soap3

      common /soapcom/ kpolya, kpolyb, sfacJSOA, csatS, cstempS, deltahS,
     &                 y_h, y_l, y_n, luse_soap3_complx, luse_soap3_simple,
     &                 luse_soap3, iksoap

c-----------------------------------------------------------------------
c     Parameters for VBS
c-----------------------------------------------------------------------
c
      integer, parameter :: NPREC_V = 10
      integer, parameter :: NVBS    = 4   ! No. of volatility bins (volatile bins only)
      integer, parameter :: NSET    = 5   ! No. of basis sets

      real, parameter :: mwPAS0 = 172.    ! MW of PAS0
      real, parameter :: mwPBS0 = 172.    ! MW of PBS0
      real, parameter :: mwPAP0 = 278.    ! MW of PAP0
      real, parameter :: mwPCP0 = 278.    ! MW of PCP0
      real, parameter :: mwPFP0 = 205.    ! MW of PFP0

      real, parameter :: mwVBS(NVBS*NSET) = 
     &                                   (/ 167., 163., 158., 153.,   ! PAS1-PAS4
     &                                      167., 163., 158., 153.,   ! PBS1-PBS4
     &                                      275., 272., 268., 266.,   ! PAP1-PAP4
     &                                      275., 272., 268., 266.,   ! PCP1-PCP4
     &                                      208., 211., 213., 215. /) ! PFP1-PFP4

      real, parameter :: csatV(NVBS*NSET) =
     &                                   (/ 1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000. /)

      real, parameter :: cstempV(NVBS*NSET) =
     &                                   (/ 298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298. /)

      real, parameter :: deltahV(NVBS*NSET) =
     &                             (/ 35000., 35000., 35000., 35000.,
     &                                35000., 35000., 35000., 35000.,
     &                                85000., 74000., 63000., 52000.,
     &                                85000., 74000., 63000., 52000.,
     &                                85000., 81000., 77000., 73000. /)

c-----------------------------------------------------------------------
c     Parameters for common arrays
c-----------------------------------------------------------------------
      integer, parameter :: NOXID = 6                     ! No. of oxidant used in precursor chemistry
      integer, parameter :: NPREC = MAX(NPREC_S, NPREC_V) ! No. of SOA precursors
      integer, parameter :: NCG = MAX(NSOAP, NVBS*NSET)   ! No. of CG species
      integer, parameter :: NNV = MAX(NNVOL, NSET)        ! No. of non-volatile products

