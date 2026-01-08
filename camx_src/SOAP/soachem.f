      subroutine soachem(iout,tempk,press,dthr,rk,
     &                   avgox,cprec,dcprec,cCG,cNV,yfac,
     &                   nfam,avgsen,sen_prec,sen_CG,sen_NV,lddm)
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine solves the oxidation reactions for SOA precursors
c     assuming first-order decay and constant oxidant concentrations
c
c
c      Copyright 1996 - 2025
c     Ramboll
c
c     Input arguments:
c        iout     standard output file unit
c        tempk    cell temperature [K]
c        press    cell pressure [mb]
c        dthr     time duration to be integrated [hr]
c        rk       gas-phase rate constants [ppm-n hr-1]
c        avgox    average oxidant concentrations [ppm]
c                 1 - O; 2 - OH; 3 - O3; 4 - NO3; 5 - NO; 6 - HO2
c        cprec    precursor concentrations [ppm]
c                 1 - BNZA; 2 - TOLA; 3 - XYLA; 4 - IVOA
c                 5 - SVOA; 6 - ISP ; 7 - TRP ; 8 - APN; 9 - SQT
c        cCG      CG species concentrations [ppm]
c                 1 - ACG2; 2 - ACG1; 3 - BCG2; 4 - BCG1
c        cNV      NV products concentration [ppm]
c                 1 - CGA (anthro); 2 - CGB (bio)
c        yfac     weighted yield factors
c        nfam     number of sensitivity families
c        avgsen   average sensitivities of oxidants [ppm]
c        sen_prec sensitivities of precursors [ppm]
c        sen_CG   sensitivities of CGs [ppm]
c        sen_NV   sensitivities of NVs [ppm]
c        lddm     logical flag for DDM sensitivities
c
c     Output arguments:
c        cprec    precursor concentrations [ppm]
c        dcprec   change in precursor concentrations [ppm]
c        cCG      CG species concentrations [ppm]
c        cNV      NV products concentration [ppm]
c        yfac     weighted yield factors
c        sen_prec sensitivities of precursors [ppm]
c        sen_CG   sensitivities of CGs [ppm]
c        sen_NV   sensitivities of NVs [ppm]
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     03/22/06  --bkoo--  Original development
c     06/20/06  --bkoo--  2-product model for isoprene
c     08/23/13  --bkoo--  Added code for DDM-PM
c     03/18/14  --bkoo--  Revised for benzene SOA
c     08/25/16  --bkoo--  Updated for new SOAP (SOAP2 begins)
c     04/22/19  --gy----  Updated hi/low NOx rate constants
c     04/05/24  --cemery-- Added SVOA (SOAP3 begins)
c     05/16/24  --lhuang-- Differentiate SOA yields from TERP and
c                          add SOA yields from APIN 
c     05/23/24  --gy----  Use rate constats from the gas-phase mechanism
c                         except IVOA & SVOA rate constants are local and
c                         RO2 rate constants are local and follow MCM
c
c-----------------------------------------------------------------------
c
      implicit none
      include 'soap.inc'
c
c     soap.inc passes pointers iksoap() to access rate constants rk()
c            1 - BNZA + OH; 2 - TOLA + OH; 3 - XYLA + OH 
c            4 - ISP + OH;  5 - ISP + O3;  6 - ISP + NO3
c            7 - TRP + OH;  8 - TRP + O3;  9 - TRP + NO3
c           10 - SQT + OH; 11 - SQT + O3; 12 - SQT + NO3
c           13 - APN + OH; 14 - APN + O3; 15 - APN + NO3
c
      integer :: iout
      real    :: tempk,press,dthr
      real    :: avgox(NOXID),cprec(NPREC_S),cCG(NSOAP),cNV(NNVOL)
      real    :: dcprec(NPREC_S)
      real    :: yfac(3,NPREC_S)

      real    :: rk(*)

      integer, parameter :: idO   = 1, idOH  = 2, idO3  = 3,
     &                      idNO3 = 4, idNO  = 5, idHO2 = 6
      integer, parameter :: idBNZA =  1, idTOLA =  2,
     &                      idXYLA =  3, idIVOA =  4,
     &                      idSVOA =  5, idISP  =  6,
     &                      idTRP  =  7, idAPN  =  8,
     &                      idSQT  =  9

      real, parameter    :: eps = 1.e-12
c
c======================== DDM Begin =======================
c
      integer :: nfam,ifam
      real    :: avgsen(nfam,NOXID)
      real    :: sen_prec(nfam,NPREC)
      real    :: sen_CG(nfam,NCG)
      real    :: sen_NV(nfam,NNV)
      real    :: srh(nfam)
      real    :: ssumk,srk1(NOXID),sdprec
      logical :: lddm
c
c======================== DDM End =======================
c
c     [1/ppm-hr] = [cm3/molecule-sec]*3600*(6.022e23*1.e-12/8.314)*(P[Pa]/T[K])
      real, parameter  :: cf0 = 2.608e14 ! 3600*(6.022e23*1.e-12/8.314)
      real, parameter  :: rkivoa = 1.34E-11 ! OH + IVOA rate constant
      real, parameter  :: rksvoa = 1.6E-11  ! OH + SVOA rate constant
      real    :: cfac
      real    :: rk1(NOXID), rh, rl, fh, fl, fn, sumk, dprec, cfrac
      real    :: dfac
      integer :: i, j, k
c
c-----Conversion factor for 2nd-order rate constants from cms to ppm-hr
c
      cfac = cf0 * (press*100./tempk) ! [1/ppm-hr] = cfac * [cm3/molecule-sec]
c
c-----Re-set NV product concentrations to zero
c
      cNV = 0.0
c
c-----NOx-dependent pathways for aromatics
c
      rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )   ! aromatic RO2 + NO rate (hr-1)
      rl = 2.39e-12*exp(1300.0/tempk)*cfac * avgox(idHO2)  ! aromatic RO2 + HO2 rate (hr-1)
      fh = rh / (rh + rl)
      fl = 1.0 - fh

      yfac(1,idBNZA) = y_h(1,idBNZA)*fh + y_l(1,idBNZA)*fl
      yfac(2,idBNZA) = y_h(2,idBNZA)*fh + y_l(2,idBNZA)*fl
      yfac(3,idBNZA) = y_h(3,idBNZA)*fh + y_l(3,idBNZA)*fl

      yfac(1,idTOLA) = y_h(1,idTOLA)*fh + y_l(1,idTOLA)*fl
      yfac(2,idTOLA) = y_h(2,idTOLA)*fh + y_l(2,idTOLA)*fl
      yfac(3,idTOLA) = y_h(3,idTOLA)*fh + y_l(3,idTOLA)*fl

      yfac(1,idXYLA) = y_h(1,idXYLA)*fh + y_l(1,idXYLA)*fl
      yfac(2,idXYLA) = y_h(2,idXYLA)*fh + y_l(2,idXYLA)*fl
      yfac(3,idXYLA) = y_h(3,idXYLA)*fh + y_l(3,idXYLA)*fl
c
c======================== DDM Begin =======================
c
      if ( lddm ) then
        do ifam = 1, nfam
          sen_NV(ifam,:) = 0.0 ! <- cNV reset
          srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
        enddo
      endif
c
c======================== DDM End =======================
c
      rk1(idOH ) = rk(iksoap(1)) * avgox(idOH ) ! BNZA + OH rate (hr-1)
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idBNZA) * cfrac
      dcprec(idBNZA) = -dprec
      if (dprec.gt.eps) then
        cprec(idBNZA) = cprec(idBNZA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idBNZA)
        cCG(2) = cCG(2) + dprec * yfac(2,idBNZA)
        cNV(1) = cNV(1) + dprec * yfac(3,idBNZA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rk(iksoap(1)) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idBNZA) * cfrac
     &             + cprec(idBNZA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idBNZA) = sen_prec(ifam,idBNZA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idBNZA)
     &                                      + dprec *
     &            ( y_h(1,idBNZA)*srh(ifam) - y_l(1,idBNZA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idBNZA)
     &                                      + dprec *
     &            ( y_h(2,idBNZA)*srh(ifam) - y_l(2,idBNZA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idBNZA)
     &                                      + dprec *
     &            ( y_h(3,idBNZA)*srh(ifam) - y_l(3,idBNZA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

      rk1(idOH ) = rk(iksoap(2)) * avgox(idOH ) ! TOLA + OH rate (hr-1)
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idTOLA) * cfrac
      dcprec(idTOLA) = -dprec
      if (dprec.gt.eps) then
        cprec(idTOLA) = cprec(idTOLA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idTOLA)
        cCG(2) = cCG(2) + dprec * yfac(2,idTOLA)
        cNV(1) = cNV(1) + dprec * yfac(3,idTOLA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rk(iksoap(2)) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idTOLA) * cfrac
     &             + cprec(idTOLA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idTOLA) = sen_prec(ifam,idTOLA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idTOLA)
     &                                      + dprec *
     &            ( y_h(1,idTOLA)*srh(ifam) - y_l(1,idTOLA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idTOLA)
     &                                      + dprec *
     &            ( y_h(2,idTOLA)*srh(ifam) - y_l(2,idTOLA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idTOLA)
     &                                      + dprec *
     &            ( y_h(3,idTOLA)*srh(ifam) - y_l(3,idTOLA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

      rk1(idOH ) = rk(iksoap(3)) * avgox(idOH ) ! XYLA + OH rate (hr-1)
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idXYLA) * cfrac
      dcprec(idXYLA) = -dprec
      if (dprec.gt.eps) then
        cprec(idXYLA) = cprec(idXYLA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idXYLA)
        cCG(2) = cCG(2) + dprec * yfac(2,idXYLA)
        cNV(1) = cNV(1) + dprec * yfac(3,idXYLA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rk(iksoap(3)) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idXYLA) * cfrac
     &             + cprec(idXYLA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idXYLA) = sen_prec(ifam,idXYLA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idXYLA)
     &                                      + dprec *
     &            ( y_h(1,idXYLA)*srh(ifam) - y_l(1,idXYLA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idXYLA)
     &                                      + dprec *
     &            ( y_h(2,idXYLA)*srh(ifam) - y_l(2,idXYLA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idXYLA)
     &                                      + dprec *
     &            ( y_h(3,idXYLA)*srh(ifam) - y_l(3,idXYLA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for isoprene
c
      rk1(idOH ) = rk(iksoap(4)) * avgox(idOH ) ! ISP + OH rate (hr-1)
      rk1(idO3 ) = rk(iksoap(5)) * avgox(idO3 ) ! ISP + O3 rate (hr-1)
      rk1(idNO3) = rk(iksoap(6)) * avgox(idNO3) ! ISP + NO3 rate (hr-1)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)

      rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )  ! isoprene RO2 + NO rate (hr-1)
      rl = 1.93e-12*exp(1300.0/tempk)*cfac * avgox(idHO2) ! isoprene RO2 + HO2 rate (hr-1)

      fn = rk1(idNO3) / sumk
      fh = (1.0 - fn) * rh / (rh + rl)
      fl = 1.0 - fh - fn

      yfac(1,idISP) = y_h(1,idISP)*fh + y_l(1,idISP)*fl + y_n(1,idISP)*fn
      yfac(2,idISP) = y_h(2,idISP)*fh + y_l(2,idISP)*fl + y_n(2,idISP)*fn
      yfac(3,idISP) = y_h(3,idISP)*fh + y_l(3,idISP)*fl + y_n(3,idISP)*fn

      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idISP) * cfrac
      dcprec(idISP) = -dprec
      if (dprec.gt.eps) then
        cprec(idISP) = cprec(idISP) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idISP)
        cCG(4) = cCG(4) + dprec * yfac(2,idISP)
        cNV(2) = cNV(2) + dprec * yfac(3,idISP)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
            ssumk = rk(iksoap(4)) * avgsen(ifam,idOH)
     &            + rk(iksoap(5)) * avgsen(ifam,idO3)
     &            + rk(iksoap(6)) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idISP) * cfrac
     &             + cprec(idISP) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idISP) = sen_prec(ifam,idISP) - sdprec

            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idISP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(1,idISP)*srh(ifam) - y_l(1,idISP)*srh(ifam) )
     &                                      + dprec * fn * y_n(1,idISP)

            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idISP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(2,idISP)*srh(ifam) - y_l(2,idISP)*srh(ifam) )
     &                                      + dprec * fn * y_n(2,idISP)

            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idISP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(3,idISP)*srh(ifam) - y_l(3,idISP)*srh(ifam) )
     &                                      + dprec * fn * y_n(3,idISP)
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for lumped monoterpenes
c
      rk1(idOH ) = rk(iksoap(7)) * avgox(idOH ) ! TRP + OH rate (hr-1)
      rk1(idO3 ) = rk(iksoap(8)) * avgox(idO3 ) ! TRP + O3 rate (hr-1)
      rk1(idNO3) = rk(iksoap(9)) * avgox(idNO3) ! TRP + NO3 rate (hr-1)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)

      rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )  ! terpene RO2 + NO rate (hr-1)
      rl = 2.60e-12*exp(1300.0/tempk)*cfac * avgox(idHO2) ! terpene RO2 + HO2 rate (hr-1)

      fn = rk1(idNO3) / sumk
      fh = (1.0 - fn) * rh / (rh + rl)
      fl = 1.0 - fh - fn

      yfac(1,idTRP) = y_h(1,idTRP)*fh + y_l(1,idTRP)*fl + y_n(1,idTRP)*fn
      yfac(2,idTRP) = y_h(2,idTRP)*fh + y_l(2,idTRP)*fl + y_n(2,idTRP)*fn
      yfac(3,idTRP) = y_h(3,idTRP)*fh + y_l(3,idTRP)*fl + y_n(3,idTRP)*fn

      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idTRP) * cfrac
      dcprec(idTRP) = -dprec
      if (dprec.gt.eps) then
        cprec(idTRP) = cprec(idTRP) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idTRP)
        cCG(4) = cCG(4) + dprec * yfac(2,idTRP)
        cNV(2) = cNV(2) + dprec * yfac(3,idTRP)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
            ssumk = rk(iksoap(7)) * avgsen(ifam,idOH)
     &            + rk(iksoap(8)) * avgsen(ifam,idO3)
     &            + rk(iksoap(9)) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idTRP) * cfrac
     &             + cprec(idTRP) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idTRP) = sen_prec(ifam,idTRP) - sdprec

            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idTRP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(1,idTRP)*srh(ifam) - y_l(1,idTRP)*srh(ifam) )
     &                                      + dprec * fn * y_n(1,idTRP)

            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idTRP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(2,idTRP)*srh(ifam) - y_l(2,idTRP)*srh(ifam) )
     &                                      + dprec * fn * y_n(2,idTRP)

            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idTRP)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(3,idTRP)*srh(ifam) - y_l(3,idTRP)*srh(ifam) )
     &                                      + dprec * fn * y_n(3,idTRP)
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for a-pinene, if present
c
      yfac(1,idAPN) = 0.
      yfac(2,idAPN) = 0.
      yfac(3,idAPN) = 0.
      if(iksoap(13) .GT. 0) then  ! check whether APN is present

        rk1(idOH ) = rk(iksoap(13)) * avgox(idOH ) ! APN + OH rate (hr-1)
        rk1(idO3 ) = rk(iksoap(14)) * avgox(idO3 ) ! APN + O3 rate (hr-1)
        rk1(idNO3) = rk(iksoap(15)) * avgox(idNO3) ! APN + NO3 rate (hr-1)
        sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)

        rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )  ! a-pinene RO2 + NO rate (hr-1)
        rl = 2.60e-12*exp(1300.0/tempk)*cfac * avgox(idHO2) ! a-pinene RO2 + HO2 rate (hr-1)

        fn = rk1(idNO3) / sumk
        fh = (1.0 - fn) * rh / (rh + rl)
        fl = 1.0 - fh - fn

        yfac(1,idAPN) = y_h(1,idAPN)*fh + y_l(1,idAPN)*fl + y_n(1,idAPN)*fn
        yfac(2,idAPN) = y_h(2,idAPN)*fh + y_l(2,idAPN)*fl + y_n(2,idAPN)*fn
        yfac(3,idAPN) = y_h(3,idAPN)*fh + y_l(3,idAPN)*fl + y_n(3,idAPN)*fn

        cfrac = 1. - exp(-sumk*dthr)
        dprec = cprec(idAPN) * cfrac
        dcprec(idAPN) = -dprec
        if (dprec.gt.eps) then
          cprec(idAPN) = cprec(idAPN) - dprec
          cCG(3) = cCG(3) + dprec * yfac(1,idAPN)
          cCG(4) = cCG(4) + dprec * yfac(2,idAPN)
          cNV(2) = cNV(2) + dprec * yfac(3,idAPN)
c
c======================== DDM Begin =======================
c
          if ( lddm ) then
            do ifam = 1, nfam
              srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
              ssumk = rk(iksoap(13)) * avgsen(ifam,idOH)
     &              + rk(iksoap(14)) * avgsen(ifam,idO3)
     &              + rk(iksoap(15)) * avgsen(ifam,idNO3)
              sdprec = sen_prec(ifam,idAPN) * cfrac
     &               + cprec(idAPN) * dthr * ssumk ! because cprec has been updated
              sen_prec(ifam,idAPN) = sen_prec(ifam,idAPN) - sdprec

              sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idAPN)
     &                                        + dprec * (1.0 - fn) *
     &               ( y_h(1,idAPN)*srh(ifam) - y_l(1,idAPN)*srh(ifam) )
     &                                        + dprec * fn * y_n(1,idAPN)

              sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idAPN)
     &                                        + dprec * (1.0 - fn) *
     &               ( y_h(2,idAPN)*srh(ifam) - y_l(2,idAPN)*srh(ifam) )
     &                                        + dprec * fn * y_n(2,idAPN)

              sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idAPN)
     &                                        + dprec * (1.0 - fn) *
     &               ( y_h(3,idAPN)*srh(ifam) - y_l(3,idAPN)*srh(ifam) )
     &                                        + dprec * fn * y_n(3,idAPN)
            enddo
          endif
c
c======================== DDM End =======================
c
        endif
      endif  ! end check on whether APN is present
c
c-----NOx-dependent pathways for sesquiterpenes
c     use local rate constants if unavailable from gas-phase mechanism
c
      if( iksoap(10) .GT. 0 ) rk1(idOH ) = rk(iksoap(10)) * avgox(idOH ) ! SQT + OH rate (hr-1)
      if( iksoap(11) .GT. 0 ) rk1(idO3 ) = rk(iksoap(11)) * avgox(idO3 ) ! SQT + O3 rate (hr-1)
      if( iksoap(12) .GT. 0 ) rk1(idNO3) = rk(iksoap(12)) * avgox(idNO3) ! SQT + NO3 rate (hr-1)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)

      rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )  ! sesquiterpene RO2 + NO rate (hr-1)
      rl = 2.75e-12*exp(1300.0/tempk)*cfac * avgox(idHO2) ! sesquiterpene RO2 + HO2 rate (hr-1)

      fn = rk1(idNO3) / sumk
      fh = (1.0 - fn) * rh / (rh + rl)
      fl = 1.0 - fh - fn

      yfac(1,idSQT) = y_h(1,idSQT)*fh + y_l(1,idSQT)*fl + y_n(1,idSQT)*fn
      yfac(2,idSQT) = y_h(2,idSQT)*fh + y_l(2,idSQT)*fl + y_n(2,idSQT)*fn
      yfac(3,idSQT) = y_h(3,idSQT)*fh + y_l(3,idSQT)*fl + y_n(3,idSQT)*fn

      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idSQT) * cfrac
      dcprec(idSQT) = -dprec
      if (dprec.gt.eps) then
        cprec(idSQT) = cprec(idSQT) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idSQT)
        cCG(4) = cCG(4) + dprec * yfac(2,idSQT)
        cNV(2) = cNV(2) + dprec * yfac(3,idSQT)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
c
            ssumk = rk(iksoap(10)) * avgsen(ifam,idOH)
     &            + rk(iksoap(11)) * avgsen(ifam,idO3)
     &            + rk(iksoap(12)) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idSQT) * cfrac
     &             + cprec(idSQT) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idSQT) = sen_prec(ifam,idSQT) - sdprec
            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idSQT)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(1,idSQT)*srh(ifam) - y_l(1,idSQT)*srh(ifam) )
     &                                      + dprec * fn * y_n(1,idSQT)

            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idSQT)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(2,idSQT)*srh(ifam) - y_l(2,idSQT)*srh(ifam) )
     &                                      + dprec * fn * y_n(2,idSQT)

            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idSQT)
     &                                      + dprec * (1.0 - fn) *
     &             ( y_h(3,idSQT)*srh(ifam) - y_l(3,idSQT)*srh(ifam) )
     &                                      + dprec * fn * y_n(3,idSQT)
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for IVOA and SVOA
c
      rh = 2.70e-12*exp(360.0/tempk)*cfac * avgox(idNO )  ! IVOA RO2 + NO rate (hr-1)
      rl = 2.75e-12*exp(1300.0/tempk)*cfac * avgox(idHO2) ! IVOA RO2 + HO2 rate (hr-1)
      fh = rh / (rh + rl)
      fl = 1.0 - fh

      yfac(1,idIVOA) = y_h(1,idIVOA)*fh + y_l(1,idIVOA)*fl
      yfac(2,idIVOA) = y_h(2,idIVOA)*fh + y_l(2,idIVOA)*fl
      yfac(3,idIVOA) = y_h(3,idIVOA)*fh + y_l(3,idIVOA)*fl
c
      yfac(1,idSVOA) = y_h(1,idSVOA)*fh + y_l(1,idSVOA)*fl
      yfac(2,idSVOA) = y_h(2,idSVOA)*fh + y_l(2,idSVOA)*fl
      yfac(3,idSVOA) = y_h(3,idSVOA)*fh + y_l(3,idSVOA)*fl
c
c======================== DDM Begin =======================
c
      if ( lddm ) then
        do ifam = 1, nfam
          srh(ifam) = (fl*rh - fh*rl) / (rh + rl)
        enddo
      endif
c
c======================== DDM End =======================
c
      rk1(idOH ) = rkivoa * cfac * avgox(idOH ) ! IVOA + OH rate (hr-1)
      cfrac = 1. - exp(-rk1(idOH )*dthr)

      dprec = cprec(idIVOA) * cfrac
      dcprec(idIVOA) = -dprec
      if (dprec.gt.eps) then
        cprec(idIVOA) = cprec(idIVOA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idIVOA)
        cCG(2) = cCG(2) + dprec * yfac(2,idIVOA)
        cNV(1) = cNV(1) + dprec * yfac(3,idIVOA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rkivoa * cfac * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idIVOA) * cfrac
     &             + cprec(idIVOA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idIVOA) = sen_prec(ifam,idIVOA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idIVOA)
     &                                      + dprec *
     &            ( y_h(1,idIVOA)*srh(ifam) - y_l(1,idIVOA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idIVOA)
     &                                      + dprec *
     &            ( y_h(2,idIVOA)*srh(ifam) - y_l(2,idIVOA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idIVOA)
     &                                      + dprec *
     &            ( y_h(3,idIVOA)*srh(ifam) - y_l(3,idIVOA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

      rk1(idOH ) = rksvoa * cfac * avgox(idOH ) ! SVOA + OH rate (hr-1)
      cfrac = 1. - exp(-rk1(idOH )*dthr)

      dprec = cprec(idSVOA) * cfrac
      dcprec(idSVOA) = -dprec
      if (dprec.gt.eps) then
        cprec(idSVOA) = cprec(idSVOA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idSVOA)
        cCG(2) = cCG(2) + dprec * yfac(2,idSVOA)
        cNV(1) = cNV(1) + dprec * yfac(3,idSVOA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rksvoa * cfac * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idSVOA) * cfrac
     &             + cprec(idSVOA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idSVOA) = sen_prec(ifam,idSVOA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idSVOA)
     &                                      + dprec *
     &            ( y_h(1,idSVOA)*srh(ifam) - y_l(1,idSVOA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idSVOA)
     &                                      + dprec *
     &            ( y_h(2,idSVOA)*srh(ifam) - y_l(2,idSVOA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idSVOA)
     &                                      + dprec *
     &            ( y_h(3,idSVOA)*srh(ifam) - y_l(3,idSVOA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----Save the weighting factors for NOx-dependent yields for PSAT
c       We assume same rates for ARO-O2 and IVOA-O2
c       If not, PSAT needs to be re-worked
c
      yfac(1,1) = fh
      yfac(2,1) = fl
c
      return
      end

