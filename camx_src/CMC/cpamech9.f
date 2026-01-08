      subroutine cpamech9(r, rk, dtfact, nr, pa, npa, npa_init, ldark )

      use filunit
      use tracer
      use procan
      implicit none
c
c----CAMx v7.32 250801
c
c     CPAMECH9 computes chemical process analysis (CPA) parameters
c
c     Copyright 1996 - 2025
c     Ramboll
c     Created by the CMC version 5.2.7
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        PASETUP
c        CHEMDRIV
c
c --- Argument definitions:
c        r        - reaction rates (ppm/hr)
c        rk       - rate constants (ppm-n hr-1)
c        dtfact   - ratio of time step to averaging interval
c        nr       - number of reactions
c        pa       - output cpa parameter values
c        npa      - dimension of cpa parameter arrays
c        npa_init - number of cpa parameters calculated below
c        ldark    - flag set true at night
c
c --- Includes:
      include "camx.prm"
c
c --- Arguments:
c
      integer  nr, npa, npa_init
      real     r(nr), rk(nr), pa(MXCPA)
      real     dtfact
      logical  ldark
c
c --- Parameters:
c
c     Cut-point between NOx and VOC sensitive O3 production
c     Default = 0.35.  Recommended range is 0.15 to 0.35
      real     ratio_cut
      parameter (ratio_cut = 0.35)
c
c     Convert PPM to PPB
      real     ppbfact
      parameter (ppbfact = 1000.0)
c
c --- Local variables:
c
      integer  n, nn
      real     sum, sun, po3, do3, HO2toNO2, HO2_loss
      real     NO2wOH, HO2wHO2, ratio_ind
      real     OH_new, HO2_new, OH_loss, RO2_loss, HOx_CL
      real     newOH_O1D, newOH_O3, newOH_phot
      real     newHO2_O3, newHO2_pht
      real     NO3wVOC, N2O5wH2O, PAN_prdNet, ON_prod
      real     NOxNet
c
c --- Entry point:
c
      nn = 0
      sun = 1.0
      if (ldark) sun = 0.0
c
c --- J(NO2) photolysis rate (per hour)
c
      nn = nn + 1
      ptname(nn)  = 'J_NO2'
      cpadesc(nn) = 'J(NO2) photolysis rate'
      cpaunit(nn) = 'hr-1'
      PA(nn) =      rk(  1)*(dtfact/ppbfact)
c
c --- J(O3O1D) photolysis rate (per hour)
c
      nn = nn + 1
      ptname(nn)  = 'J_O3O1D'
      cpadesc(nn) = 'J(O3) to O(1D) photolysis rate'
      cpaunit(nn) = 'hr-1'
      PA(nn) =      rk(  3)*(dtfact/ppbfact)
c
c --- J(HCHOr) photolysis rate (per hour)
c
      nn = nn + 1
      ptname(nn)  = 'J_HCHOr'
      cpadesc(nn) = 'J(HCHO) to 2 HO2 + CO photolysis rate'
      cpaunit(nn) = 'hr-1'
      PA(nn) =      rk( 11)*(dtfact/ppbfact)
c
c
c --- Net O3 production
c
      nn = nn + 1
      ptname(nn)  = 'PO3_net'
      cpadesc(nn) = 'Net O3 production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 38)  ! O + O2 + M =
     &       -          r(  2)  ! O3 =
     &       -          r(  3)  ! O3 =
     &       -          r( 34)  ! O3 + OH =
     &       -          r( 35)  ! O3 + HO2 =
     &       -          r( 36)  ! O3 + NO =
     &       -          r( 37)  ! O3 + NO2 =
     &       -          r( 39)  ! O + O3 =
     &       -          r(126)  ! ETE + O3 =
     &       -          r(127)  ! OLT + O3 =
     &       -          r(128)  ! OLI + O3 =
     &       -          r(129)  ! DIEN + O3 =
     &       -          r(130)  ! ISO + O3 =
     &       -          r(131)  ! API + O3 =
     &       -          r(132)  ! LIM + O3 =
     &       -          r(133)  ! MACR + O3 =
     &       -          r(134)  ! MVK + O3 =
     &       -          r(135)  ! UALD + O3 =
     &       -          r(136)  ! DCB1 + O3 =
     &       -          r(137)  ! DCB2 + O3 =
     &       -          r(138)  ! DCB3 + O3 =
     &       -          r(139)  ! EPX + O3 =
     &       -          r(140)  ! MCTO + O3 =
c
      po3 = max(0.0, PA(nn))
      do3 = min(0.0, PA(nn))
c
c --- Calculate the P(H2O2)/P(HNO3) indicator ratio and apply to PO3
c
      HO2wHO2 = + r( 45) + r( 46)
      NO2wOH = + r( 56)
      ratio_ind = min( 10., HO2wHO2/max( 1.0E-12, NO2wOH ) )*sun
c
      nn =  nn + 2
      ptname(nn-1)  = 'PO3_VOCsns'
      cpadesc(nn-1) = 'Net O3 production rate under VOC-limited condition'
      cpaunit(nn-1) = 'ppb hr-1'
      ptname(nn)  = 'PO3_NOxsns'
      cpadesc(nn) = 'Net O3 production rate under NOx-limited condition'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn-1) = 0.0
      PA(nn)   = 0.0
      if (ratio_ind .GT. 0.0) then
        if (ratio_ind .LT. ratio_cut ) then
          PA(nn-1) = PO3
        else
          PA(nn)   = PO3
        endif
      endif
c
c --- Report the P(H2O2)/P(HNO3) indicator ratio
c
      if(.NOT. lcpacum) then
        nn =  nn + 1
        ptname(nn)   = 'PH2O2_PHN3'
        cpadesc(nn) = 'Ratio of production rates for H2O2/HNO3'
        cpaunit(nn) = 'dimensionless'
        if (ldark) then
          PA(nn) = 0.0
        else
          PA(nn) = ratio_ind * (dtfact/ppbfact)
        endif
      endif
c
c --- Ozone destruction calculation from OSAT
c
      nn = nn + 1
      ptname(nn)  = 'O3_dest'
      cpadesc(nn) = 'O3 destruction rate (same implementation as OSAT)'
      cpaunit(nn) = 'ppb hr-1'
c
c  +  O1D + H2O
c
      PA(nn) =
     &       - r( 42)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
      PA(nn) = PA(nn)
     &       - r( 35)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
      HO2toNO2 = 
     &       +          r( 50)  ! NO + HO2 =
     &       + ( 0.700)*r( 59)  ! NO3 + HO2 =
c
      HO2_loss = 
     &       +          r( 35)  ! O3 + HO2 =
     &       +          r( 44)  ! OH + HO2 =
     &       + ( 2.000)*r( 45)  ! HO2 + HO2 =
     &       + ( 2.000)*r( 46)  ! HO2 + HO2 + H2O =
     &       +          r( 50)  ! NO + HO2 =
     &       +          r( 51)  ! NO + HO2 =
     &       +          r( 59)  ! NO3 + HO2 =
     &       +          r(212)  ! MO2 + HO2 =
     &       +          r(213)  ! ETHP + HO2 =
     &       +          r(214)  ! HC3P + HO2 =
     &       +          r(215)  ! HC5P + HO2 =
     &       +          r(216)  ! HC8P + HO2 =
     &       +          r(217)  ! ETEP + HO2 =
     &       +          r(218)  ! OLTP + HO2 =
     &       +          r(219)  ! OLIP + HO2 =
     &       +          r(220)  ! BENP + HO2 =
     &       +          r(221)  ! TLP1 + HO2 =
     &       +          r(222)  ! TOLP + HO2 =
     &       +          r(223)  ! PER1 + HO2 =
     &       +          r(224)  ! XYL1 + HO2 =
     &       +          r(225)  ! XYLP + HO2 =
     &       +          r(226)  ! PER2 + HO2 =
     &       +          r(227)  ! XYOP + HO2 =
     &       +          r(228)  ! ISOP + HO2 =
     &       +          r(229)  ! APIP + HO2 =
     &       +          r(230)  ! LIMP + HO2 =
     &       +          r(231)  ! ACO3 + HO2 =
     &       +          r(232)  ! RCO3 + HO2 =
     &       +          r(233)  ! ACTP + HO2 =
     &       +          r(234)  ! MEKP + HO2 =
     &       +          r(235)  ! KETP + HO2 =
     &       +          r(236)  ! MACP + HO2 =
     &       +          r(237)  ! MCP + HO2 =
     &       +          r(238)  ! MVKP + HO2 =
     &       +          r(239)  ! UALP + HO2 =
     &       +          r(240)  ! ADDC + HO2 =
     &       +          r(241)  ! CHO + HO2 =
     &       +          r(242)  ! MCTP + HO2 =
     &       +          r(243)  ! ORAP + HO2 =
     &       +          r(244)  ! OLNN + HO2 =
     &       +          r(245)  ! OLND + HO2 =
     &       +          r(246)  ! ADCN + HO2 =
     &       +          r(247)  ! XO2 + HO2 =
c
      PA(nn) = PA(nn) - (
     &       + r( 34)  ! O3 + OH =
     &                     ) * (HO2_loss-HO2toNO2)/max( 1.0E-12, HO2_loss )
c
c  +  O3 + VOC
c
      PA(nn) = PA(nn)
     &       - r(126)  ! ETE + O3 =
     &       - r(127)  ! OLT + O3 =
     &       - r(128)  ! OLI + O3 =
     &       - r(129)  ! DIEN + O3 =
     &       - r(130)  ! ISO + O3 =
     &       - r(131)  ! API + O3 =
     &       - r(132)  ! LIM + O3 =
     &       - r(133)  ! MACR + O3 =
     &       - r(134)  ! MVK + O3 =
     &       - r(135)  ! UALD + O3 =
     &       - r(136)  ! DCB1 + O3 =
     &       - r(137)  ! DCB2 + O3 =
     &       - r(138)  ! DCB3 + O3 =
     &       - r(139)  ! EPX + O3 =
     &       - r(140)  ! MCTO + O3 =
c
c  +  O(3P) + VOC
c
      PA(nn) = PA(nn)
c
      PA(nn) = min(PA(nn), do3)
c
c --- Net NO production
c
      nn = nn + 1
      ptname(nn)  = 'PNO_net'
      cpadesc(nn) = 'Net NO production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  1)  ! NO2 =
     &       +          r(  5)  ! NO3 =
     &       +          r(  7)  ! HONO =
     &       +          r( 54)  ! NO2 + O =
     &       +          r( 61)  ! NO3 + NO2 =
     &       -          r( 36)  ! O3 + NO =
     &       -          r( 48)  ! NO + O =
     &       -          r( 49)  ! NO + OH =
     &       -          r( 50)  ! NO + HO2 =
     &       -          r( 51)  ! NO + HO2 =
     &       + ( 2.000)*r( 52)  ! NO + NO + O2 =
     &       -          r( 60)  ! NO3 + NO =
     &       -          r(172)  ! MO2 + NO =
     &       -          r(173)  ! ETHP + NO =
     &       -          r(174)  ! HC3P + NO =
     &       -          r(175)  ! HC5P + NO =
     &       -          r(176)  ! HC8P + NO =
     &       -          r(177)  ! ETEP + NO =
     &       -          r(178)  ! OLTP + NO =
     &       -          r(179)  ! OLIP + NO =
     &       -          r(180)  ! BENP + NO =
     &       -          r(181)  ! TLP1 + NO =
     &       -          r(182)  ! TOLP + NO =
     &       -          r(183)  ! PER1 + NO =
     &       -          r(184)  ! XYL1 + NO =
     &       -          r(185)  ! XYLP + NO =
     &       -          r(186)  ! PER2 + NO =
     &       -          r(187)  ! XYOP + NO =
     &       -          r(188)  ! ISOP + NO =
     &       -          r(189)  ! APIP + NO =
     &       -          r(190)  ! LIMP + NO =
     &       -          r(191)  ! ACO3 + NO =
     &       -          r(192)  ! RCO3 + NO =
     &       -          r(193)  ! ACTP + NO =
     &       -          r(194)  ! MEKP + NO =
     &       -          r(195)  ! KETP + NO =
     &       -          r(196)  ! MACP + NO =
     &       -          r(197)  ! MCP + NO =
     &       -          r(198)  ! MVKP + NO =
     &       -          r(199)  ! UALP + NO =
     &       -          r(200)  ! BALP + NO =
     &       -          r(201)  ! BAL1 + NO =
     &       -          r(202)  ! ADDC + NO =
     &       -          r(203)  ! MCTP + NO =
     &       -          r(204)  ! ORAP + NO =
     &       -          r(205)  ! OLNN + NO =
     &       -          r(206)  ! OLND + NO =
     &       -          r(207)  ! ADCN + NO =
     &       -          r(208)  ! XO2 + NO =

      NOxNet = PA(nn)
c
c --- Net NO2 production
c
      nn = nn + 1
      ptname(nn)  = 'PNO2_net'
      cpadesc(nn) = 'Net NO2 production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  6)  ! NO3 =
     &       +          r(  8)  ! HNO3 =
     &       + ( 0.800)*r(  9)  ! HNO4 =
     &       +          r( 31)  ! ONIT =
     &       +          r( 32)  ! PAN =
     &       +          r( 36)  ! O3 + NO =
     &       +          r( 48)  ! NO + O =
     &       +          r( 50)  ! NO + HO2 =
     &       + ( 2.000)*r( 52)  ! NO + NO + O2 =
     &       +          r( 53)  ! HONO + OH =
     &       +          r( 58)  ! NO3 + OH =
     &       + ( 0.700)*r( 59)  ! NO3 + HO2 =
     &       + ( 2.000)*r( 60)  ! NO3 + NO =
     &       + ( 0.000)*r( 61)  ! NO3 + NO2 =
     &       + ( 2.000)*r( 62)  ! NO3 + NO3 =
     &       +          r( 64)  ! N2O5 =
     &       +          r( 67)  ! HNO4 =
     &       +          r( 68)  ! HNO4 + OH =
     &       +          r(122)  ! MPAN + OH =
     &       +          r(123)  ! ONIT + OH =
     &       +          r(124)  ! NALD + OH =
     &       + ( 0.680)*r(151)  ! MACR + NO3 =
     &       + ( 0.500)*r(157)  ! EPX + NO3 =
     &       +          r(159)  ! MPAN + NO3 =
     &       +          r(167)  ! PAN =
     &       +          r(169)  ! PPN =
     &       +          r(171)  ! MPAN =
     &       +          r(172)  ! MO2 + NO =
     &       +          r(173)  ! ETHP + NO =
     &       + ( 0.935)*r(174)  ! HC3P + NO =
     &       + ( 0.864)*r(175)  ! HC5P + NO =
     &       + ( 0.739)*r(176)  ! HC8P + NO =
     &       +          r(177)  ! ETEP + NO =
     &       + ( 0.970)*r(178)  ! OLTP + NO =
     &       + ( 0.950)*r(179)  ! OLIP + NO =
     &       + ( 0.918)*r(180)  ! BENP + NO =
     &       +          r(181)  ! TLP1 + NO =
     &       + ( 0.950)*r(182)  ! TOLP + NO =
     &       + ( 0.950)*r(183)  ! PER1 + NO =
     &       +          r(184)  ! XYL1 + NO =
     &       + ( 0.950)*r(185)  ! XYLP + NO =
     &       + ( 0.950)*r(186)  ! PER2 + NO =
     &       + ( 0.950)*r(187)  ! XYOP + NO =
     &       + ( 0.880)*r(188)  ! ISOP + NO =
     &       + ( 0.820)*r(189)  ! APIP + NO =
     &       +          r(190)  ! LIMP + NO =
     &       +          r(191)  ! ACO3 + NO =
     &       +          r(192)  ! RCO3 + NO =
     &       +          r(193)  ! ACTP + NO =
     &       +          r(194)  ! MEKP + NO =
     &       +          r(195)  ! KETP + NO =
     &       +          r(196)  ! MACP + NO =
     &       +          r(197)  ! MCP + NO =
     &       +          r(198)  ! MVKP + NO =
     &       +          r(199)  ! UALP + NO =
     &       +          r(200)  ! BALP + NO =
     &       +          r(201)  ! BAL1 + NO =
     &       +          r(202)  ! ADDC + NO =
     &       +          r(203)  ! MCTP + NO =
     &       +          r(204)  ! ORAP + NO =
     &       +          r(205)  ! OLNN + NO =
     &       + ( 2.000)*r(206)  ! OLND + NO =
     &       + ( 2.000)*r(207)  ! ADCN + NO =
     &       +          r(208)  ! XO2 + NO =
     &       +          r(273)  ! MCP + MO2 =
     &       + ( 0.500)*r(282)  ! OLND + MO2 =
     &       + ( 0.700)*r(283)  ! ADCN + MO2 =
     &       +          r(309)  ! MCP + ACO3 =
     &       +          r(318)  ! OLND + ACO3 =
     &       + ( 0.700)*r(319)  ! ADCN + ACO3 =
     &       +          r(322)  ! MO2 + NO3 =
     &       +          r(323)  ! ETHP + NO3 =
     &       +          r(324)  ! HC3P + NO3 =
     &       +          r(325)  ! HC5P + NO3 =
     &       +          r(326)  ! HC8P + NO3 =
     &       +          r(327)  ! ETEP + NO3 =
     &       +          r(328)  ! OLTP + NO3 =
     &       +          r(329)  ! OLIP + NO3 =
     &       +          r(330)  ! BENP + NO3 =
     &       +          r(331)  ! TLP1 + NO3 =
     &       +          r(332)  ! TOLP + NO3 =
     &       +          r(333)  ! PER1 + NO3 =
     &       +          r(334)  ! XYL1 + NO3 =
     &       +          r(335)  ! XYLP + NO3 =
     &       +          r(336)  ! PER2 + NO3 =
     &       +          r(337)  ! XYOP + NO3 =
     &       +          r(338)  ! ISOP + NO3 =
     &       +          r(339)  ! APIP + NO3 =
     &       +          r(340)  ! LIMP + NO3 =
     &       +          r(341)  ! ACO3 + NO3 =
     &       +          r(342)  ! RCO3 + NO3 =
     &       +          r(343)  ! ACTP + NO3 =
     &       +          r(344)  ! MEKP + NO3 =
     &       +          r(345)  ! KETP + NO3 =
     &       +          r(346)  ! MACP + NO3 =
     &       +          r(347)  ! MCP + NO3 =
     &       +          r(348)  ! MVKP + NO3 =
     &       +          r(349)  ! UALP + NO3 =
     &       +          r(350)  ! BALP + NO3 =
     &       +          r(351)  ! BAL1 + NO3 =
     &       +          r(352)  ! ADDC + NO3 =
     &       +          r(353)  ! MCTP + NO3 =
     &       +          r(354)  ! ORAP + NO3 =
     &       +          r(355)  ! OLNN + NO3 =
     &       + ( 2.000)*r(356)  ! OLND + NO3 =
     &       + ( 2.000)*r(357)  ! ADCN + NO3 =
     &       + ( 0.500)*r(359)  ! OLNN + OLND =
     &       +          r(360)  ! OLND + OLND =
     &       +          r(361)  ! XO2 + NO3 =
     &       -          r(  1)  ! NO2 =
     &       -          r( 37)  ! O3 + NO2 =
     &       -          r( 54)  ! NO2 + O =
     &       -          r( 55)  ! NO2 + O =
     &       -          r( 56)  ! NO2 + OH =
     &       + ( 0.000)*r( 61)  ! NO3 + NO2 =
     &       -          r( 63)  ! NO3 + NO2 =
     &       -          r( 66)  ! NO2 + HO2 =
     &       -          r(166)  ! ACO3 + NO2 =
     &       -          r(168)  ! RCO3 + NO2 =
     &       -          r(170)  ! MACP + NO2 =
     &       -          r(209)  ! BAL2 + NO2 =
     &       -          r(210)  ! CHO + NO2 =
     &       -          r(211)  ! MCTO + NO2 =

      NOxNet = NOxNet + PA(nn)
c
c --- Net NOx production
c
      nn = nn + 1
      ptname(nn)  = 'PNOx_net'
      cpadesc(nn) = 'Net NOx production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = NOxNet

c
c --- Total OH production
c
      nn = nn + 1
      ptname(nn)  = 'OH_prod'
      cpadesc(nn) = 'Total OH production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r(  4)  ! H2O2 =
     &       +          r(  7)  ! HONO =
     &       +          r(  8)  ! HNO3 =
     &       + ( 0.200)*r(  9)  ! HNO4 =
     &       + ( 0.340)*r( 19)  ! MACR =
     &       +          r( 28)  ! OP1 =
     &       +          r( 29)  ! OP2 =
     &       +          r( 30)  ! PAA =
     &       +          r( 35)  ! O3 + HO2 =
     &       + ( 2.000)*r( 42)  ! O1D + H2O =
     &       +          r( 50)  ! NO + HO2 =
     &       + ( 0.700)*r( 59)  ! NO3 + HO2 =
     &       + ( 0.080)*r(126)  ! ETE + O3 =
     &       + ( 0.220)*r(127)  ! OLT + O3 =
     &       + ( 0.460)*r(128)  ! OLI + O3 =
     &       + ( 0.280)*r(129)  ! DIEN + O3 =
     &       + ( 0.250)*r(130)  ! ISO + O3 =
     &       + ( 0.850)*r(131)  ! API + O3 =
     &       + ( 0.850)*r(132)  ! LIM + O3 =
     &       + ( 0.190)*r(133)  ! MACR + O3 =
     &       + ( 0.160)*r(134)  ! MVK + O3 =
     &       + ( 0.100)*r(135)  ! UALD + O3 =
     &       + ( 0.050)*r(136)  ! DCB1 + O3 =
     &       + ( 0.050)*r(137)  ! DCB2 + O3 =
     &       + ( 0.050)*r(138)  ! DCB3 + O3 =
     &       + ( 0.050)*r(139)  ! EPX + O3 =
     &       + ( 0.500)*r(157)  ! EPX + NO3 =
     &       + ( 0.280)*r(160)  ! TR2 =
     &       + ( 0.490)*r(161)  ! TOLP =
     &       + ( 0.158)*r(162)  ! XY2 =
     &       + ( 0.390)*r(163)  ! XYLP =
     &       + ( 0.158)*r(164)  ! XYO2 =
     &       + ( 0.390)*r(165)  ! XYOP =
     &       + ( 0.440)*r(231)  ! ACO3 + HO2 =
     &       + ( 0.440)*r(232)  ! RCO3 + HO2 =
     &       + ( 0.150)*r(233)  ! ACTP + HO2 =
c
c --- OH from O(1D)
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O1D'
      cpadesc(nn) = 'OH production rate from O(1D) + H2O'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = ( 2.000)*r( 42)  ! O1D + H2O =
      newOH_O1D = PA(nn)
c
c --- OH from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O3'
      cpadesc(nn) = 'OH production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.080)*r(126)  ! ETE + O3 =
     &       + ( 0.220)*r(127)  ! OLT + O3 =
     &       + ( 0.460)*r(128)  ! OLI + O3 =
     &       + ( 0.280)*r(129)  ! DIEN + O3 =
     &       + ( 0.250)*r(130)  ! ISO + O3 =
     &       + ( 0.850)*r(131)  ! API + O3 =
     &       + ( 0.850)*r(132)  ! LIM + O3 =
     &       + ( 0.190)*r(133)  ! MACR + O3 =
     &       + ( 0.160)*r(134)  ! MVK + O3 =
     &       + ( 0.100)*r(135)  ! UALD + O3 =
     &       + ( 0.050)*r(136)  ! DCB1 + O3 =
     &       + ( 0.050)*r(137)  ! DCB2 + O3 =
     &       + ( 0.050)*r(138)  ! DCB3 + O3 =
     &       + ( 0.050)*r(139)  ! EPX + O3 =
      newOH_O3 = PA(nn)
c
c --- OH directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newOH_phot'
      cpadesc(nn) = 'OH production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r(  4)  ! H2O2 =
     &       +          r(  7)  ! HONO =
     &       +          r(  8)  ! HNO3 =
     &       + ( 0.340)*r( 19)  ! MACR =
     &       +          r( 28)  ! OP1 =
     &       +          r( 29)  ! OP2 =
     &       +          r( 30)  ! PAA =
      newOH_phot = PA(nn)
c
c --- New OH
c
      nn = nn + 1
      ptname(nn)  = 'OH_new'
      cpadesc(nn) = 'New OH, i.e., newOH_O1D + newOH_O3 + newOH_phot'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = newOH_O1D + newOH_O3 + newOH_phot
      OH_new = PA(nn)
c
c --- OH from HONO
c
      nn = nn + 1
      ptname(nn)  = 'newOH_HONO'
      cpadesc(nn) = 'OH production rate from HONO photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =          r(  7)  ! HONO =
c
c --- OH loss
c
      nn = nn + 1
      ptname(nn)  = 'OH_loss'
      cpadesc(nn) = 'OH loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 34)  ! O3 + OH =
     &       +          r( 43)  ! OH + H2 =
     &       +          r( 44)  ! OH + HO2 =
     &       +          r( 47)  ! H2O2 + OH =
     &       +          r( 49)  ! NO + OH =
     &       +          r( 53)  ! HONO + OH =
     &       +          r( 56)  ! NO2 + OH =
     &       +          r( 57)  ! HNO3 + OH =
     &       +          r( 58)  ! NO3 + OH =
     &       +          r( 68)  ! HNO4 + OH =
     &       +          r( 69)  ! SO2 + OH =
     &       +          r( 70)  ! CO + OH =
     &       +          r( 71)  ! CH4 + OH =
     &       +          r( 72)  ! ETH + OH =
     &       +          r( 73)  ! HC3 + OH =
     &       +          r( 74)  ! HC5 + OH =
     &       +          r( 75)  ! HC8 + OH =
     &       +          r( 76)  ! ETE + OH =
     &       +          r( 77)  ! OLT + OH =
     &       +          r( 78)  ! OLI + OH =
     &       +          r( 79)  ! DIEN + OH =
     &       + ( 0.350)*r( 80)  ! ACE + OH =
     &       +          r( 81)  ! BEN + OH =
     &       +          r( 82)  ! TOL + OH =
     &       +          r( 83)  ! XYM + OH =
     &       +          r( 84)  ! XYP + OH =
     &       +          r( 85)  ! XYO + OH =
     &       +          r( 86)  ! ISO + OH =
     &       +          r( 87)  ! API + OH =
     &       +          r( 88)  ! LIM + OH =
     &       +          r( 89)  ! HCHO + OH =
     &       +          r( 90)  ! ACD + OH =
     &       +          r( 91)  ! ALD + OH =
     &       +          r( 92)  ! ACT + OH =
     &       +          r( 93)  ! MEK + OH =
     &       +          r( 94)  ! KET + OH =
     &       +          r( 95)  ! HKET + OH =
     &       +          r( 96)  ! MACR + OH =
     &       +          r( 97)  ! MVK + OH =
     &       +          r( 98)  ! UALD + OH =
     &       +          r( 99)  ! GLY + OH =
     &       +          r(100)  ! MGLY + OH =
     &       +          r(101)  ! DCB1 + OH =
     &       +          r(102)  ! DCB2 + OH =
     &       +          r(103)  ! DCB3 + OH =
     &       +          r(104)  ! BALD + OH =
     &       +          r(105)  ! PHEN + OH =
     &       +          r(106)  ! CSL + OH =
     &       +          r(107)  ! EPX + OH =
     &       +          r(108)  ! MCT + OH =
     &       +          r(109)  ! MOH + OH =
     &       +          r(110)  ! EOH + OH =
     &       +          r(111)  ! ROH + OH =
     &       +          r(112)  ! ETEG + OH =
     &       + ( 0.650)*r(113)  ! OP1 + OH =
     &       + ( 0.990)*r(114)  ! OP2 + OH =
     &       +          r(116)  ! MAHP + OH =
     &       +          r(117)  ! ORA1 + OH =
     &       +          r(118)  ! ORA2 + OH =
     &       + ( 0.650)*r(119)  ! PAA + OH =
     &       +          r(120)  ! PAN + OH =
     &       +          r(121)  ! PPN + OH =
     &       +          r(122)  ! MPAN + OH =
     &       +          r(123)  ! ONIT + OH =
     &       +          r(124)  ! NALD + OH =
     &       +          r(125)  ! ISON + OH =
     &       +          r(364)  ! ACT + OH =
     &       +          r(365)  ! ECH4 + OH =
     &       +          r(366)  ! DMS + OH =
     &       +          r(367)  ! DMS + OH + O2 =
      OH_loss = PA(nn)
c
c --- OH with CO
c
      nn = nn + 1
      ptname(nn)  = 'OHwCO'
      cpadesc(nn) = 'OH reaction rate with CO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 70)  ! CO + OH =
c
c --- OH with ECH4
c
      nn = nn + 1
      ptname(nn)  = 'OHwECH4'
      cpadesc(nn) = 'OH reaction rate with ECH4 (emitted CH4)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(365)  ! ECH4 + OH =
c
c --- OH with isoprene
c
      nn = nn + 1
      ptname(nn)  = 'OHwISOP'
      cpadesc(nn) = 'OH reaction rate with isoprene'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 86)  ! ISO + OH =
c
c --- OH with VOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwVOC'
      cpadesc(nn) = 'OH reaction rate with all VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 72)  ! ETH + OH =
     &       +          r( 73)  ! HC3 + OH =
     &       +          r( 74)  ! HC5 + OH =
     &       +          r( 75)  ! HC8 + OH =
     &       +          r( 76)  ! ETE + OH =
     &       +          r( 77)  ! OLT + OH =
     &       +          r( 78)  ! OLI + OH =
     &       +          r( 79)  ! DIEN + OH =
     &       +          r( 80)  ! ACE + OH =
     &       +          r( 81)  ! BEN + OH =
     &       +          r( 82)  ! TOL + OH =
     &       +          r( 83)  ! XYM + OH =
     &       +          r( 84)  ! XYP + OH =
     &       +          r( 85)  ! XYO + OH =
     &       +          r( 86)  ! ISO + OH =
     &       +          r( 87)  ! API + OH =
     &       +          r( 88)  ! LIM + OH =
     &       +          r( 89)  ! HCHO + OH =
     &       +          r( 90)  ! ACD + OH =
     &       +          r( 91)  ! ALD + OH =
     &       +          r( 92)  ! ACT + OH =
     &       +          r( 93)  ! MEK + OH =
     &       +          r( 94)  ! KET + OH =
     &       +          r(109)  ! MOH + OH =
     &       +          r(110)  ! EOH + OH =
     &       +          r(111)  ! ROH + OH =
     &       +          r(112)  ! ETEG + OH =
     &       +          r(364)  ! ACT + OH =
c
c --- OH with HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwHRVOC'
      cpadesc(nn) = 'OH reaction rate with HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 76)  ! ETE + OH =
     &       +          r( 77)  ! OLT + OH =
     &       +          r( 78)  ! OLI + OH =
     &       +          r( 79)  ! DIEN + OH =
c
c --- OH with Aromatics
c
      nn = nn + 1
      ptname(nn)  = 'OHwArom'
      cpadesc(nn) = 'OH reaction rate with aromatic VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 81)  ! BEN + OH =
     &       +          r( 82)  ! TOL + OH =
     &       +          r( 83)  ! XYM + OH =
     &       +          r( 84)  ! XYP + OH =
     &       +          r( 85)  ! XYO + OH =
c
c --- OH with Alkanes (except methane)
c
      nn = nn + 1
      ptname(nn)  = 'OHwAlkane'
      cpadesc(nn) = 'OH reaction rate with alkanes (except methane)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 72)  ! ETH + OH =
     &       +          r( 73)  ! HC3 + OH =
     &       +          r( 74)  ! HC5 + OH =
     &       +          r( 75)  ! HC8 + OH =
c
c --- Total HCHO production
c
      nn = nn + 1
      ptname(nn)  = 'HCHO_prod'
      cpadesc(nn) = 'Total HCHO production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.350)*r( 15)  ! UALD =
     &       +          r( 18)  ! HKET =
     &       + ( 0.670)*r( 19)  ! MACR =
     &       +          r( 22)  ! GLY =
     &       +          r( 28)  ! OP1 =
     &       +          r(109)  ! MOH + OH =
     &       + ( 0.350)*r(113)  ! OP1 + OH =
     &       + ( 0.350)*r(119)  ! PAA + OH =
     &       +          r(120)  ! PAN + OH =
     &       +          r(121)  ! PPN + OH =
     &       + ( 0.070)*r(125)  ! ISON + OH =
     &       +          r(126)  ! ETE + O3 =
     &       + ( 0.560)*r(127)  ! OLT + O3 =
     &       + ( 0.090)*r(128)  ! OLI + O3 =
     &       + ( 0.900)*r(129)  ! DIEN + O3 =
     &       + ( 0.580)*r(130)  ! ISO + O3 =
     &       + ( 0.040)*r(132)  ! LIM + O3 =
     &       + ( 0.100)*r(134)  ! MVK + O3 =
     &       + ( 0.080)*r(135)  ! UALD + O3 =
     &       + ( 0.050)*r(136)  ! DCB1 + O3 =
     &       + ( 0.050)*r(137)  ! DCB2 + O3 =
     &       + ( 0.680)*r(151)  ! MACR + NO3 =
     &       + ( 0.332)*r(152)  ! UALD + NO3 =
     &       +          r(172)  ! MO2 + NO =
     &       + ( 0.018)*r(175)  ! HC5P + NO =
     &       + ( 1.600)*r(177)  ! ETEP + NO =
     &       + ( 0.780)*r(178)  ! OLTP + NO =
     &       + ( 0.200)*r(188)  ! ISOP + NO =
     &       + ( 0.230)*r(189)  ! APIP + NO =
     &       + ( 0.430)*r(190)  ! LIMP + NO =
     &       +          r(193)  ! ACTP + NO =
     &       + ( 0.330)*r(194)  ! MEKP + NO =
     &       + ( 0.650)*r(196)  ! MACP + NO =
     &       + ( 0.500)*r(197)  ! MCP + NO =
     &       + ( 0.300)*r(198)  ! MVKP + NO =
     &       + ( 0.030)*r(199)  ! UALP + NO =
     &       + ( 0.287)*r(206)  ! OLND + NO =
     &       + ( 0.150)*r(233)  ! ACTP + HO2 =
     &       + ( 1.370)*r(248)  ! MO2 + MO2 =
     &       + ( 0.750)*r(249)  ! ETHP + MO2 =
     &       + ( 0.827)*r(250)  ! HC3P + MO2 =
     &       + ( 0.777)*r(251)  ! HC5P + MO2 =
     &       + ( 0.750)*r(252)  ! HC8P + MO2 =
     &       + ( 1.950)*r(253)  ! ETEP + MO2 =
     &       + ( 1.500)*r(254)  ! OLTP + MO2 =
     &       + ( 0.750)*r(255)  ! OLIP + MO2 =
     &       +          r(256)  ! BENP + MO2 =
     &       +          r(257)  ! TLP1 + MO2 =
     &       +          r(258)  ! TOLP + MO2 =
     &       +          r(259)  ! PER1 + MO2 =
     &       +          r(260)  ! XYL1 + MO2 =
     &       +          r(261)  ! XYLP + MO2 =
     &       +          r(262)  ! PER2 + MO2 =
     &       +          r(263)  ! XYOP + MO2 =
     &       + ( 1.310)*r(264)  ! ISOP + MO2 =
     &       + ( 0.750)*r(265)  ! APIP + MO2 =
     &       + ( 1.040)*r(266)  ! LIMP + MO2 =
     &       +          r(267)  ! ACO3 + MO2 =
     &       +          r(268)  ! RCO3 + MO2 =
     &       + ( 1.500)*r(269)  ! ACTP + MO2 =
     &       +          r(270)  ! MEKP + MO2 =
     &       + ( 0.750)*r(271)  ! KETP + MO2 =
     &       + ( 1.660)*r(272)  ! MACP + MO2 =
     &       + ( 1.500)*r(273)  ! MCP + MO2 =
     &       + ( 1.500)*r(274)  ! MVKP + MO2 =
     &       + ( 0.773)*r(275)  ! UALP + MO2 =
     &       +          r(276)  ! BALP + MO2 =
     &       +          r(277)  ! BAL1 + MO2 =
     &       +          r(278)  ! ADDC + MO2 =
     &       +          r(279)  ! MCTP + MO2 =
     &       +          r(280)  ! ORAP + MO2 =
     &       +          r(281)  ! OLNN + MO2 =
     &       + ( 0.965)*r(282)  ! OLND + MO2 =
     &       +          r(283)  ! ADCN + MO2 =
     &       +          r(284)  ! XO2 + MO2 =
     &       + ( 0.130)*r(286)  ! HC3P + ACO3 =
     &       + ( 0.042)*r(287)  ! HC5P + ACO3 =
     &       + ( 1.600)*r(289)  ! ETEP + ACO3 =
     &       +          r(290)  ! OLTP + ACO3 =
     &       + ( 0.750)*r(300)  ! ISOP + ACO3 =
     &       + ( 0.385)*r(302)  ! LIMP + ACO3 =
     &       +          r(305)  ! ACTP + ACO3 =
     &       + ( 0.330)*r(306)  ! MEKP + ACO3 =
     &       +          r(308)  ! MACP + ACO3 =
     &       +          r(309)  ! MCP + ACO3 =
     &       +          r(310)  ! MVKP + ACO3 =
     &       + ( 0.030)*r(311)  ! UALP + ACO3 =
     &       + ( 0.287)*r(318)  ! OLND + ACO3 =
     &       +          r(322)  ! MO2 + NO3 =
     &       + ( 0.024)*r(325)  ! HC5P + NO3 =
     &       + ( 1.600)*r(327)  ! ETEP + NO3 =
     &       + ( 0.790)*r(328)  ! OLTP + NO3 =
     &       + ( 0.750)*r(338)  ! ISOP + NO3 =
     &       + ( 0.385)*r(340)  ! LIMP + NO3 =
     &       +          r(343)  ! ACTP + NO3 =
     &       + ( 0.330)*r(344)  ! MEKP + NO3 =
     &       +          r(346)  ! MACP + NO3 =
     &       +          r(347)  ! MCP + NO3 =
     &       + ( 0.300)*r(348)  ! MVKP + NO3 =
     &       + ( 0.030)*r(349)  ! UALP + NO3 =
     &       + ( 0.287)*r(356)  ! OLND + NO3 =
     &       + ( 0.202)*r(359)  ! OLNN + OLND =
     &       + ( 0.504)*r(360)  ! OLND + OLND =
     &       +          r(366)  ! DMS + OH =
     &       +          r(368)  ! DMS + NO3 =
c
c --- HCHO production from HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_HRV'
      cpadesc(nn) = 'HCHO production rate from HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(126)  ! ETE + O3 =
     &       + ( 0.560)*r(127)  ! OLT + O3 =
     &       + ( 0.090)*r(128)  ! OLI + O3 =
     &       + ( 0.900)*r(129)  ! DIEN + O3 =
c
c --- HCHO production from Isoprene at first generation
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_ISP'
      cpadesc(nn) = 'HCHO production rate from isoprene at first generation'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.580)*r(130)  ! ISO + O3 =
     &       + ( 0.200)*r(188)  ! ISOP + NO =
     &       + ( 1.310)*r(264)  ! ISOP + MO2 =
     &       + ( 0.750)*r(300)  ! ISOP + ACO3 =
     &       + ( 0.750)*r(338)  ! ISOP + NO3 =
c
c --- Total HO2 production, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'HO2_prod'
      cpadesc(nn) = 'Total HO2 production rate, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 11)  ! HCHO =
     &       +          r( 12)  ! ACD =
     &       +          r( 13)  ! ALD =
     &       + ( 1.220)*r( 15)  ! UALD =
     &       +          r( 18)  ! HKET =
     &       + ( 0.660)*r( 19)  ! MACR =
     &       + ( 2.000)*r( 23)  ! GLY =
     &       +          r( 24)  ! MGLY =
     &       + ( 1.500)*r( 25)  ! DCB1 =
     &       + ( 1.500)*r( 26)  ! DCB2 =
     &       +          r( 27)  ! BALD =
     &       +          r( 28)  ! OP1 =
     &       +          r( 29)  ! OP2 =
     &       +          r( 31)  ! ONIT =
     &       +          r( 34)  ! O3 + OH =
     &       +          r( 43)  ! OH + H2 =
     &       +          r( 47)  ! H2O2 + OH =
     &       +          r( 58)  ! NO3 + OH =
     &       +          r( 69)  ! SO2 + OH =
     &       +          r( 70)  ! CO + OH =
     &       + ( 0.049)*r( 75)  ! HC8 + OH =
     &       + ( 0.350)*r( 80)  ! ACE + OH =
     &       + ( 0.648)*r( 81)  ! BEN + OH =
     &       + ( 0.177)*r( 82)  ! TOL + OH =
     &       + ( 0.177)*r( 83)  ! XYM + OH =
     &       + ( 0.177)*r( 84)  ! XYP + OH =
     &       + ( 0.177)*r( 85)  ! XYO + OH =
     &       +          r( 89)  ! HCHO + OH =
     &       +          r( 95)  ! HKET + OH =
     &       +          r( 99)  ! GLY + OH =
     &       + ( 0.520)*r(101)  ! DCB1 + OH =
     &       + ( 0.520)*r(102)  ! DCB2 + OH =
     &       + ( 0.560)*r(103)  ! DCB3 + OH =
     &       + ( 0.730)*r(105)  ! PHEN + OH =
     &       + ( 0.730)*r(106)  ! CSL + OH =
     &       +          r(107)  ! EPX + OH =
     &       +          r(109)  ! MOH + OH =
     &       +          r(110)  ! EOH + OH =
     &       +          r(111)  ! ROH + OH =
     &       +          r(112)  ! ETEG + OH =
     &       +          r(117)  ! ORA1 + OH =
     &       + ( 0.150)*r(126)  ! ETE + O3 =
     &       + ( 0.320)*r(127)  ! OLT + O3 =
     &       + ( 0.070)*r(128)  ! OLI + O3 =
     &       + ( 0.300)*r(129)  ! DIEN + O3 =
     &       + ( 0.250)*r(130)  ! ISO + O3 =
     &       + ( 0.100)*r(131)  ! API + O3 =
     &       + ( 0.100)*r(132)  ! LIM + O3 =
     &       + ( 0.140)*r(133)  ! MACR + O3 =
     &       + ( 0.110)*r(134)  ! MVK + O3 =
     &       + ( 0.072)*r(135)  ! UALD + O3 =
     &       +          r(136)  ! DCB1 + O3 =
     &       +          r(137)  ! DCB2 + O3 =
     &       +          r(138)  ! DCB3 + O3 =
     &       + ( 1.500)*r(139)  ! EPX + O3 =
     &       +          r(148)  ! HCHO + NO3 =
     &       +          r(152)  ! UALD + NO3 =
     &       +          r(153)  ! GLY + NO3 =
     &       + ( 1.500)*r(157)  ! EPX + NO3 =
     &       + ( 0.290)*r(160)  ! TR2 =
     &       + ( 0.010)*r(161)  ! TOLP =
     &       + ( 0.308)*r(162)  ! XY2 =
     &       + ( 0.010)*r(163)  ! XYLP =
     &       + ( 0.308)*r(164)  ! XYO2 =
     &       + ( 0.010)*r(165)  ! XYOP =
     &       +          r(172)  ! MO2 + NO =
     &       +          r(173)  ! ETHP + NO =
     &       + ( 0.660)*r(174)  ! HC3P + NO =
     &       + ( 0.200)*r(175)  ! HC5P + NO =
     &       + ( 0.606)*r(176)  ! HC8P + NO =
     &       +          r(177)  ! ETEP + NO =
     &       + ( 0.780)*r(178)  ! OLTP + NO =
     &       + ( 0.830)*r(179)  ! OLIP + NO =
     &       + ( 0.918)*r(180)  ! BENP + NO =
     &       + ( 0.950)*r(182)  ! TOLP + NO =
     &       + ( 0.500)*r(183)  ! PER1 + NO =
     &       + ( 0.950)*r(185)  ! XYLP + NO =
     &       + ( 0.950)*r(186)  ! PER2 + NO =
     &       + ( 0.950)*r(187)  ! XYOP + NO =
     &       + ( 0.880)*r(188)  ! ISOP + NO =
     &       + ( 0.820)*r(189)  ! APIP + NO =
     &       +          r(190)  ! LIMP + NO =
     &       + ( 0.670)*r(194)  ! MEKP + NO =
     &       + ( 0.770)*r(195)  ! KETP + NO =
     &       + ( 0.500)*r(197)  ! MCP + NO =
     &       + ( 0.300)*r(198)  ! MVKP + NO =
     &       +          r(199)  ! UALP + NO =
     &       +          r(202)  ! ADDC + NO =
     &       +          r(204)  ! ORAP + NO =
     &       +          r(205)  ! OLNN + NO =
     &       + ( 0.740)*r(248)  ! MO2 + MO2 =
     &       +          r(249)  ! ETHP + MO2 =
     &       + ( 0.894)*r(250)  ! HC3P + MO2 =
     &       + ( 0.842)*r(251)  ! HC5P + MO2 =
     &       + ( 0.910)*r(252)  ! HC8P + MO2 =
     &       +          r(253)  ! ETEP + MO2 =
     &       +          r(254)  ! OLTP + MO2 =
     &       +          r(255)  ! OLIP + MO2 =
     &       + ( 1.600)*r(256)  ! BENP + MO2 =
     &       +          r(257)  ! TLP1 + MO2 =
     &       + ( 2.000)*r(258)  ! TOLP + MO2 =
     &       + ( 2.000)*r(259)  ! PER1 + MO2 =
     &       +          r(260)  ! XYL1 + MO2 =
     &       + ( 2.000)*r(261)  ! XYLP + MO2 =
     &       + ( 2.000)*r(262)  ! PER2 + MO2 =
     &       + ( 2.000)*r(263)  ! XYOP + MO2 =
     &       +          r(264)  ! ISOP + MO2 =
     &       +          r(265)  ! APIP + MO2 =
     &       +          r(266)  ! LIMP + MO2 =
     &       + ( 0.900)*r(267)  ! ACO3 + MO2 =
     &       + ( 0.900)*r(268)  ! RCO3 + MO2 =
     &       + ( 0.500)*r(269)  ! ACTP + MO2 =
     &       + ( 0.834)*r(270)  ! MEKP + MO2 =
     &       +          r(271)  ! KETP + MO2 =
     &       + ( 0.500)*r(272)  ! MACP + MO2 =
     &       +          r(273)  ! MCP + MO2 =
     &       +          r(274)  ! MVKP + MO2 =
     &       +          r(275)  ! UALP + MO2 =
     &       +          r(276)  ! BALP + MO2 =
     &       +          r(277)  ! BAL1 + MO2 =
     &       + ( 2.000)*r(278)  ! ADDC + MO2 =
     &       +          r(279)  ! MCTP + MO2 =
     &       +          r(280)  ! ORAP + MO2 =
     &       + ( 2.000)*r(281)  ! OLNN + MO2 =
     &       + ( 0.500)*r(282)  ! OLND + MO2 =
     &       +          r(283)  ! ADCN + MO2 =
     &       +          r(284)  ! XO2 + MO2 =
     &       + ( 0.500)*r(285)  ! ETHP + ACO3 =
     &       + ( 0.394)*r(286)  ! HC3P + ACO3 =
     &       + ( 0.342)*r(287)  ! HC5P + ACO3 =
     &       + ( 0.303)*r(288)  ! HC8P + ACO3 =
     &       + ( 0.500)*r(289)  ! ETEP + ACO3 =
     &       + ( 0.500)*r(290)  ! OLTP + ACO3 =
     &       + ( 0.500)*r(291)  ! OLIP + ACO3 =
     &       + ( 0.600)*r(292)  ! BENP + ACO3 =
     &       +          r(294)  ! TOLP + ACO3 =
     &       +          r(295)  ! PER1 + ACO3 =
     &       +          r(297)  ! XYLP + ACO3 =
     &       +          r(298)  ! PER2 + ACO3 =
     &       +          r(299)  ! XYOP + ACO3 =
     &       + ( 0.500)*r(300)  ! ISOP + ACO3 =
     &       + ( 0.500)*r(301)  ! APIP + ACO3 =
     &       + ( 0.500)*r(302)  ! LIMP + ACO3 =
     &       + ( 0.330)*r(306)  ! MEKP + ACO3 =
     &       + ( 0.500)*r(307)  ! KETP + ACO3 =
     &       + ( 0.500)*r(309)  ! MCP + ACO3 =
     &       + ( 0.500)*r(310)  ! MVKP + ACO3 =
     &       + ( 0.500)*r(311)  ! UALP + ACO3 =
     &       + ( 2.000)*r(314)  ! ADDC + ACO3 =
     &       +          r(315)  ! MCTP + ACO3 =
     &       +          r(317)  ! OLNN + ACO3 =
     &       +          r(319)  ! ADCN + ACO3 =
     &       +          r(322)  ! MO2 + NO3 =
     &       +          r(323)  ! ETHP + NO3 =
     &       + ( 0.254)*r(324)  ! HC3P + NO3 =
     &       + ( 0.488)*r(325)  ! HC5P + NO3 =
     &       + ( 0.820)*r(326)  ! HC8P + NO3 =
     &       +          r(327)  ! ETEP + NO3 =
     &       + ( 0.790)*r(328)  ! OLTP + NO3 =
     &       + ( 0.860)*r(329)  ! OLIP + NO3 =
     &       +          r(330)  ! BENP + NO3 =
     &       +          r(332)  ! TOLP + NO3 =
     &       + ( 0.500)*r(333)  ! PER1 + NO3 =
     &       +          r(335)  ! XYLP + NO3 =
     &       +          r(336)  ! PER2 + NO3 =
     &       +          r(337)  ! XYOP + NO3 =
     &       +          r(338)  ! ISOP + NO3 =
     &       +          r(339)  ! APIP + NO3 =
     &       +          r(340)  ! LIMP + NO3 =
     &       + ( 0.670)*r(344)  ! MEKP + NO3 =
     &       +          r(345)  ! KETP + NO3 =
     &       +          r(347)  ! MCP + NO3 =
     &       + ( 0.300)*r(348)  ! MVKP + NO3 =
     &       +          r(349)  ! UALP + NO3 =
     &       +          r(352)  ! ADDC + NO3 =
     &       +          r(354)  ! ORAP + NO3 =
     &       +          r(355)  ! OLNN + NO3 =
     &       +          r(358)  ! OLNN + OLNN =
     &       + ( 0.500)*r(359)  ! OLNN + OLND =
c
c --- HO2 from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_O3'
      cpadesc(nn) = 'HO2 production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.150)*r(126)  ! ETE + O3 =
     &       + ( 0.320)*r(127)  ! OLT + O3 =
     &       + ( 0.070)*r(128)  ! OLI + O3 =
     &       + ( 0.300)*r(129)  ! DIEN + O3 =
     &       + ( 0.250)*r(130)  ! ISO + O3 =
     &       + ( 0.100)*r(131)  ! API + O3 =
     &       + ( 0.100)*r(132)  ! LIM + O3 =
     &       + ( 0.140)*r(133)  ! MACR + O3 =
     &       + ( 0.110)*r(134)  ! MVK + O3 =
     &       + ( 0.072)*r(135)  ! UALD + O3 =
     &       +          r(136)  ! DCB1 + O3 =
     &       +          r(137)  ! DCB2 + O3 =
     &       +          r(138)  ! DCB3 + O3 =
     &       + ( 1.500)*r(139)  ! EPX + O3 =
      newHO2_O3 = PA(nn)
c
c --- HO2 directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_pht'
      cpadesc(nn) = 'HO2 production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 11)  ! HCHO =
     &       +          r( 12)  ! ACD =
     &       +          r( 13)  ! ALD =
     &       + ( 1.220)*r( 15)  ! UALD =
     &       +          r( 18)  ! HKET =
     &       + ( 0.660)*r( 19)  ! MACR =
     &       + ( 2.000)*r( 23)  ! GLY =
     &       +          r( 24)  ! MGLY =
     &       + ( 1.500)*r( 25)  ! DCB1 =
     &       + ( 1.500)*r( 26)  ! DCB2 =
     &       +          r( 27)  ! BALD =
     &       +          r( 28)  ! OP1 =
     &       +          r( 29)  ! OP2 =
     &       +          r( 31)  ! ONIT =
      newHO2_pht = PA(nn)
c
c --- New HO2
c
      nn = nn + 1
      ptname(nn)  = 'HO2_new'
      cpadesc(nn) = 'New HO2 production rate, i.e. newHO2_O3 + newHO2_pht'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = newHO2_O3 + newHO2_pht
      OH_new = PA(nn)
c
c --- HO2 from HCHO photolysis
c
      nn = nn + 1
      ptname(nn)  = 'nwHO2_HCHO'
      cpadesc(nn) = 'HO2 production rate from HCHO photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 11)  ! HCHO =
c
c --- HO2 loss (except to HNO4) - computed above
c
      nn = nn + 1
      ptname(nn)  = 'HO2_loss'
      cpadesc(nn) = 'HO2 loss (except to HNO4)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = HO2_loss
c
c --- HO2wHO2 - computed above
c
      nn = nn + 1
      ptname(nn)  = 'HO2wHO2'
      cpadesc(nn) = 'HO2 + HO2 self-reaction'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = HO2wHO2
c
c --- HOx chain length = HOx reacted / new HOx
c
      if(.NOT. lcpacum) then
        nn = nn + 1
        ptname(nn)  = 'HOx_CL'
        cpadesc(nn) = 'HOx chain length = HOx reacted / new HOx'
        cpaunit(nn) = 'dimensionless'
        PA(nn) = ((OH_loss + HO2_loss) / max( 1.0E-12, (OH_new + HO2_new)))
     &           *sun*(dtfact/ppbfact)
      endif
c
c --- NO3 production
c
      nn = nn + 1
      ptname(nn)  = 'NO3_prod'
      cpadesc(nn) = 'NO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.200)*r(  9)  ! HNO4 =
     &       +          r( 33)  ! PAN =
     &       +          r( 37)  ! O3 + NO2 =
     &       +          r( 55)  ! NO2 + O =
     &       +          r( 57)  ! HNO3 + OH =
     &       +          r( 64)  ! N2O5 =
     &       +          r(120)  ! PAN + OH =
     &       +          r(121)  ! PPN + OH =
c
c --- NO3 production from N2O5
c
      nn = nn + 1
      ptname(nn)  = 'N2O5toNO3'
      cpadesc(nn) = 'NO3 production from N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 64)  ! N2O5 =
c
c --- NO3 loss
c
      nn = nn + 1
      ptname(nn)  = 'NO3_loss'
      cpadesc(nn) = 'NO3 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  5)  ! NO3 =
     &       +          r(  6)  ! NO3 =
     &       +          r( 58)  ! NO3 + OH =
     &       +          r( 59)  ! NO3 + HO2 =
     &       +          r( 60)  ! NO3 + NO =
     &       +          r( 61)  ! NO3 + NO2 =
     &       + ( 2.000)*r( 62)  ! NO3 + NO3 =
     &       +          r( 63)  ! NO3 + NO2 =
     &       +          r(141)  ! ETE + NO3 =
     &       +          r(142)  ! OLT + NO3 =
     &       +          r(143)  ! OLI + NO3 =
     &       +          r(144)  ! DIEN + NO3 =
     &       +          r(145)  ! ISO + NO3 =
     &       +          r(146)  ! API + NO3 =
     &       +          r(147)  ! LIM + NO3 =
     &       +          r(148)  ! HCHO + NO3 =
     &       +          r(149)  ! ACD + NO3 =
     &       +          r(150)  ! ALD + NO3 =
     &       +          r(151)  ! MACR + NO3 =
     &       +          r(152)  ! UALD + NO3 =
     &       +          r(153)  ! GLY + NO3 =
     &       +          r(154)  ! MGLY + NO3 =
     &       +          r(155)  ! PHEN + NO3 =
     &       +          r(156)  ! CSL + NO3 =
     &       +          r(157)  ! EPX + NO3 =
     &       +          r(158)  ! MCT + NO3 =
     &       +          r(159)  ! MPAN + NO3 =
     &       +          r(322)  ! MO2 + NO3 =
     &       +          r(323)  ! ETHP + NO3 =
     &       +          r(324)  ! HC3P + NO3 =
     &       +          r(325)  ! HC5P + NO3 =
     &       +          r(326)  ! HC8P + NO3 =
     &       +          r(327)  ! ETEP + NO3 =
     &       +          r(328)  ! OLTP + NO3 =
     &       +          r(329)  ! OLIP + NO3 =
     &       +          r(330)  ! BENP + NO3 =
     &       +          r(331)  ! TLP1 + NO3 =
     &       +          r(332)  ! TOLP + NO3 =
     &       +          r(333)  ! PER1 + NO3 =
     &       +          r(334)  ! XYL1 + NO3 =
     &       +          r(335)  ! XYLP + NO3 =
     &       +          r(336)  ! PER2 + NO3 =
     &       +          r(337)  ! XYOP + NO3 =
     &       +          r(338)  ! ISOP + NO3 =
     &       +          r(339)  ! APIP + NO3 =
     &       +          r(340)  ! LIMP + NO3 =
     &       +          r(341)  ! ACO3 + NO3 =
     &       +          r(342)  ! RCO3 + NO3 =
     &       +          r(343)  ! ACTP + NO3 =
     &       +          r(344)  ! MEKP + NO3 =
     &       +          r(345)  ! KETP + NO3 =
     &       +          r(346)  ! MACP + NO3 =
     &       +          r(347)  ! MCP + NO3 =
     &       +          r(348)  ! MVKP + NO3 =
     &       +          r(349)  ! UALP + NO3 =
     &       +          r(350)  ! BALP + NO3 =
     &       +          r(351)  ! BAL1 + NO3 =
     &       +          r(352)  ! ADDC + NO3 =
     &       +          r(353)  ! MCTP + NO3 =
     &       +          r(354)  ! ORAP + NO3 =
     &       +          r(355)  ! OLNN + NO3 =
     &       +          r(356)  ! OLND + NO3 =
     &       +          r(357)  ! ADCN + NO3 =
     &       +          r(361)  ! XO2 + NO3 =
     &       +          r(368)  ! DMS + NO3 =
c
c --- NO3 to N2O5
c
      nn = nn + 1
      ptname(nn)  = 'NO3toN2O5'
      cpadesc(nn) = 'NO3 conversion to N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 63)  ! NO3 + NO2 =
c
c --- Total production of ON compounds
c
      nn = nn + 1
      ptname(nn)  = 'ON_prod'
      cpadesc(nn) = 'Total production of organic nitrate (ON)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(152)  ! UALD + NO3 =
     &       + ( 0.065)*r(174)  ! HC3P + NO =
     &       + ( 0.136)*r(175)  ! HC5P + NO =
     &       + ( 0.261)*r(176)  ! HC8P + NO =
     &       + ( 0.030)*r(178)  ! OLTP + NO =
     &       + ( 0.050)*r(179)  ! OLIP + NO =
     &       + ( 0.082)*r(180)  ! BENP + NO =
     &       + ( 0.050)*r(182)  ! TOLP + NO =
     &       + ( 0.050)*r(183)  ! PER1 + NO =
     &       + ( 0.050)*r(185)  ! XYLP + NO =
     &       + ( 0.050)*r(186)  ! PER2 + NO =
     &       + ( 0.050)*r(187)  ! XYOP + NO =
     &       + ( 0.180)*r(189)  ! APIP + NO =
     &       +          r(205)  ! OLNN + NO =
     &       +          r(209)  ! BAL2 + NO2 =
     &       +          r(210)  ! CHO + NO2 =
     &       +          r(211)  ! MCTO + NO2 =
     &       +          r(244)  ! OLNN + HO2 =
     &       +          r(245)  ! OLND + HO2 =
     &       +          r(281)  ! OLNN + MO2 =
     &       + ( 0.500)*r(282)  ! OLND + MO2 =
     &       + ( 0.300)*r(283)  ! ADCN + MO2 =
     &       +          r(317)  ! OLNN + ACO3 =
     &       + ( 0.300)*r(319)  ! ADCN + ACO3 =
     &       +          r(355)  ! OLNN + NO3 =
     &       + ( 2.000)*r(358)  ! OLNN + OLNN =
     &       + ( 1.500)*r(359)  ! OLNN + OLND =
     &       +          r(360)  ! OLND + OLND =
     &       +          r(145)  ! ISO + NO3 =
     &       + ( 0.120)*r(188)  ! ISOP + NO =
      ON_prod = PA(nn)
c
c --- INTR production
c
      nn = nn + 1
      ptname(nn)  = 'INTR_prod'
      cpadesc(nn) = 'Total production of isoprene nitrates'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(145)  ! ISO + NO3 =
     &       + ( 0.120)*r(188)  ! ISOP + NO =
c
c --- ONIT production
c
      nn = nn + 1
      ptname(nn)  = 'ONIT_prod'
      cpadesc(nn) = 'Total production of the species ONIT'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(152)  ! UALD + NO3 =
     &       + ( 0.065)*r(174)  ! HC3P + NO =
     &       + ( 0.136)*r(175)  ! HC5P + NO =
     &       + ( 0.261)*r(176)  ! HC8P + NO =
     &       + ( 0.030)*r(178)  ! OLTP + NO =
     &       + ( 0.050)*r(179)  ! OLIP + NO =
     &       + ( 0.082)*r(180)  ! BENP + NO =
     &       + ( 0.050)*r(182)  ! TOLP + NO =
     &       + ( 0.050)*r(183)  ! PER1 + NO =
     &       + ( 0.050)*r(185)  ! XYLP + NO =
     &       + ( 0.050)*r(186)  ! PER2 + NO =
     &       + ( 0.050)*r(187)  ! XYOP + NO =
     &       + ( 0.180)*r(189)  ! APIP + NO =
     &       +          r(205)  ! OLNN + NO =
     &       +          r(209)  ! BAL2 + NO2 =
     &       +          r(210)  ! CHO + NO2 =
     &       +          r(211)  ! MCTO + NO2 =
     &       +          r(244)  ! OLNN + HO2 =
     &       +          r(245)  ! OLND + HO2 =
     &       +          r(281)  ! OLNN + MO2 =
     &       + ( 0.500)*r(282)  ! OLND + MO2 =
     &       + ( 0.300)*r(283)  ! ADCN + MO2 =
     &       +          r(317)  ! OLNN + ACO3 =
     &       + ( 0.300)*r(319)  ! ADCN + ACO3 =
     &       +          r(355)  ! OLNN + NO3 =
     &       + ( 2.000)*r(358)  ! OLNN + OLNN =
     &       + ( 1.500)*r(359)  ! OLNN + OLND =
     &       +          r(360)  ! OLND + OLND =
c
c --- HNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'HNO3_prod'
      cpadesc(nn) = 'HNO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 51)  ! NO + HO2 =
     &       +          r( 56)  ! NO2 + OH =
     &       + ( 0.300)*r( 59)  ! NO3 + HO2 =
     &       + ( 2.000)*r( 65)  ! N2O5 + H2O =
     &       +          r(148)  ! HCHO + NO3 =
     &       +          r(149)  ! ACD + NO3 =
     &       +          r(150)  ! ALD + NO3 =
     &       + ( 0.320)*r(151)  ! MACR + NO3 =
     &       +          r(153)  ! GLY + NO3 =
     &       +          r(154)  ! MGLY + NO3 =
     &       + ( 0.500)*r(155)  ! PHEN + NO3 =
     &       + ( 0.500)*r(156)  ! CSL + NO3 =
     &       + ( 0.500)*r(157)  ! EPX + NO3 =
     &       +          r(158)  ! MCT + NO3 =
     &       +          r(368)  ! DMS + NO3 =
     &       +          r(371)  ! ONIT =
     &       +          r(372)  ! ISON =
c
c --- HNO3 production from NO2 with OH - computed above
c
      nn = nn + 1
      ptname(nn)  = 'NO2wOH'
      cpadesc(nn) = 'HNO3 production from OH + NO2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = NO2wOH
c
c --- HNO3 production from N2O5 with H2O
c
      nn = nn + 1
      ptname(nn)  = 'N2O5wH2O'
      cpadesc(nn) = 'HNO3 production from N2O5 + H2O'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 65)  ! N2O5 + H2O =
      N2O5wH2O = PA(nn)
c
c --- NO2 or NO3 recycled from HNO3
c
      nn = nn + 1
      ptname(nn)  = 'NO2rcyHNO3'
      cpadesc(nn) = 'NO2 or NO3 recycled from HNO3 and organic nitrates'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  8)  ! HNO3 =
     &       +          r( 57)  ! HNO3 + OH =
c
c --- NO2 or NO3 recycled from organic nitrates
c
      nn = nn + 1
      ptname(nn)  = 'NO2rcyNTR'
      cpadesc(nn) = 'NO2 or NO3 recycled from organic nitrates'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 31)  ! ONIT =
     &       +          r(123)  ! ONIT + OH =
     &       +          r(124)  ! NALD + OH =
c
c --- Net production of PAN compounds
c
      nn = nn + 1
      ptname(nn)  = 'PAN_prdNet'
      cpadesc(nn) = 'Net production of PAN species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(166)  ! ACO3 + NO2 =
     &       +          r(168)  ! RCO3 + NO2 =
     &       +          r(170)  ! MACP + NO2 =
     &       -          r( 32)  ! PAN =
     &       -          r(122)  ! MPAN + OH =
     &       -          r(159)  ! MPAN + NO3 =
     &       -          r(167)  ! PAN =
     &       -          r(169)  ! PPN =
     &       -          r(171)  ! MPAN =
      PAN_prdNet = PA(nn)
c
c --- NO3 with VOC
c
      nn = nn + 1
      ptname(nn)  = 'NO3wVOC'
      cpadesc(nn) = 'NO3 reaction with VOC species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(141)  ! ETE + NO3 =
     &       +          r(142)  ! OLT + NO3 =
     &       +          r(143)  ! OLI + NO3 =
     &       +          r(144)  ! DIEN + NO3 =
     &       +          r(145)  ! ISO + NO3 =
     &       +          r(146)  ! API + NO3 =
     &       +          r(147)  ! LIM + NO3 =
     &       +          r(148)  ! HCHO + NO3 =
     &       +          r(149)  ! ACD + NO3 =
     &       +          r(150)  ! ALD + NO3 =
      NO3wVOC = PA(nn)
c
c --- LN/Q calculation
c
      if(.NOT. lcpacum) then
        nn = nn + 1
        ptname(nn)  = 'LNoQ'
        cpadesc(nn) = 'LN/Q metric for ozone sensitivity'
        cpaunit(nn) = 'dimensionless'
        PA(nn) = ( NO2wOH + NO3wVOC + N2O5wH2O +
     &             max( 0., PAN_prdNet ) + ON_prod )
     &           / max( 1.0E-12, (OH_new + HO2_new) )
        PA(nn) = max( 0., min( 10., PA(nn) ))*sun*(dtfact/ppbfact)
      endif
c
c --- Convert to ppb/hr
c
      Do n = 1,nn
        PA(n) = ppbfact*PA(n)
      End do
c
c --- Set num of outputs on first dummy call from pasetup
c
      npa_init = nn
c
      return
      end

