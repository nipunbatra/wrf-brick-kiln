c*** O3PRDSA
c
      subroutine o3prdsa(rrxn_irr,convfac,delo3,prdo3n,prdo3v,desto3)
      use filunit
      use chmstry
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c     This routine calculates the production/destruction of ozone
c     based on the chemical mechanism.  The production is in two
c     parts NOx-attributed and VOC-attributed.
c
c      Copyright 1996 - 2025
c     Ramboll
c
c   Argument descriptions:
c     Inputs:
c       rrxn_irr   R  array of reactions from last chemistry step
c       delo3      R  change in Ozone
c       convfac    R  conversion factor (Mmoles to ug/m^3)
c     Outputs:
c       prdo3n     R  ozone production attributed to NOx
c       prdo3v     R  ozone production attributed to VOC
c       desto3     R  ozone destruction
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     07/20/04   --gwilson--  Original development
c     02/21/05   --gyarwood-  Add SAPRC mechanism
c     10/06/05   --gyarwood-  Removed mechanism 2
c     01/08/06   --bkoo--     Added Mechanism 6 (CB05)
c     08/22/06   --gyarwood-  Revised ozone destruction calculation
c     06/24/16   --bkoo--     Updated Mechanism 4 (CB6r4)
c     11/28/16   --bkoo--     Revised desto3 using automated CPA
c     11/27/17   --bkoo--     Updated based on new CPAMECH routines (CMC v5.2.6)
c     02/21/20   --rlb--      Added Mechanism 1 (CB6r5)
c     09/22/20   --gy--       Added Mechanism 3 (CB6r5h)
c     03/24/21   --gy--       Added Mechanism 7 (CB7)
c     09/02/21   --gy--       Added Mechanism 9 (RACM2s21)
c     11/03/23   --cemery--   Removed CB05 (Mech 6)
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
      real rrxn_irr(*)
      real delo3
      real convfac
      real prdo3n
      real prdo3v
      real desto3
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      real prodo3, phno3, pho2h
      real HO2toNO2, HO2_loss
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- Estimate daytime ozone destruction due to:
c       (1) O1D + water 
c       (2) HO2 + O3
c       (3) OH + O3 --> HO2 + O2 (subtracting HO2 reacting with NO)
c       (4) O(3P) + VOC 
c       (5) O3 + VOC
c     Get the rates for OH + NO2 and HO2 + HO2
c
c     1. Copy from CPAMECH?.F
c     2. Replace 'PA(nn)' with 'desto3'
c     3. Replace 'r(' with 'rrxn_irr('
c
      if( idmech .EQ. 3 ) then
c
c  +  O1D + H2O
c
        desto3 =
     &       - rrxn_irr( 11)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
        desto3 = desto3
     &       - rrxn_irr( 13)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
        HO2toNO2 =
     &       +          rrxn_irr( 25)  ! HO2 + NO =
     &       +          rrxn_irr( 33)  ! NO3 + HO2 =
c
        HO2_loss =
     &       +          rrxn_irr( 13)  ! O3 + HO2 =
     &       +          rrxn_irr( 15)  ! HO2 + O =
     &       +          rrxn_irr( 18)  ! OH + HO2 =
     &       + ( 2.000)*rrxn_irr( 19)  ! HO2 + HO2 =
     &       + ( 2.000)*rrxn_irr( 20)  ! HO2 + HO2 + H2O =
     &       +          rrxn_irr( 25)  ! HO2 + NO =
     &       +          rrxn_irr( 33)  ! NO3 + HO2 =
     &       +          rrxn_irr( 57)  ! C2O3 + HO2 =
     &       +          rrxn_irr( 65)  ! CXO3 + HO2 =
     &       +          rrxn_irr( 72)  ! MEO2 + HO2 =
     &       +          rrxn_irr( 76)  ! XO2H + HO2 =
     &       +          rrxn_irr( 80)  ! XO2 + HO2 =
     &       +          rrxn_irr( 84)  ! XO2N + HO2 =
     &       +          rrxn_irr(100)  ! FORM + HO2 =
     &       + ( 0.800)*rrxn_irr(103)  ! HCO3 + HO2 =
     &       + ( 0.880)*rrxn_irr(145)  ! ISO2 + HO2 =
     &       + ( 0.175)*rrxn_irr(159)  ! EPX2 + HO2 =
     &       +          rrxn_irr(170)  ! BZO2 + HO2 =
     &       +          rrxn_irr(175)  ! TO2 + HO2 =
     &       +          rrxn_irr(179)  ! XLO2 + HO2 =
     &       +          rrxn_irr(185)  ! CRO + HO2 =
     &       +          rrxn_irr(202)  ! OPO3 + HO2 =
     &       +          rrxn_irr(214)  ! IO + HO2 =
     &       +          rrxn_irr(239)  ! CL + HO2 =
     &       +          rrxn_irr(243)  ! CLO + HO2 =
     &       +          rrxn_irr(288)  ! BR + HO2 =
     &       +          rrxn_irr(297)  ! BRO + HO2 =
     &       +          rrxn_irr(322)  ! I + HO2 =
c
        desto3 = desto3 - (
     &       + rrxn_irr( 12)  ! O3 + OH =
     &                     ) * (HO2_loss-HO2toNO2)/max( 1.0E-12, HO2_loss )
c
c  +  O3 + VOC
c
        desto3 = desto3
     &       - rrxn_irr(135)  ! ETH + O3 =
     &       - rrxn_irr(138)  ! OLE + O3 =
     &       - rrxn_irr(141)  ! IOLE + O3 =
     &       - rrxn_irr(149)  ! ISOP + O3 =
     &       - rrxn_irr(152)  ! ISPD + O3 =
     &       - rrxn_irr(165)  ! TERP + O3 =
c
c  +  O(3P) + VOC
c
!        desto3 = desto3
c
c  +  Halogen catalytic destruction (the limiting reactions)
c
        desto3 = desto3
     &       -          rrxn_irr(214)  ! IO + HO2 =
     &       - ( 0.400)*rrxn_irr(213)  ! IO + IO =
     &       -          rrxn_irr(242)  ! CLO + IO =
     &       -          rrxn_irr(296)  ! BRO + IO =
     &       -          rrxn_irr(216)  ! IO + NO2 =
     &       -          rrxn_irr(297)  ! BRO + HO2 =
     &       - ( 1.700)*rrxn_irr(294)  ! BRO + BRO =
     &       -          rrxn_irr(295)  ! BRO + CLO =
     &       -          rrxn_irr(296)  ! BRO + IO =
     &       -          rrxn_irr(299)  ! BRO + NO2 =
     &       -          rrxn_irr(243)  ! CLO + HO2 =
     &       - ( 1.400)*rrxn_irr(241)  ! CLO + CLO =
     &       -          rrxn_irr(242)  ! CLO + IO =
     &       -          rrxn_irr(295)  ! BRO + CLO =
c
        phno3 = rrxn_irr(45) + rrxn_irr(234)
        pho2h = rrxn_irr(19) + rrxn_irr(20)
c
      elseif( idmech .EQ. 1 .OR. idmech .EQ. 4) then
c
c --- CB6r5 Mechanism 1 and CB6r4 Mechanism 4 ---
c
c  +  O1D + H2O
c
        desto3 =
     &         - rrxn_irr( 11)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
        desto3 = desto3
     &         - rrxn_irr( 13)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
        HO2toNO2 = 
     &         +          rrxn_irr( 25)  ! HO2 + NO =
     &         +          rrxn_irr( 33)  ! NO3 + HO2 =
c
        HO2_loss = 
     &         +          rrxn_irr( 13)  ! O3 + HO2 =
     &         +          rrxn_irr( 15)  ! HO2 + O =
     &         +          rrxn_irr( 18)  ! OH + HO2 =
     &         + ( 2.000)*rrxn_irr( 19)  ! HO2 + HO2 =
     &         + ( 2.000)*rrxn_irr( 20)  ! HO2 + HO2 + H2O =
     &         +          rrxn_irr( 25)  ! HO2 + NO =
     &         +          rrxn_irr( 33)  ! NO3 + HO2 =
     &         +          rrxn_irr( 57)  ! C2O3 + HO2 =
     &         +          rrxn_irr( 65)  ! CXO3 + HO2 =
     &         +          rrxn_irr( 72)  ! MEO2 + HO2 =
     &         +          rrxn_irr( 76)  ! XO2H + HO2 =
     &         +          rrxn_irr( 80)  ! XO2 + HO2 =
     &         +          rrxn_irr( 84)  ! XO2N + HO2 =
     &         +          rrxn_irr(100)  ! FORM + HO2 =
     &         + ( 0.800)*rrxn_irr(103)  ! HCO3 + HO2 =
     &         + ( 0.880)*rrxn_irr(145)  ! ISO2 + HO2 =
     &         + ( 0.175)*rrxn_irr(159)  ! EPX2 + HO2 =
     &         +          rrxn_irr(170)  ! BZO2 + HO2 =
     &         +          rrxn_irr(175)  ! TO2 + HO2 =
     &         +          rrxn_irr(179)  ! XLO2 + HO2 =
     &         +          rrxn_irr(185)  ! CRO + HO2 =
     &         +          rrxn_irr(202)  ! OPO3 + HO2 =
     &         +          rrxn_irr(214)  ! IO + HO2 =
c
        desto3 = desto3 - (
     &         + rrxn_irr( 12)  ! O3 + OH =
     &                       ) * (HO2_loss-HO2toNO2)/HO2_loss
c
c  +  O3 + VOC
c
        desto3 = desto3
     &         - rrxn_irr(135)  ! ETH + O3 =
     &         - rrxn_irr(138)  ! OLE + O3 =
     &         - rrxn_irr(141)  ! IOLE + O3 =
     &         - rrxn_irr(149)  ! ISOP + O3 =
     &         - rrxn_irr(152)  ! ISPD + O3 =
     &         - rrxn_irr(165)  ! TERP + O3 =
c
c  +  O(3P) + VOC
c
!        desto3 = desto3
c
c  +  Halogen catalytic destruction (the limiting reactions)
c
        desto3 = desto3
     &         -          rrxn_irr(214)  ! IO + HO2 =
     &         - ( 0.400)*rrxn_irr(213)  ! IO + IO =
     &         -          rrxn_irr(216)  ! IO + NO2 =
c
        phno3 = rrxn_irr(45) 
        pho2h = rrxn_irr(19) + rrxn_irr(20) 
c
      elseif( idmech .EQ. 5 ) then
c
c --- SAPRC07TC Mechanism 5 ---
c
c  +  O1D + H2O
c
        desto3 =
     &         - rrxn_irr( 20)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
        desto3 = desto3
     &         - rrxn_irr( 36)  ! HO2 + O3 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
        HO2toNO2 = 
     &         +          rrxn_irr( 31)  ! HO2 + NO =
     &         + ( 0.800)*rrxn_irr( 39)  ! NO3 + HO2 =
c
        HO2_loss = 
     &         +          rrxn_irr( 31)  ! HO2 + NO =
     &         +          rrxn_irr( 36)  ! HO2 + O3 =
     &         + ( 2.000)*rrxn_irr( 37)  ! HO2 + HO2 =
     &         + ( 2.000)*rrxn_irr( 38)  ! HO2 + HO2 + H2O =
     &         +          rrxn_irr( 39)  ! NO3 + HO2 =
     &         +          rrxn_irr( 43)  ! OH + HO2 =
     &         +          rrxn_irr( 47)  ! MEO2 + HO2 =
     &         +          rrxn_irr( 48)  ! MEO2 + HO2 =
     &         +          rrxn_irr( 53)  ! RO2C + HO2 =
     &         +          rrxn_irr( 58)  ! RO2X + HO2 =
     &         +          rrxn_irr( 67)  ! MCO3 + HO2 =
     &         +          rrxn_irr( 77)  ! RCO3 + HO2 =
     &         +          rrxn_irr( 88)  ! BZC3 + HO2 =
     &         +          rrxn_irr(100)  ! MAC3 + HO2 =
     &         +          rrxn_irr(112)  ! BZO + HO2 =
c
        desto3 = desto3 - (
     &         + rrxn_irr( 30)  ! OH + O3 =
     &                       ) * (HO2_loss-HO2toNO2)/HO2_loss
c
c  +  O3 + VOC
c
        desto3 = desto3
     &         - rrxn_irr(247)  ! AFG1 + O3 =
     &         - rrxn_irr(250)  ! AFG2 + O3 =
     &         - rrxn_irr(253)  ! AFG3 + O3 =
     &         - rrxn_irr(255)  ! MACR + O3 =
     &         - rrxn_irr(260)  ! MVK + O3 =
     &         - rrxn_irr(264)  ! IPRD + O3 =
     &         - rrxn_irr(275)  ! ACRO + O3 =
     &         - rrxn_irr(515)  ! ETHE + O3 =
     &         - rrxn_irr(519)  ! PRPE + O3 =
     &         - rrxn_irr(523)  ! BD13 + O3 =
     &         - rrxn_irr(527)  ! ISOP + O3 =
     &         - rrxn_irr(531)  ! APIN + O3 =
     &         - rrxn_irr(535)  ! ACYE + O3 =
     &         - rrxn_irr(549)  ! OLE1 + O3 =
     &         - rrxn_irr(553)  ! OLE2 + O3 =
     &         - rrxn_irr(559)  ! TERP + O3 =
     &         - rrxn_irr(563)  ! SESQ + O3 =
c
c  +  O(3P) + VOC
c
        desto3 = desto3
     &         - rrxn_irr(257)  ! MACR + O3P =
     &         - rrxn_irr(261)  ! MVK + O3P =
     &         - rrxn_irr(277)  ! ACRO + O3P =
     &         - rrxn_irr(517)  ! ETHE + O3P =
     &         - rrxn_irr(521)  ! PRPE + O3P =
     &         - rrxn_irr(525)  ! BD13 + O3P =
     &         - rrxn_irr(529)  ! ISOP + O3P =
     &         - rrxn_irr(533)  ! APIN + O3P =
     &         - rrxn_irr(551)  ! OLE1 + O3P =
     &         - rrxn_irr(555)  ! OLE2 + O3P =
     &         - rrxn_irr(561)  ! TERP + O3P =
c
        phno3 = rrxn_irr(25)
        pho2h = rrxn_irr(37) + rrxn_irr(38)
c
      elseif( idmech .EQ. 7) then
c
c --- CB7 Mechanism 7 ---
c
c  +  O1D + H2O
c
        desto3 =
     &         - rrxn_irr( 11)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
        desto3 = desto3
     &         - rrxn_irr( 13)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
      HO2toNO2 =
     &       +          rrxn_irr( 25)  ! NO + HO2 =
     &       +          rrxn_irr( 32)  ! NO3 + HO2 =
c
      HO2_loss =
     &       +          rrxn_irr( 13)  ! O3 + HO2 =
     &       +          rrxn_irr( 15)  ! HO2 + O =
     &       +          rrxn_irr( 18)  ! OH + HO2 =
     &       + ( 2.000)*rrxn_irr( 19)  ! HO2 + HO2 =
     &       + ( 2.000)*rrxn_irr( 20)  ! HO2 + HO2 + H2O =
     &       +          rrxn_irr( 25)  ! NO + HO2 =
     &       +          rrxn_irr( 32)  ! NO3 + HO2 =
     &       +          rrxn_irr( 60)  ! C2O3 + HO2 =
     &       + ( 0.500)*rrxn_irr( 67)  ! CXO3 + HO2 =
     &       +          rrxn_irr( 73)  ! OPO3 + HO2 =
     &       +          rrxn_irr( 79)  ! MEO2 + HO2 =
     &       +          rrxn_irr( 85)  ! XO2H + HO2 =
     &       +          rrxn_irr( 88)  ! XO2 + HO2 =
     &       +          rrxn_irr( 91)  ! XO2N + HO2 =
     &       +          rrxn_irr(150)  ! BZO2 + HO2 =
     &       +          rrxn_irr(154)  ! TO2 + HO2 =
     &       +          rrxn_irr(158)  ! XLO2 + HO2 =
     &       +          rrxn_irr(171)  ! CRO + HO2 =
     &       + ( 0.940)*rrxn_irr(179)  ! ISO2 + HO2 =
     &       + ( 0.520)*rrxn_irr(197)  ! APO2 + HO2 =
     &       + ( 0.930)*rrxn_irr(203)  ! TPO2 + HO2 =
     &       + ( 0.900)*rrxn_irr(209)  ! SQO2 + HO2 =
     &       +          rrxn_irr(218)  ! IO + HO2 =
c
        desto3 = desto3 - (
     &         + rrxn_irr( 12)  ! O3 + OH =
     &                       ) * (HO2_loss-HO2toNO2)/HO2_loss
c
c  +  O3 + VOC
c
        desto3 = desto3
     &       - rrxn_irr(140)  ! ETH + O3 =
     &       - rrxn_irr(143)  ! OLE + O3 =
     &       - rrxn_irr(146)  ! IOLE + O3 =
     &       - rrxn_irr(162)  ! OPEN + O3 =
     &       - rrxn_irr(166)  ! XOPN + O3 =
     &       - rrxn_irr(182)  ! ISOP + O3 =
     &       - rrxn_irr(199)  ! APIN + O3 =
     &       - rrxn_irr(205)  ! TERP + O3 =
     &       - rrxn_irr(211)  ! SQT + O3 =
c
c  +  Halogen catalytic destruction (the limiting reactions)
c
        desto3 = desto3
     &       -          rrxn_irr(218)  ! IO + HO2 =
     &       - ( 0.400)*rrxn_irr(217)  ! IO + IO =
     &       -          rrxn_irr(220)  ! IO + NO2 =
c
        phno3 = rrxn_irr(41) + rrxn_irr(42)
        pho2h = rrxn_irr(19) + rrxn_irr(20)
c
      elseif( idmech .EQ. 9) then
c
c --- RACM2s21 Mechanism 9 ---
c
c  +  O1D + H2O
c
        desto3 =
     &         - rrxn_irr( 42)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
        desto3 = desto3
     &         - rrxn_irr( 35)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
      HO2toNO2 =
     &       +          rrxn_irr( 50)  ! NO + HO2 =
     &       + ( 0.700)*rrxn_irr( 59)  ! NO3 + HO2 =

c
      HO2_loss =
     &       +          rrxn_irr( 35)  ! O3 + HO2 =
     &       +          rrxn_irr( 44)  ! OH + HO2 =
     &       + ( 2.000)*rrxn_irr( 45)  ! HO2 + HO2 =
     &       + ( 2.000)*rrxn_irr( 46)  ! HO2 + HO2 + H2O =
     &       +          rrxn_irr( 50)  ! NO + HO2 =
     &       +          rrxn_irr( 51)  ! NO + HO2 =
     &       +          rrxn_irr( 59)  ! NO3 + HO2 =
     &       +          rrxn_irr(212)  ! MO2 + HO2 =
     &       +          rrxn_irr(213)  ! ETHP + HO2 =
     &       +          rrxn_irr(214)  ! HC3P + HO2 =
     &       +          rrxn_irr(215)  ! HC5P + HO2 =
     &       +          rrxn_irr(216)  ! HC8P + HO2 =
     &       +          rrxn_irr(217)  ! ETEP + HO2 =
     &       +          rrxn_irr(218)  ! OLTP + HO2 =
     &       +          rrxn_irr(219)  ! OLIP + HO2 =
     &       +          rrxn_irr(220)  ! BENP + HO2 =
     &       +          rrxn_irr(221)  ! TLP1 + HO2 =
     &       +          rrxn_irr(222)  ! TOLP + HO2 =
     &       +          rrxn_irr(223)  ! PER1 + HO2 =
     &       +          rrxn_irr(224)  ! XYL1 + HO2 =
     &       +          rrxn_irr(225)  ! XYLP + HO2 =
     &       +          rrxn_irr(226)  ! PER2 + HO2 =
     &       +          rrxn_irr(227)  ! XYOP + HO2 =
     &       +          rrxn_irr(228)  ! ISOP + HO2 =
     &       +          rrxn_irr(229)  ! APIP + HO2 =
     &       +          rrxn_irr(230)  ! LIMP + HO2 =
     &       +          rrxn_irr(231)  ! ACO3 + HO2 =
     &       +          rrxn_irr(232)  ! RCO3 + HO2 =
     &       +          rrxn_irr(233)  ! ACTP + HO2 =
     &       +          rrxn_irr(234)  ! MEKP + HO2 =
     &       +          rrxn_irr(235)  ! KETP + HO2 =
     &       +          rrxn_irr(236)  ! MACP + HO2 =
     &       +          rrxn_irr(237)  ! MCP + HO2 =
     &       +          rrxn_irr(238)  ! MVKP + HO2 =
     &       +          rrxn_irr(239)  ! UALP + HO2 =
     &       +          rrxn_irr(240)  ! ADDC + HO2 =
     &       +          rrxn_irr(241)  ! CHO + HO2 =
     &       +          rrxn_irr(242)  ! MCTP + HO2 =
     &       +          rrxn_irr(243)  ! ORAP + HO2 =
     &       +          rrxn_irr(244)  ! OLNN + HO2 =
     &       +          rrxn_irr(245)  ! OLND + HO2 =
     &       +          rrxn_irr(246)  ! ADCN + HO2 =
     &       +          rrxn_irr(247)  ! XO2 + HO2 =
c
        desto3 = desto3 - (
     &         + rrxn_irr( 12)  ! O3 + OH =
     &                       ) * (HO2_loss-HO2toNO2)/HO2_loss
c
c  +  O3 + VOC
c
        desto3 = desto3
     &       - rrxn_irr(126)  ! ETE + O3 =
     &       - rrxn_irr(127)  ! OLT + O3 =
     &       - rrxn_irr(128)  ! OLI + O3 =
     &       - rrxn_irr(129)  ! DIEN + O3 =
     &       - rrxn_irr(130)  ! ISO + O3 =
     &       - rrxn_irr(131)  ! API + O3 =
     &       - rrxn_irr(132)  ! LIM + O3 =
     &       - rrxn_irr(133)  ! MACR + O3 =
     &       - rrxn_irr(134)  ! MVK + O3 =
     &       - rrxn_irr(135)  ! UALD + O3 =
     &       - rrxn_irr(136)  ! DCB1 + O3 =
     &       - rrxn_irr(137)  ! DCB2 + O3 =
     &       - rrxn_irr(138)  ! DCB3 + O3 =
     &       - rrxn_irr(139)  ! EPX + O3 =
     &       - rrxn_irr(140)  ! MCTO + O3 =
c
        phno3 = rrxn_irr(56)
        pho2h = rrxn_irr(45) + rrxn_irr(46)
c
      else
        write(iout,'(//,a)') 'ERROR in O3PRDSA:'
        write(iout,'(/,1X,A,I10)')
     &             'Unknown chemical mechanism ID number: ',idmech
        write(iout,'(/,1X,2A)')
     &             'Ozone source apportionment is not available ',
     &             'for this chemical mechanism'
        call camxerr()
      endif
c
c --- Assign ozone production to VOC or NOx ---
c
      desto3 = AMIN1(desto3*convfac,delo3)
      prodo3 = delo3 - desto3
      if( pho2h/phno3 .LE. 0.35) then
        prdo3n = 0.0
        prdo3v = prodo3
      else
        prdo3n = prodo3
        prdo3v = 0.0
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
