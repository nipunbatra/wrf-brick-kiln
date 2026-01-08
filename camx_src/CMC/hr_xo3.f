      subroutine hr_xo3(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.32 250801
c
c     HR_XO3 solves halogen oxides by EBI iterations
c
c     Copyright 1996 - 2025
c     Ramboll
c     Created by the CMC version 5.3.2
c
c     Mechanism is CB6r5h plus GLY/MGLY SOA
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        EBISOLV
c
c --- Argument definitions:
c        y0  - initial Conc for this step (ppm)
c        yh  - current Conc iteration (ppm)
c        y1  - next Conc iteration (ppm)
c        H2O - water vapor Conc (ppm)
c        M   - total gas Conc (ppm)
c        O2  - oxygen Conc (ppm)
c        CH4 - methane Conc (ppm)
c        H2  - hydrogen Conc (ppm)
c        ny  - dimension of y0, y1 and yh
c        rk  - rate constants (ppm-n hr-1)
c        r   - reaction rates (hr-1)
c        nr  - dimension of rk and r
c        dt  - time step (hr)
c
c --- Includes:
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
c --- Arguments:
      integer ny, nr
      real y0(ny+1), y1(ny+1), yh(ny+1)
      real rk(nr), r(nr)
      real H2O, M, O2, CH4, H2, dt
c
c --- Local variables:
      real N2
      integer iter, nit
      real TINY
      parameter (TINY = 1.0E-20)
      real gainI, lossI
      real gainIO, lossIO
      real gainOIO, lossOIO
      real gainCL, lossCL
      real gainCLO, lossCLO
      real gainBR, lossBR
      real gainBRO, lossBRO
      real r209, r210, r211, r212, r213, r214, r215, r216, r217, r218
      real r219, r220, r221, r223, r235, r236, r237, r238, r239, r240
      real r241, r242, r243, r244, r245, r246, r247, r248, r249, r250
      real r253, r254, r255, r256, r257, r258, r259, r260, r261, r262
      real r263, r265, r268, r269, r270, r271, r272, r273, r274, r275
      real r276, r277, r278, r279, r280, r281, r282, r283, r284, r285
      real r286, r287, r288, r289, r290, r291, r292, r293, r294, r295
      real r296, r297, r298, r299, r300, r301, r303, r304, r305, r306
      real r307, r308, r309, r310, r312, r313, r314, r315, r316, r317
      real r318, r319, r320, r321, r322, r323, r324, r325, r327
c
c --- Entry Point
c
      N2 = M - O2
      r209 = r(209)
      r210 = r(210)
      r221 = r(221)
      r223 = r(223)
      r235 = r(235)
      r236 = r(236)
      r237 = r(237)
      r247 = r(247)
      r248 = r(248)
      r249 = r(249)
      r250 = r(250)
      r265 = r(265)
      r283 = r(283)
      r284 = r(284)
      r285 = r(285)
      r286 = r(286)
      r291 = r(291)
      r292 = r(292)
      r300 = r(300)
      r301 = r(301)
      r312 = r(312)
      r313 = r(313)
      r314 = r(314)
      r315 = r(315)
      r316 = r(316)
      r317 = r(317)
      r318 = r(318)
      r319 = r(319)
      r320 = r(320)
      r321 = r(321)
      r323 = r(323)
      r325 = r(325)
c
      nit = 2
c
      do iter = 1,nit
c
c --- Update reaction rates that change with iteration of y1
c
      r211 = rk(211)*y1(lI)*yh(lO3)
      r212 = rk(212)*y1(lIO)
      r213 = rk(213)*y1(lIO)*y1(lIO)
      r214 = rk(214)*y1(lIO)*yh(lHO2)
      r215 = rk(215)*y1(lIO)*yh(lNO)
      r216 = rk(216)*y1(lIO)*yh(lNO2)
      r217 = rk(217)*y1(lOIO)
      r218 = rk(218)*y1(lOIO)*yh(lOH)
      r219 = rk(219)*y1(lOIO)*y1(lIO)
      r220 = rk(220)*y1(lOIO)*yh(lNO)
      r238 = rk(238)*y1(lCL)*yh(lO3)
      r239 = rk(239)*y1(lCL)*yh(lHO2)
      r240 = rk(240)*y1(lCL)*H2
      r241 = rk(241)*y1(lCLO)*y1(lCLO)
      r242 = rk(242)*y1(lCLO)*y1(lIO)
      r243 = rk(243)*y1(lCLO)*yh(lHO2)
      r244 = rk(244)*y1(lCLO)*yh(lMEO2)
      r245 = rk(245)*y1(lCLO)*yh(lNO)
      r246 = rk(246)*y1(lCLO)*yh(lNO2)
      r253 = rk(253)*yh(lFORM)*y1(lCL)
      r254 = rk(254)*yh(lALD2)*y1(lCL)
      r255 = rk(255)*yh(lALDX)*y1(lCL)
      r256 = rk(256)*yh(lGLY)*y1(lCL)
      r257 = rk(257)*yh(lGLYD)*y1(lCL)
      r258 = rk(258)*yh(lMGLY)*y1(lCL)
      r259 = rk(259)*yh(lACET)*y1(lCL)
      r260 = rk(260)*yh(lKET)*y1(lCL)
      r261 = rk(261)*yh(lMEOH)*y1(lCL)
      r262 = rk(262)*yh(lETOH)*y1(lCL)
      r263 = rk(263)*yh(lISPD)*y1(lCL)
      r268 = rk(268)*CH4*y1(lCL)
      r269 = rk(269)*yh(lECH4)*y1(lCL)
      r270 = rk(270)*yh(lETHA)*y1(lCL)
      r271 = rk(271)*yh(lPRPA)*y1(lCL)
      r272 = rk(272)*yh(lPAR)*y1(lCL)
      r273 = rk(273)*yh(lETHY)*y1(lCL)
      r274 = rk(274)*yh(lETH)*y1(lCL)
      r275 = rk(275)*yh(lOLE)*y1(lCL)
      r276 = rk(276)*yh(lIOLE)*y1(lCL)
      r277 = rk(277)*yh(lISOP)*y1(lCL)
      r278 = rk(278)*yh(lTERP)*y1(lCL)
      r279 = rk(279)*yh(lTOL)*y1(lCL)
      r280 = rk(280)*yh(lXYL)*y1(lCL)
      r281 = rk(281)*yh(lCRES)*y1(lCL)
      r282 = rk(282)*yh(lDMS)*y1(lCL)
      r287 = rk(287)*y1(lBR)*yh(lO3)
      r288 = rk(288)*y1(lBR)*yh(lHO2)
      r289 = rk(289)*y1(lBR)*yh(lNO2)
      r290 = rk(290)*y1(lBR)*yh(lNO3)
      r293 = rk(293)*y1(lBRO)
      r294 = rk(294)*y1(lBRO)*y1(lBRO)
      r295 = rk(295)*y1(lBRO)*y1(lCLO)
      r296 = rk(296)*y1(lBRO)*y1(lIO)
      r297 = rk(297)*y1(lBRO)*yh(lHO2)
      r298 = rk(298)*y1(lBRO)*yh(lNO)
      r299 = rk(299)*y1(lBRO)*yh(lNO2)
      r303 = rk(303)*y1(lBR)*yh(lFORM)
      r304 = rk(304)*y1(lBR)*yh(lALD2)
      r305 = rk(305)*y1(lBR)*yh(lALDX)
      r306 = rk(306)*y1(lBR)*yh(lETH)
      r307 = rk(307)*y1(lBR)*yh(lOLE)
      r308 = rk(308)*y1(lBR)*yh(lIOLE)
      r309 = rk(309)*y1(lBR)*yh(lISOP)
      r310 = rk(310)*y1(lBR)*yh(lTERP)
      r322 = rk(322)*y1(lI)*yh(lHO2)
      r324 = rk(324)*y1(lI)*yh(lNO2)
      r327 = rk(327)*y1(lBR)*yh(lBRN2)
c
c --- Gain and loss terms for the XO species
c
      gainI = 0.0
     &       + ( 2.000)*r209
     &       +          r210
     &       +          r212
     &       + ( 0.400)*r213
     &       +          r215
     &       +          r217
     &       +          r221
     &       +          r223
     &       +          r236
     &       +          r242
     &       +          r284
     &       +          r296
     &       +          r312
     &       + ( 2.000)*r313
     &       +          r314
     &       +          r315
     &       +          r323
     &       +          r325
c
      lossI = 0.0
     &       +          r211
     &       +          r322
     &       +          r324
c
      gainIO = 0.0
     &       +          r211
     &       +          r220
c
      lossIO = 0.0
     &       +          r212
     &       + ( 2.000)*r213
     &       +          r214
     &       +          r215
     &       +          r216
     &       +          r219
     &       +          r242
     &       +          r296
c
      gainOIO = 0.0
     &       + ( 0.400)*r213
     &       +          r221
c
      lossOIO = 0.0
     &       +          r217
     &       +          r218
     &       +          r219
     &       +          r220
c
      gainCL = 0.0
     &       + ( 2.000)*r235
     &       +          r236
     &       +          r237
     &       + ( 1.400)*r241
     &       +          r242
     &       +          r245
     &       +          r249
     &       +          r250
     &       +          r265
     &       + ( 0.210)*r273
     &       +          r285
     &       +          r295
     &       +          r315
c
      lossCL = 0.0
     &       +          r238
     &       +          r239
     &       +          r240
     &       +          r253
     &       +          r254
     &       +          r255
     &       +          r256
     &       +          r257
     &       +          r258
     &       +          r259
     &       +          r260
     &       +          r261
     &       +          r262
     &       +          r263
     &       +          r268
     &       +          r269
     &       +          r270
     &       +          r271
     &       +          r272
     &       +          r273
     &       +          r274
     &       +          r275
     &       +          r276
     &       +          r277
     &       +          r278
     &       +          r279
     &       +          r280
     &       +          r281
     &       +          r282
c
      gainCLO = 0.0
     &       +          r238
     &       + ( 0.220)*r239
     &       +          r247
     &       +          r248
c
      lossCLO = 0.0
     &       + ( 2.000)*r241
     &       +          r242
     &       +          r243
     &       +          r244
     &       +          r245
     &       +          r246
     &       +          r295
c
      gainBR = 0.0
     &       + ( 2.000)*r283
     &       +          r284
     &       +          r285
     &       +          r286
     &       +          r291
     &       +          r292
     &       +          r293
     &       + ( 1.700)*r294
     &       +          r295
     &       +          r296
     &       +          r298
     &       +          r300
     &       + ( 0.850)*r301
     &       +          r314
     &       + ( 3.000)*r316
     &       + ( 3.000)*r317
     &       + ( 2.000)*r318
     &       +          r319
     &       +          r320
     &       +          r321
c
      lossBR = 0.0
     &       +          r287
     &       +          r288
     &       +          r289
     &       +          r290
     &       +          r303
     &       +          r304
     &       +          r305
     &       +          r306
     &       +          r307
     &       +          r308
     &       +          r309
     &       +          r310
     &       +          r327
c
      gainBRO = 0.0
     &       +          r287
     &       +          r290
     &       + ( 0.150)*r301
c
      lossBRO = 0.0
     &       +          r293
     &       + ( 2.000)*r294
     &       +          r295
     &       +          r296
     &       +          r297
     &       +          r298
     &       +          r299
c
c --- EBI solution for the XO species
c
      y1(lI) = MIN(1.0, MAX(TINY,(y0(lI) + gainI*dt)
     &                         / (1.0 + lossI*dt/y1(lI))))
c
      y1(lIO) = MIN(1.0, MAX(TINY,(y0(lIO) + gainIO*dt)
     &                         / (1.0 + lossIO*dt/y1(lIO))))
c
      y1(lOIO) = MIN(1.0, MAX(TINY,(y0(lOIO) + gainOIO*dt)
     &                         / (1.0 + lossOIO*dt/y1(lOIO))))
c
      y1(lCL) = MIN(1.0, MAX(TINY,(y0(lCL) + gainCL*dt)
     &                         / (1.0 + lossCL*dt/y1(lCL))))
c
      y1(lCLO) = MIN(1.0, MAX(TINY,(y0(lCLO) + gainCLO*dt)
     &                         / (1.0 + lossCLO*dt/y1(lCLO))))
c
      y1(lBR) = MIN(1.0, MAX(TINY,(y0(lBR) + gainBR*dt)
     &                         / (1.0 + lossBR*dt/y1(lBR))))
c
      y1(lBRO) = MIN(1.0, MAX(TINY,(y0(lBRO) + gainBRO*dt)
     &                         / (1.0 + lossBRO*dt/y1(lBRO))))
c
      enddo
c
      return
      end

