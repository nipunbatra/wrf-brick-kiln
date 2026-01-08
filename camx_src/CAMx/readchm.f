      subroutine readchm
      use filunit
      use grid
      use chmstry
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c     READCHM reads the CAMx chemistry parameter file, which defines
c     the chemical system to be simulated
c
c      Copyright 1996 - 2025
c     Ramboll
c  
c     Modifications:
c        10/26/99  Added check for zero Henry's law constants to avoid a
c                  divide by 0 later in the code
c        10/20/00  Added CAMx version as first record on chemistry parameters
c                  file
c        1/9/02    Aerosol size cut points and density now defined on
c                  chemistry parameters file; removed conversion of aerosol
c                  BDNL values from ug/m3 to umol/m3
c        1/18/02   Added EC, PFIN, and PCRS to mechanism 4 species list
c        12/12/02  Expanded species list for Mechanism 4
c        3/26/03   Added surface resistance scaling factor to gas params
c        4/21/04   Incorporated sectional PM (Mech 4 CMU)
c        10/14/04  Added Mechanism 10 (user defined)
c        10/05/05  Split gas-phase mechanism ID and aerosol option string
c                  to separate records in the chemistry parameters file
c        12/29/06 -bkoo-     Expanded species list for the updated SOA scheme
c        01/08/07 -bkoo-     Added Mechanism 6 (CB05)
c                            Revised the code to re-order fast species (now use FASTS list)
c        01/10/07 -bkoo-     PM modules now linked to mechanisms 4 thru 6
c        07/04/07 -bkoo-     Added code to set pointer to hydrolysis of N2O5
c        12/15/08 -gwilson-  Added code to handle averaging of radicals
c        01/29/09 -bkoo-     De-activated mechanism 1
c        04/22/10 -gwilson-  Removed the code that eliminted PM species for
c                            SAPRC
c        11/20/10 -gwilson-  Consolidated radical and concentration arrays
c        12/21/10 -bkoo-     Added Mechanism 7 (CB6)
c        01/20/11 -gwilson-  De-activated Mechanisms 3 and 4
c        03/29/11 -cemery-   Revised to allow for inert PM with gas-phase
c                            chemistry and to support in-line TUV with
c                            aerosol optical depth
c        05/20/11 -bkoo-     Added Mechanisms 8 and 9 (CB05 and SAPRC99 extensions)
c        04/12/12 -cemery-   Added T and P adjustments to photolysis rates
c        07/20/12 -ou-       Added Mechanism 1 (CB6 with Iodine chemistry)
c        10/08/12 -jjung-    Added Mechanism 2 (CB6r1)
c        12/12/12 -bkoo-     Set I,IO,OIO as IEH steady-state radicals (mech1,8,9)
c        01/15/13 -cemery-   Replaced gas "diffrat" with "Molwt"; diffrat now
c                            calculated from Molwt when read
c        04/12/13 -cemery-   Added surface model species/reactions
c        06/28/13 -bkoo-     Added mapping for new GREASD PiG chemistry
c        10/14/13 -bkoo-     Added code to set pointer to hydrolysis of organic nitrate
c        10/24/13 -bkoo-     Revised Mechanism 2 (CB6r1 -> CB6r2)
c        03/18/14 -bkoo-     Added species for benzene SOA
c        06/27/14 -bkoo-     Added Mechanism 3 (CB6r2h)
c        08/25/14 -cemery-   Added snow effects to surface model
c        11/11/14 -bkoo-     Added Mechanism 4 (CB6r3)
c        12/07/14 -bkoo-     Added VBS option
c        01/08/16 -bkoo-     Updated for SAPRC99->SAPRC07TC
c                                        Revised VBS
c                                        Removal of IEH solver
c        02/29/16 -bkoo-     Added pointer to heterogeneous rxn of N2O5 + HCL
c        06/24/16 -bkoo-     Updated Mechanism 4 (CB6r4)
c        07/20/16 -bkoo-     Added pointer to heterogeneous hydrolysis of INTR
c        08/25/16 -bkoo-     Updated for new SOAP
c        09/02/16 -bkoo-     Added pointer to het rxn of SO2; revised het rxn pointers
c        11/28/16 -bkoo-     Revised VBS precursors
c        10/18/17 -cemery-   Restructured chemparam file header for explicit
c                            choice of PM size and chemistry options (removed
c                            original SOAP)
c        11/22/17 -bkoo-     Added PFE/PMN/PK/PCA/PMG
c        01/12/18 -bkoo-     Removed BNZA/TOLA/XYLA/ISP/TRP
c        07/23/18 -bkoo-     Check NH3 if Bi-Di NH3 is on
c        01/06/19 -cemery-   Added PAL/PSI/PTI
c        01/09/19 -cemery-   Added DMS to CB6r4 (remains Mech 4)
c        01/30/20 -cemery-   Removed CB6r2 (Mech 2)
c        02/20/20 -rlb-      Added CB6r5 (Mech 1)
c        11/17/20 -cemery-   Moved definition of TUV photolysis file labels here
c        03/24/21 -gy-       Added CB7 (Mech 7)
c        04/15/21 -gwilson-  Now puts molwt into global array to be used later
c        08/10/21 -gy-       Added RACM2s21 (Mech 9) abbreviated as RACM2 or RC2
c        08/10/21 -gy-       Removed check ngas < mchgas because unhelpful
c        09/08/21 -cemery-   Enforce 8-character max on species names
c                            (to accomodate deposition prefixes)
c        01/20/22 -gwilson-  Added "blank" species names for traditional tracers
c        04/14/22 -gy-       Update CB7 to CB7r1 (Mech 7)
c        10/02/22 -gy-       Replace CB6r2h with CB6r5h (Mech 2)
c        03/31/23 -gy-       Fix bug molwt_spc not re-ordered in species sort
c        11/03/23 -cemery-   Removed CB05 (Mech 6)
c        11/03/23 -cemery-   Removed CMU PM option (only CB05 was plugged into CMU)
c        04/05/24 -cemery-   SOAP3: rename species, POA->HOA (unevaporated POA)
c                                   add SVOA (evaporated POA)
c                                   updated yields
c        05/16/24 -cemery-   All gas mechs include heterogeneous (M)GLY->CGLY for
c                            SOAP BOA0
c                            Differentiate TERP+NO3 reaction from TERP+OH/O3
c                            Add SOA yields from APIN
c
c     Input arguments: 
c        none 
c 
c     Output arguments: 
c        none 
c            
c     Routines called: 
c        ALLOC_CHMSTRY
c        AEROSET
c        EXPTBL
c        KTHERM
c            
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'flags.inc'
      include 'ncf_iodat.inc'
      include 'ddmchm.inc'
      include 'lsbox.inc'
      include 'soap.inc'
      include 'vbs.inc'
      include 'soap3.inc'
c
      integer istrln
c
      integer, parameter :: NRADNM = 65
c
      character*180    record, recpht1
      character*10     nametmp, splist(NSPNAM), radlist(NRADNM,10)
      character*10     blank,camxv,camxvin
      character*10     gasopt,siaopt
      character*20     soaopt
      character*10     tmpnam, tmpnam1
      integer          mchactv(10), mchgas(10), mchaero(10), mchrad(10)
      integer          mchrxn(10), mchphot(10)
c
      integer          jno2s07, jno2cb6, jno2cb7, jno2rc2
      integer          jo3s07, jo3cb6, jo3cb7, jo3rc2
      integer          jhcho1s07, jhcho1cb6, jhcho1cb7, 
     &                 jhcho1rc2
      integer          jhcho2s07, jhcho2cb6, jhcho2cb7, 
     &                 jhcho2rc2
      integer          jch3chos07, jch3chocb6, jch3chocb7, 
     &                 jch3chorc2
      integer          ipigs07, ipigcb6, ipigcb7, ipigrc2
      integer, dimension(NRPIGCHEM) :: lrmaps07, lrmapcb6r4,
     &                                 lrmapcb7, lrmaprc2
      integer          ihet01_s07, ihet01_cb6, ihet01_cb7,
     &                 ihet01_rc2
      integer          ihet02_cb6r5h
      integer          ihet03_cb6r4, ihet03_cb7, ihet03_rc2
      integer          ihet04_cb6r4, ihet04_rc2
      integer          ihet05_cb6r4, ihet05_cb7, ihet05_rc2
      integer          ihet06_cb6r4, ihet06_cb6r5, ihet06_cb6r5h,
     &                 ihet06_cb7, ihet06_s07, ihet06_rc2
      integer          ihet07_cb6r4, ihet07_cb6r5, ihet07_cb6r5h,
     &                 ihet07_cb7, ihet07_s07, ihet07_rc2
      integer, dimension(NSOAPRXN) :: iksoap_cb6, iksoap_cb7, iksoap_rc2,
     &                                iksoap_s07
c
      integer          npar(7)
      integer          nsec_c,  nphot, nerr, nhit, iaero, ierr, iref
      integer          ncms, nppm, num, iunits, isec
      integer          i, j, l, nn, n, m, nl
      integer          ibin
      integer          rhadj_tmp
      real             tdum(3), pdum(3), tmp
      real             bdnl_tmp,roprt_tmp,bext_tmp,ssa_tmp,molwt_spc
      real*8           dsec_i(MXSECT+1)
c
      double precision rxnpar(MXREACT,12)
      real             kdum(MXREACT,3)
      integer          rxntyp(MXREACT)
      integer          rxnord(MXREACT)
c
c-----Data that define the mechanism/solver options:
c
      data camxv  /'VERSION7.3'/
      data blank  /'BLANK     '/
c
      data (splist(j),j=1,270)
     &            /'APO2      ','AUTX      ','BZO2      ',
     &             'C2O3      ','CRO       ','CXO3      ',
     &             'DDO2      ','EPX2      ','HO2       ',
     &             'I         ','IO        ','ISO2      ',
     &             'MEO2      ','O         ','O1D       ',
     &             'OH        ','OIO       ','OPO3      ',
     &             'RO2       ','ROR       ','TO2       ',
     &             'TPO2      ','XLO2      ','XO2       ',
     &             'XO2H      ','XO2N      ','XPAR      ',
     &             'XPRP      ','BR        ','BRO       ',
     &             'CL        ','CLO       ','HCO3      ',
     &             'ARPX      ','DDCA      ','EDOH      ',
     &             'ETHR      ','HACT      ','HKET      ',
     &             'IBTA      ','IPOH      ','NPOH      ',
     &             'PDOH      ','ROH       ','SXD5      ',
     &             'TPRD      ','CGLY      ','AACD      ',
     &             'ACET      ','ALD2      ','ALDX      ',
     &             'BENZ      ','BRCL      ','CAT1      ',
     &             'CLAD      ','CLAO      ','CO        ',
     &             'CRES      ','CRON      ','DMS       ',
     &             'ECH4      ','EPOX      ','ETH       ',
     &             'ETHA      ','ETHY      ','ETOH      ',
     &             'FACD      ','FORM      ','GLY       ',
     &             'GLYD      ','H2O2      ','HCL       ',
     &             'HI        ','HIO3      ','HNO3      ',
     &             'HOI       ','HONO      ','HPLD      ',
     &             'I2        ','I2O2      ','IBR       ',
     &             'ICL       ','INO2      ','INO3      ',
     &             'INTR      ','IOLE      ','ISOP      ',
     &             'ISPD      ','ISPX      ','IXOY      ',
     &             'KET       ','MEOH      ','MEPX      ',
     &             'MGLY      ','N2O5      ','NH3       ',
     &             'NO        ','NO2       ','NO3       ',
     &             'NTR1      ','NTR2      ','O3        ',
     &             'OLE       ','OPAN      ','OPEN      ',
     &             'PACD      ','PAN       ','PANX      ',
     &             'PAR       ','PNA       ','PRPA      ',
     &             'ROOH      ','SO2       ','SQT       ',
     &             'SULF      ','TERP      ','TOL       ',
     &             'XOPN      ','XYL       ','CCRS      ',
     &             'ACG1      ','ACG2      ','BCG1      ',
     &             'BCG2      ','CPRM      ','FCRS      ',
     &             'FPRM      ','IVOA      ','SVOA      ',
     &             'NA        ','PAL       ','PCA       ',
     &             'PCL       ','PEC       ','PFE       ',
     &             'PH2O      ','PK        ','PMG       ',
     &             'PMN       ','PNH4      ','PNO3      ',
     &             'HOA       ','PSI       ','PSO4      ',
     &             'PTI       ','AOA1      ','AOA2      ',
     &             'BOA1      ','BOA2      ','AOA0      ',
     &             'BOA0      ','BR2       ','BRN2      ',
     &             'BRN3      ','CH3I      ','CL2       ',
     &             'CLN2      ','CLN3      ','FMBR      ',
     &             'FMCL      ','HBR       ','HOBR      ',
     &             'HOCL      ','IALK      ','MB2       ',
     &             'MB2C      ','MB3       ','MBC       ',
     &             'MBC2      ','MI2       ','MIB       ',
     &             'MIC       ','SSBR      ','SSCL      ',
     &             'SSN3      ','HG0       ','HG2       ',
     &             'HGIIP     ','HGIIPC    ','HGP       ',
     &             'IVOB      ','IVOD      ','IVOG      ',
     &             'PAP0      ','PAP1      ','PAP2      ',
     &             'PAP3      ','PAP4      ','PAS0      ',
     &             'PAS1      ','PAS2      ','PAS3      ',
     &             'PAS4      ','PBS0      ','PBS1      ',
     &             'PBS2      ','PBS3      ','PBS4      ',
     &             'PCP0      ','PCP1      ','PCP2      ',
     &             'PCP3      ','PCP4      ','PFP0      ',
     &             'PFP1      ','PFP2      ','PFP3      ',
     &             'PFP4      ','VAP1      ','VAP2      ',
     &             'VAP3      ','VAP4      ','VAS1      ',
     &             'VAS2      ','VAS3      ','VAS4      ',
     &             'VBS1      ','VBS2      ','VBS3      ',
     &             'VBS4      ','VCP1      ','VCP2      ',
     &             'VCP3      ','VCP4      ','VFP1      ',
     &             'VFP2      ','VFP3      ','VFP4      ',
     &             'ACO3      ','ACTP      ','ADCN      ',
     &             'ADDC      ','APIP      ','BAL1      ',
     &             'BAL2      ','BALP      ','BENP      ',
     &             'CHO       ','ETEP      ','ETHP      ',
     &             'HC3P      ','HC5P      ','HC8P      ',
     &             'KETP      ','LIMP      ','MACP      ',
     &             'MCP       ','MCTO      ','MCTP      ',
     &             'MEKP      ','MO2       ','MVKP      ',
     &             'OLIP      ','OLND      ','OLNN      ',
     &             'OLTP      ','ORAP      ','PER1      ',
     &             'PER2      ','TLP1      ','TOLP      ',
     &             'TR2       ','UALP      ','XY2       ',
     &             'XYL1      ','XYLP      ','XYO2      ',
     &             'XYOP      ','ACD       ','ACE       '/
      data (splist(j),j=271,NSPNAM)
     &            /'ACT       ','ALD       ','API       ',
     &             'BEN       ','CSL       ','DCB1      ',
     &             'DCB2      ','DCB3      ','DIEN      ',
     &             'EOH       ','EPX       ','ETE       ',
     &             'ETEG      ','HC3       ','HC5       ',
     &             'HC8       ','HNO4      ','ISHP      ',
     &             'ISO       ','ISON      ','LIM       ',
     &             'MAHP      ','MCT       ','MOH       ',
     &             'NALD      ','OLI       ','OLT       ',
     &             'ONIT      ','OP1       ','OP2       ',
     &             'ORA1      ','ORA2      ','PAA       ',
     &             'PHEN      ','PPN       ','UALD      ',
     &             'XYM       ','XYO       ','XYP       ',
     &             'BZC3      ','BZO       ','MAC3      ',
     &             'MCO3      ','O3P       ','RCO3      ',
     &             'RO2C      ','RO2X      ','TBUO      ',
     &             'XACE      ','XACR      ','XAF1      ',
     &             'XAF2      ','XAF3      ','XBAC      ',
     &             'XBAL      ','XCCH      ','XCO       ',
     &             'XGLD      ','XGLY      ','XHCH      ',
     &             'XHO2      ','XIPR      ','XMA3      ',
     &             'XMAC      ','XMC3      ','XMEK      ',
     &             'XMEO      ','XMGL      ','XMVK      ',
     &             'XNO2      ','XOH       ','XPD2      ',
     &             'XRC3      ','XRCH      ','XRN3      ',
     &             'XTBU      ','Y6PX      ','YAPX      ',
     &             'YRPX      ','ZRN3      ','ACRO      ',
     &             'ACYE      ','AFG1      ','AFG2      ',
     &             'AFG3      ','ALK1      ','ALK2      ',
     &             'ALK3      ','ALK4      ','ALK5      ',
     &             'APIN      ','ARO1      ','ARO2      ',
     &             'B124      ','BACL      ','BALD      ',
     &             'BD13      ','CCHO      ','CO3H      ',
     &             'COOH      ','ETHE      ','HCHO      ',
     &             'IPRD      ','MACR      ','MEK       ',
     &             'MPAN      ','MVK       ','MXYL      ',
     &             'NPHE      ','OLE1      ','OLE2      ',
     &             'OXYL      ','PAN2      ','PBZN      ',
     &             'PRD2      ','PRPE      ','PXYL      ',
     &             'R6PX      ','RAPX      ','RCHO      ',
     &             'RNO3      ','RO3H      ','SESQ      ',
     &             'TOLU      ','XN        ','TR01      ',
     &             'TR02      ','TR03      ','TR04      ',
     &             'TR05      ','TR06      ','TR07      ',
     &             'TR08      ','TR09      '/
c
      data (radlist(j,1),j=1,26)
     &            /'O1D       ','O         ','I         ',
     &             'IO        ','OIO       ','OH        ',
     &             'HO2       ','C2O3      ','XO2       ',
     &             'XO2N      ','CXO3      ','MEO2      ',
     &             'TO2       ','ROR       ','HCO3      ',
     &             'CRO       ','BZO2      ','EPX2      ',
     &             'ISO2      ','OPO3      ','RO2       ',
     &             'XLO2      ','XO2H      ','XPRP      ',
     &             'XPAR      ','CGLY      '/
      data (radlist(j,3),j=1,28)
     &            /'O1D       ','O         ','CL        ',
     &             'CLO       ','BR        ','BRO       ',
     &             'I         ','IO        ','OIO       ',
     &             'OH        ','HO2       ','C2O3      ',
     &             'XO2       ','XO2N      ','CXO3      ',
     &             'MEO2      ','TO2       ','ROR       ',
     &             'HCO3      ','CRO       ','BZO2      ',
     &             'EPX2      ','ISO2      ','OPO3      ',
     &             'RO2       ','XLO2      ','XO2H      ',
     &             'CGLY      '/
      data (radlist(j,4),j=1,26)
     &            /'O1D       ','O         ','I         ',
     &             'IO        ','OIO       ','OH        ',
     &             'HO2       ','C2O3      ','XO2       ',
     &             'XO2N      ','CXO3      ','MEO2      ',
     &             'TO2       ','ROR       ','HCO3      ',
     &             'CRO       ','BZO2      ','EPX2      ',
     &             'ISO2      ','OPO3      ','RO2       ',
     &             'XLO2      ','XO2H      ','XPRP      ',
     &             'XPAR      ','CGLY      '/
      data (radlist(j,5),j=1,46)
     &            /'O1D       ','O3P       ','OH        ',
     &             'HO2       ','MCO3      ','RCO3      ',
     &             'MEO2      ','TBUO      ','BZO       ',
     &             'BZC3      ','MAC3      ','RO2C      ',
     &             'RO2X      ','XHO2      ','XOH       ',
     &             'XACE      ','XACR      ','XAF1      ',
     &             'XAF2      ','XAF3      ','XBAC      ',
     &             'XBAL      ','XCCH      ','XCO       ',
     &             'XGLD      ','XGLY      ','XHCH      ',
     &             'XIPR      ','XMA3      ','XMAC      ',
     &             'XMC3      ','XMEK      ','XMEO      ',
     &             'XMGL      ','XMVK      ','XNO2      ',
     &             'XPD2      ','XRC3      ','XRCH      ',
     &             'XRN3      ','XTBU      ','YRPX      ',
     &             'Y6PX      ','YAPX      ','ZRN3      ',
     &             'CGLY      '/
      data (radlist(j,7),j=1,27)
     &            /'O1D       ', 'O         ', 'OH        ',
     &             'HO2       ', 'C2O3      ', 'CXO3      ',
     &             'OPO3      ', 'MEO2      ', 'XO2       ',
     &             'XO2H      ', 'XO2N      ', 'ISO2      ',
     &             'APO2      ', 'TPO2      ', 'CGLY      ',
     &             'EPX2      ', 'BZO2      ', 'TO2       ',
     &             'XLO2      ', 'XPRP      ', 'RO2       ',
     &             'XPAR      ', 'CRO       ', 'ROR       ',
     &             'I         ', 'IO        ', 'OIO       '/
      data (radlist(j,9),j=1,48)
     &            /'O1D       ', 'O         ', 'ADCN      ',
     &             'ADDC      ', 'APIP      ', 'BAL1      ',
     &             'BAL2      ', 'BALP      ', 'BENP      ',
     &             'CHO       ', 'ETEP      ', 'ETHP      ',
     &             'HC3P      ', 'HC5P      ', 'HC8P      ',
     &             'HO2       ', 'ISOP      ', 'KETP      ',
     &             'LIMP      ', 'MACP      ', 'MCP       ',
     &             'MCTO      ', 'MCTP      ', 'MEKP      ',
     &             'MO2       ', 'MVKP      ', 'ACO3      ',
     &             'ACTP      ', 'OH        ', 'OLIP      ',
     &             'OLND      ', 'OLNN      ', 'OLTP      ', 
     &             'ORAP      ', 'PER1      ', 'PER2      ',
     &             'RCO3      ', 'TLP1      ', 'TOLP      ',
     &             'TR2       ', 'UALP      ', 'XO2       ',
     &             'XY2       ', 'XYL1      ', 'XYLP      ',
     &             'XYO2      ', 'XYOP      ', 'CGLY      '/
c
      data (soap2_spec_names(j),j=1,NUM_SOAP_SPECIES)
     &            /'POA       ',
     &             'SOPA      ',
     &             'SOPB      ',
     &             'SOA2      ',
     &             'SOA1      ',
     &             'SOA4      ',
     &             'SOA3      ',
     &             'CG1       ',
     &             'CG2       ',
     &             'CG3       ',
     &             'CG4       '/
      data (soap3_spec_names(j),j=1,NUM_SOAP_SPECIES)
     &            /'HOA       ',
     &             'AOA0      ',
     &             'BOA0      ',
     &             'AOA1      ',
     &             'AOA2      ',
     &             'BOA1      ',
     &             'BOA2      ',
     &             'ACG2      ',
     &             'ACG1      ',
     &             'BCG2      ',
     &             'BCG1      '/
      data soap_emiss_names(IDX_SOAP_IVOC)   /'IVOC'/
      data soap_emiss_names(IDX_SOAP_IVOD)   /'IVOD'/
      data soap_emiss_names(IDX_SOAP_IVOA)   /'IVOA'/
      data soap_emiss_names(IDX_SOAP_IVOB)   /'IVOB'/
      data soap_emiss_names(IDX_SOAP_POA_OP) /'POA_OP'/
      data soap_emiss_names(IDX_SOAP_POA_GV) /'POA_GV'/
      data soap_emiss_names(IDX_SOAP_POA_DV) /'POA_DV'/
      data soap_emiss_names(IDX_SOAP_POA_MC) /'POA_MC'/
      data soap_emiss_names(IDX_SOAP_POA_IC) /'POA_IC'/
      data soap_emiss_names(IDX_SOAP_POA_AV) /'POA_AV'/
      data soap_emiss_names(IDX_SOAP_POA_BB) /'POA_BB'/
c
c-----NOTE: mchgas is max among gas mech, VBS + 2 optional Hg specs
c           mchaero is max among PM mech, VBS + 3 optional Hg specs
c     data mchactv  /  1,  0,  1,  1,  1,  0,  1,  0,  1,  0 /
c     data mchgas   /116,  0,128,116,126,  0,101,  0,125,999 /
c     data mchaero  / 47,  0, 47, 47, 47,  0, 47,  0, 47,999 /
      data mchrad   / 25,  0, 27, 25, 45,  0, 26,  0, 47,  0 /
      data mchrxn   /239,  0,332,238,567,  0,232,  0,374,  0 /
      data mchphot  / 34,  0, 54, 34, 41,  0, 34,  0, 34,  0 /
c
c-----Point to the NO + NO + O2 reaction within main mechanisms
c
      data ipigs07,ipigcb6,ipigcb7,ipigrc2
     &                              / 10, 24, 24, 52 /
c
c-----Pointers to heterogeneous reactions within main mechanisms
c
      data ihet01_s07,ihet01_cb6,ihet01_cb7,ihet01_rc2
     &                                  / 13, 39, 37, 65 / ! N2O5 + H2O -> 2 HNO3
      data ihet02_cb6r5h                           / 237 / ! N2O5 + HCL -> CLN2 + HNO3
      data ihet03_cb6r4,ihet03_cb7,ihet03_rc2
     &                                    / 207, 97, 371 / ! NTR2 -> HNO3
      data ihet04_cb6r4,ihet04_rc2              / 229, 0 / ! INTR -> HNO3
      data ihet05_cb6r4,ihet05_cb7,ihet05_rc2
     &                                    / 230, 52, 369 / ! SO2 -> SULF
      data ihet06_cb6r4,ihet06_cb6r5,ihet06_cb6r5h,ihet06_CB7,ihet06_s07,ihet06_rc2
     &                   / 234, 235, 328, 231, 566, 373 /  ! GLY -> CGLY
      data ihet07_cb6r4,ihet07_cb6r5,ihet07_cb6r5h,ihet07_CB7,ihet07_s07,ihet07_rc2
     &                   / 235, 236, 329, 232, 567, 374 /  ! MGLY -> CGLY
c
c-----Pointers to photolysis reactions with pressure adjustments
c
      data jno2s07,jno2cb6,jno2cb7,jno2rc2
     &                                  /  1, 1, 1, 1/
      data jo3s07,jo3cb6,jo3cb7,jo3rc2
     &                                  / 18, 9, 9, 3/
      data jhcho1s07,jhcho1cb6,jhcho1cb7,jhcho1rc2
     &                              /204, 97, 101, 11/
      data jhcho2s07,jhcho2cb6,jhcho2cb7,jhcho2rc2
     &                              /205, 98, 102, 10/
      data jch3chos07,jch3chocb6,jch3chocb7,jch3chorc2
     &                             /209, 106, 106, 12/
c
c-----Cross-reference PiG reactions to main mechanisms
c
      data lrmaps07  /  1,  2,  7, 10,  0, 18, 21, 20, 23, 25,
     &                  8, 17, 16,  9, 15, 11, 12, 13, 44, 29,
     &                204,205, 31 /
      data lrmapcb6r4 / 1,  2,  3, 24, 41,  9, 10, 11, 43, 45,
     &                 26, 27, 28, 29, 30, 36, 37, 39, 52,120,
     &                 97, 98, 25 /
      data lrmapcb7   / 1,  2,  3, 24,  0,  9, 10, 11, 39, 41,
     &                 26, 27, 28, 29, 30, 34, 35, 37, 51, 50,
     &                 101, 102, 25 /
      data lrmaprc2   / 1, 38, 36, 52,  0,  3, 41, 42,  7, 56,
     &                 37,  5,  6, 60, 61, 63, 64, 65, 69, 70,
     &                 11, 10, 50 /
c
c-----Cross-reference SOAP reactions to main mechanisms
c
c 15 reactions in this order
c
c 1       BNZA + OH
c 2       TOLA + OH
c 3       XYLA + OH
c 4       ISP + OH
c 5       ISP + O3
c 6       ISP + NO3
c 7       TRP + OH
c 8       TRP + O3
c 9       TRP + NO3
c 10      SQT + OH
c 11      SQT + O3
c 12      SQT + NO3
c 13      APN + OH
c 14      APN + O3
c 15      APN + NO3
c
      data iksoap_cb6 /167, 172, 177, 143, 149, 150, 164, 165, 166,
     &                 237, 238, 239,  -9,  -9,  -9 /
      data iksoap_cb7 /148, 152, 156, 178, 183, 184, 202, 206, 207,
     &                 208, 209, 210, 196, 200, 201 /
      data iksoap_rc2 / 81,  82,  84,  86, 130, 145, 88,  132, 147,
     &                  -9,  -9,  -9,  87, 131, 146 /
      data iksoap_s07 /536, 537, 540, 526, 527, 528, 558, 559, 560,
     &                 562, 563, 564, 530, 531, 532 /
c
c-----Parameters for rate constant expressions and diagnostics
c
      data npar     /  1,  2,  4, 10,  5, 12,  8 /
      data tdum     / 298.,  273.,  273. /
      data pdum     / 1013., 1013., 500. /
c
c     Many rules about names, number and order of species are 
c     enforced here unless LCHEM is false.
c
c     Arrays MCHGAS, MCHAERO, MCHRXN and MCHPHOT allow for up to 10
c     gas-phase mechanisms to be called by RADDRIVR and CHEMDRIV.
c     The current mechansims are:
c      1 = CB6r5
c      2 = Inactive
c      3 = CB6r5h (includes halogens)
c      4 = CB6r4 (includes DMS and iodine chemistry)
c      5 = SAPRC07TC
c      6 = Inactive
c      7 = CB7r1
c      8 = Inactive
c      9 = RACM2
c     10 = User-defined chemistry and species list
c
c     Radical species will be re-ordered according to the order defined
c     in RADLIST.
c     Only the species that are named in SPLIST will be allowed.
c     Update CHMDAT.INC & DDMCHM.INC if new species are added to SPLIST
c
c-----Entry point
c
      idmech = 0
      nphot1 = 0
      nphot2 = 0
      read(ichem,'(a)') record
      camxvin = record(21:30)
      call jstlft(camxvin)
      call toupper(camxvin)
      if (camxvin.ne.camxv) then
        write(iout,'(/,a)') ' CAMx version in CHEMPARAM file is INVALID'
        write(iout,'(a,a)') ' Expecting: ',camxv
        write(iout,'(a,a)') '     Found: ',camxvin
        goto 910
      endif
c
c-----Gas Mechanism ID
c
      read(ichem,'(a)') record
      read(record(21:80),'(a)',ERR=110,END=110) gasopt
      call jstlft(gasopt)
      call toupper(gasopt)
 110  continue
      if (lchem) then
        if (gasopt.NE.'CB6R5H'    .and.
     &      gasopt.NE.'CB6R4'     .and.
     &      gasopt.NE.'CB6R5'     .and.
     &      gasopt.NE.'CB7R1'     .and.
     &      gasopt.NE.'RACM2'     .and.
     &      gasopt.NE.'SAPRC07TC' .and.
     &      gasopt.NE.'MECH10') then
          write(iout,'(/,3a)') ' Invalid gas chemistry mechanism',
     &                       ' in CHEMPARAM input file: ',gasopt
          if( gasopt .EQ. ' ' ) then
            write(iout,'(a)') ' You must specify the gas mechanism.'
          endif
          write(iout,'(a,10(/,10x,a))') 'Acceptable options are: ',
     &                           'CB6r5H   ',
     &                           'CB6r4    ',
     &                           'CB6r5    ',
     &                           'CB7r1    ',
     &                           'RACM2    ',
     &                           'SAPRC07TC',
     &                           'MECH10'
          goto 910
        endif
        if (gasopt.EQ.'CB6R5') then
          idmech = 1 
          tuvlabl = 'TUV4.8CAMx7.30_CB6r5'
        elseif (gasopt.EQ.'CB6R5H') then
          idmech = 3
          tuvlabl = 'TUV4.8CAMx7.30_CB6r5'
        elseif (gasopt.EQ.'CB6R4') then
          idmech = 4
          tuvlabl = 'TUV4.8CAMx7.30_CB6r4'
        elseif (gasopt.EQ.'SAPRC07TC') then
          idmech = 5
          tuvlabl = 'TUV4.8CAMx7.30_SPR07'
        elseif (gasopt.EQ.'CB7R1') then
          idmech = 7
          tuvlabl = 'TUV4.8CAMx7.30_CB7r1'
        elseif (gasopt.EQ.'RACM2') then
          idmech = 9
          tuvlabl = 'TUV4.8CAMx7.30_RACM2'
        elseif (gasopt.EQ.'MECH10') then
          idmech = 10
          tuvlabl = 'TUV4.8CAMx7.30_OSPM'
        endif

        write(iout,'(/,2a)') 'Gas Chemistry Mechanism      :',gasopt
        write(idiag,'(a,i4)')'Mechanism ID                 :',idmech
 
        nrad = 0
        if (idmech.eq.10) then
          write(iout,'(/,a,/,a)') ' You selected Mechamism 10:',
     &         ' This will use your customized CAMx/chem10.f subroutine'
        elseif ((idmech.eq.1 .OR. idmech.eq.3 .OR. idmech.eq.4 .OR.
     &        idmech.eq.7) .AND. .NOT.lixemis .AND. .NOT.lixbypass) then
          write(iout,'(/,a)') ' WARNING in READCHM:'
          write(iout,'(a)')
     &      ' You selected CB6r5, CB6r5h, CB6r4 or CB7r1,'
          write(iout,'(a)')' BUT did not turn on in-line Ix emissions.'
          write(iout,'(2a)')' Be sure to provide I2 and HOI emissions',
     &                     ' externally.'
        endif
        nrad = mchrad(idmech)      ! need to know NRAD for advection, etc.
        if (nrad.GT.NRADNM) then
          write(iout,'(/,2a)')'Number of radical species is ',
     &                        'greater than internal dimension.'
          write(iout,'(a,i3)') 'Number of radical species  : ',nrad
          write(iout,'(a,i3)') 'Maximum dimension (NRADNM) : ',NRADNM
          write(iout,'(2a)')'Increase the parameter NRADNM ',
     &                                   'in READCHM.F and recompile.'
          goto 910
        endif
      endif
c
c-----Options for the aerosol treatment: 'NONE','INERT','CF'
c
      read(ichem,'(a)') record
      read(record(21:80),'(a)',ERR=111,END=111) aeropt
      call jstlft(aeropt)
      call toupper(aeropt)
 111  continue
      if (aeropt.NE.'CF'     .and. 
     &    aeropt.NE.'NONE'   .and.
     &    aeropt.NE.'INERT') then
        write(iout,'(/,3a)') ' Invalid option for aerosol treatment',
     &                       ' in CHEMPARAM input file: ',aeropt
        if( aeropt .EQ. ' ' ) then
           write(iout,'(a)') ' You must specify the aerosol treatment.'
        endif
        write(iout,'(a,4(/,10x,a))') 'Acceptable options are: ',
     &                           'NONE     - None (NAERO must = 0)',
     &                           'INERT    - Inert PM',
     &                           'CF       - Static Coarse/Fine scheme'
        goto 910
      endif
      if (idmech.EQ.10 .and. (aeropt.EQ.'CF')) then
        write(iout,'(/,3a)') ' Invalid option for aerosol treatment',
     &                       ' in CHEMPARAM input file: ',aeropt
        write(iout,'(a,2(/,10x,a))') 'Acceptable options for MECH10: ',
     &      'NONE  - None (NAERO must = 0)',
     &      'INERT - Inert PM (PM chemistry performed in CAMx/chem10.f)'
        goto 910
      endif
      write(idiag,'(2a)')  'Aerosol Treatment            :',aeropt
c
c-----Options for Inorganic PM chemistry: 'ISORROPIA' or 'EQSAM'
c
      read(ichem,'(a)') record
      read(record(21:80),'(a)',ERR=112,END=112) siaopt
      call jstlft(siaopt)
      call toupper(siaopt)
 112  if (aeropt.NE.'NONE' .and. aeropt.NE.'INERT') then
        if (siaopt.NE.'ISORROPIA'.and. siaopt.NE.'EQSAM') then
          write(iout,'(/,3a)') ' Invalid option for inorganic PM',
     &                         ' chemistry in CHEMPARAM input file: ',
     &                         siaopt
          if( siaopt .EQ. ' ' ) then
             write(iout,'(2a)') ' You must specify the inorganic PM',
     &                          ' chemistry.'
          endif
          write(iout,'(2a,2(/,10x,a))') 'Acceptable inorganic PM',
     &                                  ' chemistry options are: ',
     &                           'ISORROPIA - ISORROPIA scheme',
     &                           'EQSAM     - EQSAM4clim scheme'
          goto 910
        endif
      else
        write(iout,'(/,2a)') ' Inorganic PM chemistry is ignored for',
     &                       ' aerosol options NONE and INERT'
      endif
      leqsam = .false.
      if (siaopt.EQ.'EQSAM') leqsam = .true.
      write(idiag,'(2a)')  'Inorganic PM Chemistry       :',siaopt
c
c-----Options for SOA chemistry: 'SOAP3.0' or 'VBS'
c
      read(ichem,'(a)') record
      read(record(21:80),'(a)',ERR=113,END=113) soaopt
      call jstlft(soaopt)
      call toupper(soaopt)
 113  if (aeropt.NE.'NONE' .and. aeropt.NE.'INERT') then
        if(  soaopt.NE.'SOAP3_COMPLX' .and.
     &       soaopt.NE.'SOAP3_COMPLX_NOEVAP' .and.
     &       soaopt.NE.'SOAP3_SIMPLE' .and. 
     &       soaopt.NE.'SOAP3_SIMPLE_NOEVAP' .and.
     &       soaopt.NE.'1.5DVBS') then
          write(iout,'(/,3a)') ' Invalid option for SOA chemistry',
     &                         ' in CHEMPARAM input file: ',soaopt
          if( soaopt .EQ. ' ' ) then
             write(iout,'(a)') ' You must specify the SOA chemistry.'
          endif
          write(iout,'(a,5(/,10x,a))') 'Acceptable SOA options are: ',
     &       'SOAP3_SIMPLE         - SOAP3 non-volatile SOA, POA evap',
     &       'SOAP3_SIMPLE_NOEVAP  - SOAP3 non-volatile SOA, no POA evap',
     &       'SOAP3_COMPLX         - SOAP3 volatile SOA, POA evap',
     &       'SOAP3_COMPLX_NOEVAP  - SOAP3 volatile SOA, no POA evap',
     &       '1.5DVBS              - 1.5D Volatility Basis Set'
          goto 910
        endif
      else
        write(iout,'(/,2a)') ' SOA chemistry is ignored for aerosol',
     &                       ' options NONE and INERT'
      endif
      lvbs   = .false.
      luse_soap3 = .FALSE.
      luse_soap3_simple = .FALSE.
      luse_soap3_complx = .FALSE.
      if (soaopt.EQ.'1.5DVBS') lvbs = .true.
      if( soaopt(1:12) .EQ. 'SOAP3_SIMPLE' ) then
         luse_soap3_simple = .TRUE.
         luse_soap3 = .TRUE.
      endif
      if( soaopt(1:12) .EQ. 'SOAP3_COMPLX' ) then
         luse_soap3_complx = .TRUE.
         luse_soap3 = .TRUE.
      endif
      if( soaopt(1:5) .EQ. 'SOAP3') then
         if( soaopt(13:) .EQ. '_NOEVAP') then
           do n = IDX_SOAP_POA_OP, IDX_SOAP_POA_BB
             do i = 1,2
               soap3_evap_factors(i,n) = 0.
             enddo
             soap3_evap_factors(3,n) = 1.
           enddo
         else
           do n = IDX_SOAP_POA_OP, IDX_SOAP_POA_BB
             do i = 1,3
               soap3_evap_factors(i,n) = soap3_evap_factors_complx(i,n)
             enddo
           enddo
         endif
         if (idmech.eq.1 .OR. idmech.eq.3 .OR. idmech.eq.4) then
           iksoap = iksoap_cb6
         elseif (idmech.eq.5) then
           iksoap = iksoap_s07
         elseif (idmech.eq.7) then
           iksoap = iksoap_cb7
         elseif (idmech.eq.9) then
           iksoap = iksoap_rc2
         elseif (idmech.eq.10) then
           write(iout,'(/,a)')'SOAP cannot be run with Mechanism 10'
           goto 910
         else
           write(iout,'(/,a,i3)')
     &     'Dont know how to set SOAP rxns for mech #', idmech
           goto 910
         endif
      endif
      write(idiag,'(2a)')  'Organic PM Chemistry         :',soaopt
c
      read(ichem,'(a)') record
      write(idiag,'(2a)') 'Description                  :',
     &                    record(21:istrln(record))
      ncf_chemistry = record(21:istrln(record))
c
c-----Number of species
c
      read(ichem,'(a)') record
      read(record(21:80),*) ngas
      read(ichem,'(a)') record
      read(record(21:80),*) naero
      nspec = ngas + naero
      if (nspec.lt.1) then
        write(iout,'(/,a,i5,a)')
     &    ' Number of GAS + AERO species on CHEMPARAM file =', nspec,
     &    ' is less than 1'
          goto 910
      endif
      if (naero.EQ.0 .AND. (aeropt.EQ.'CF' .OR.
     &                      aeropt.EQ.'INERT')) then
        write(iout,'(//A)')
     &  ' Number of aerosols = 0 for CF or INERT PM options'
        write(iout,'(a)') ' Set aerosol treatment to NONE'
        goto 910
      endif
      if (naero.GT.0) then
        if (aeropt.EQ.'NONE') then
          write(iout,'(//A)')
     &      ' Aerosol species are requested with aerosol option = NONE'
          write(iout,'(a,2(/,10x,a))') 'Acceptable aerosol options: ',
     &                            'INERT  - Inert PM',
     &                            'CF     - Static Coarse/Fine scheme'
          goto 910
        endif
        read(record(21:),*,ERR=903,END=903) naero,dtaero,nsec_c,
     &                      (dsec_i(i),i=1,nsec_c+1)
        nbin = nsec_c
        write(idiag,'(a,i4)')   'Number of Aerosol Species    :',naero
        write(idiag,'(a,f4.0)') 'Aerosol Coupling Freq (min)  :',dtaero
        write(idiag,'(a,i4)')   'Number of Aerosol Size bins  :',nbin
      endif
      if( nspec .GT. MXSPEC ) then
         write(iout,'(//,A)') 'ERROR in READCHM:'
         write(iout,'(2a)') 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
         write(iout,*) 'Please change the value for parameter MXSPEC: ',MXSPEC
         write(iout,*) 'It should be set to a value of at least: ',nspec
         goto 910
      endif
c
c-----Read number of reactions and check that the dimension parameter is large enough
c
      read(ichem,'(a)') record
      read(record(21:80),*) nreact
      if( nreact .GT. MXREACT ) then
         write(iout,'(//,A)') 'ERROR in READCHM:'
         write(iout,'(2a)') 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
         write(iout,*) 'Please change the value for parameter MXREACT: ',MXREACT
         write(iout,*) 'It should be set to a value of at least: ',nreact
         goto 910
      endif
c
c-----Check number of species and number of reactions against chosen mechanism
c 
      if (lchem) then
c        if (ngas.gt.mchgas(idmech)) then
c          write(iout,'(/,a,i5,/,a,i5)')
c     &     ' Number of GAS species in CHEMPARAM file =', ngas,
c     &     ' greater than ngas for this mechanism =', mchgas(idmech)
c          goto 910
c        endif
c
c       if (aeropt.NE.'NONE' .AND. aeropt.NE.'INERT' .AND.
c    &      naero.gt.mchaero(idmech)) then
c         write(iout,'(/,a,i5,/,a,i5)')
c    &     ' Number of AERO species on CHEMPARAM file =', naero,
c    &     ' Max allowed for this mechanism =', mchaero(idmech)
c         goto 910
c       endif
c
        if (nreact.ne.mchrxn(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' Number of reactions on CHEMPARAM file =',nreact,
     &     ' not equal to reactions for this mechanism =',mchrxn(idmech)
          goto 910
        endif
      endif
c
c-----Read primary photolysis ID record
c
      read(ichem,'(a)') recpht1
      if( lchem ) read(recpht1(21:80),*) nphot1
c
c-----Read secondary (scaled) photolysis ID records
c
      read(ichem,'(a)') record
      read(record(21:80),*) nphot2
      if( lchem ) then
         if (nphot2.gt.0) then
           if (nphot1.eq.0) then
             write(iout,'(/,a)')
     &           'Need at least one primary photolysis reaction'
             goto 910
           endif
         endif
      endif
      nphot = nphot1 + nphot2
c
c-----Call routine to allocate the arrays
c
      call alloc_chmstry(ngrid,nspec,nemiss_files,npoint_files,
     &                                      nreact,nphot1,nphot2)
      mole_weight = 1.0
c
c-----Check photolysis ID records
c
      if( nphot .GT. 0 ) then
c
        if( idmech .GT. 0 ) then
           if(nphot.ne.mchphot(idmech)) then
             write(iout,'(/,a,i5,/,a,i5)')
     &        ' Chemistry mechanism requires', mchphot(idmech),
     &        ' photolysis reactions, but CHEMPARAM file has' , nphot
             goto 910
           endif
        endif
c
        if( lchem ) then
           if (nphot.ne.mchphot(idmech)) then
             write(iout,'(/,a,i5,/,a,i5)')
     &           ' Chemistry mechanism requires', mchphot(idmech),
     &        ' photolysis reactions, but CHEMPARAM file has' , nphot
             goto 910
           endif
        endif
c
        if (nphot1.gt.MXPHT1) then
          write(iout,'(2a)')'Number of photolysis reactions is ',
     &                              'greater than internal dimensions.'
          write(iout,'(a,i3)') 'Number of reactions needed  : ',nphot1
          write(iout,'(a,i3)') 'Maximum dimension (MXPHT1)  : ',MXPHT1
          write(iout,'(2a)')'Increase the parameter MXPHT1 ',
     &                                   'in CAMx.PRM and recompile.'
          goto 910
        endif
c
        if (nphot2.gt.MXPHT2) then
          write(iout,'(/,a,i5)') 
          write(iout,'(2a)')'Number of photolysis reactions is ',
     &                              'greater than internal dimensions.'
          write(iout,'(a,i3)') 'Number of reactions needed  : ',nphot2
          write(iout,'(a,i3)') 'Maximum dimension (MXPHT2)  : ',MXPHT2
          write(iout,'(2a)')'Increase the parameter MXPHT2 ',
     &                                   'in CAMx.PRM and recompile.'
          goto 910
        endif
c
c-----Read primary photolysis ID record
c
        if (nphot1.gt.0) then
          read(recpht1(21:180),*) nphot1,(idphot1(n),n=1,nphot1)
          write(idiag,'(a)') 'The primary photolysis reactions are'
          write(idiag,'(i6)') (idphot1(n),n=1,nphot1)
        endif
c
c-----Read secondary (scaled) photolysis ID records
c
        if (nphot2.gt.0) then
          do n = 1,nphot2
            read(ichem,'(a)') record
            if( lchem ) read(record(21:80),*) idphot2(n),idphot3(n),phtscl(n)
          enddo
          if( lchem ) then
             write(idiag,'(a)') 'The secondary photolysis reactions are'
             write(idiag,'(i6,a,i6,a,1pe10.3)')
     &            (idphot2(n),' =',idphot3(n),' *',phtscl(n),n=1,nphot2)
           endif
        endif
c
        if( lchem ) then
           nerr = 0
           do n = 1,nphot1
             if (idphot1(n).gt.nreact) nerr = 1
           enddo
c
           do n = 1,nphot2
             if (idphot2(n).gt.nreact) nerr = 1
             if (idphot3(n).gt.nreact) nerr = 1
             nhit = 0
             do nn = 1,nphot1
               if (idphot3(n).eq.idphot1(nn)) nhit = 1
             enddo
             if (nhit.eq.0) nerr = 1
           enddo
c
           if (nerr.eq.1) then
              write(iout,'(/,a)') 'ERROR in the CHEMPARAM file:'
              write(iout,'(a)')
     &        ' Bad reaction number in one of the photolysis reaction IDs.'
             goto 910
           endif
        endif
      endif
c
c-----Surface model flags
c
      read(ichem,'(a)') record
      read(record(21:80),*) nsmspc,nsmrxn
      if (lsrfmod) then
         if (nsmspc.eq.0 .or. nsmrxn.eq.0) then
           write(iout,'(/,a)') 'ERROR in READCHM'
           write(iout,'(3a)') 'Surface Model is invoked but chemistry',
     &                        ' parameters file specifies zero surface',
     &                        ' model species or reactions.'
           write(iout,'(a)') 'Review your chemistry parameters file'
           goto 910
         endif
         if (nsmspc.gt.MXSMSPC) then
           write(iout,'(/,a)') 'ERROR in READCHM'
           write(iout,'(2a)') 'Number of Surface Model species exceeds',
     &                        ' maximum'
           write(iout,'(a,i5)') 'Reset MXSMSPC in CAMx.prm to ',nsmspc
           write(iout,'(a)') 'Then re-compile CAMx.'
           goto 910
         endif
         call alloc_srfmod(nsmspc,nsmrxn)
      else
         nsmspc = 0
         nsmrxn = 0
      endif
c
c-----Species records, gases come first
c
c      write(idiag,'(a)') 'The state species are'
      write(idiag,'(a,i5)')
     &     ' Number of GAS species in CHEMPARAM file =', ngas
      read(ichem,'(a)') record
      if (ngas.gt.0) then
        read(ichem,'(a)') record
        write(idiag,'(a)') record(:istrln(record))
        do l=1,ngas
          read(ichem,'(5x,a10,2e10.0,4f10.0)')
     &               spname(l),bdnl(l),henry0(l),tfact(l),
     &               molwt_spc,f0(l),rscale(l)
          diffrat(l) = sqrt(molwt_spc/18.)
          rscale(l) = min(1.,rscale(l))
          rscale(l) = max(0.,rscale(l))
          mole_weight(l) = max(1.,molwt_spc)
          write(idiag,'(i3,2x,a10,2e10.2,4f10.2)')
     &               l,spname(l),bdnl(l),henry0(l),tfact(l),
     &               mole_weight(l),f0(l),rscale(l)
          if (istrln(spname(l)).gt.8) goto 907
        enddo
      endif
c
c-----Check over the gas phase species for deposition calculations
c
      if ((ldry .or. lwet) .and. ngas.gt.0) then
        do l = 1,ngas
          if (henry0(l).eq.0.) then
            write(iout,'(/,a,i5)')
     &      'The Henry0 value must be non-zero for species ',l
            goto 910
          endif
        enddo
      endif
c
c-----Reorder radicals to come first
c
      if (lchem .and. ngas.gt.0 .and. idmech.NE.10) then
        if (nrad.gt.0) then
          do i=1,nrad
            ierr = 1
            do j=1,ngas
              if (radlist(i,idmech).EQ.spname(j)) then
                nametmp = spname(i)
                spname(i) = spname(j)
                spname(j) = nametmp
                tmp = bdnl(i)
                bdnl(i) = bdnl(j)
                bdnl(j) = tmp
                tmp = henry0(i)
                henry0(i) = henry0(j)
                henry0(j) = tmp
                tmp = tfact(i)
                tfact(i) = tfact(j)
                tfact(j) = tmp
                tmp = diffrat(i)
                diffrat(i) = diffrat(j)
                diffrat(j) = tmp
                tmp = f0(i)
                f0(i) = f0(j)
                f0(j) = tmp
                tmp = rscale(i)
                rscale(i) = rscale(j)
                rscale(j) = tmp
                tmp = mole_weight(i)
                mole_weight(i) = mole_weight(j)
                mole_weight(j) = tmp
                ierr = 0
                EXIT
              endif
            enddo
            if (ierr.eq.1) then
              write(iout,'(/,2a)') 'You must include all of the ',
     &                'following radical species in CHEMPARAM file:'
              write(iout,'(6a10)') (radlist(l,idmech),l=1,nrad)
              goto 910
            endif
          enddo
        endif
      endif
c
c-----Read the aero species, if any
c
      if (naero.ne.0) then
        read(ichem,'(a)') record
        write(idiag,'(a)') record(:istrln(record))
        do l = ngas+1,nspec
          read(ichem,'(5x,a10,e10.0,2f10.0,i10,f10.0,i10)')
     &       spname(l),bdnl(l),roprt(l),bext(l),rhadj(l),ssa(l),ibin
          if (ibin .lt. 1 .OR. ibin.gt.nbin) then
            write(iout,'(/,3a)') 'PM size bin for species: ',
     &                           spname(l)(:istrln(spname(l))),
     &                       ' is outside range from CHEMPARAM header.'
            write(iout,'(a,i5)') 'PM size bin read:       ',ibin
            write(iout,'(a,i5)') 'Number of PM size bins: ',nbin
            goto 910
          endif
          dcut(l,1) = REAL(dsec_i(ibin))
          dcut(l,2) = REAL(dsec_i(ibin+1))
          write(idiag,'(i3,2x,a10,e10.2,2f10.2,i10,3f10.2)')
     &       l,spname(l),bdnl(l),roprt(l),bext(l),rhadj(l),
     &       ssa(l),(dcut(l,m),m=1,2)
          bext(l) = bext(l)*1.e-6
          roprt(l) = roprt(l)*1.e6
          if (istrln(spname(l)).gt.8) goto 907
        enddo
      endif
c
c-----Map species names to the internal name list.  The default setting
c     to nspec+1 allows species that are in the chem solvers to be
c     omitted from the species list for this run.  Testing if a named
c     pointer is set to nspec+1 is used to identify species that are not
c     in the run.  Named pointers (e.g., kno) are set by equivalence
c     in CHMDAT.INC
c
      if (lchem .and. idmech.NE.10) then
        do j = 1,NSPNAM
          kmap(j) = nspec+1
        enddo
c
        nn = nspec
        if (aeropt.eq.'INERT') nn = ngas
        do 10 j = 1,nn
          do i = 1,NSPNAM
            if (splist(i).eq.spname(j)) then
              kmap(i) = j
              goto 10
            endif
          enddo
          write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file.'
          write(iout,'(3a)') 'Species in chemparam file ',
     &     spname(j)(:istrln(spname(j))),' is not in the internal list.'
          write(iout,'(2a)') 'Please make sure you are using the ',
     &               'correct chemparam file for this version of CAMx.'
          goto 910
  10    continue
        spname(nspec+1) = blank
        bdnl(nspec+1) = 0.
        write(idiag,'(a)')'Internal species order is:'
        write(idiag,'(i4,1x,a10)') (i,spname(i),i=1,nspec)
c
c-----Check NH3 if Bi-Di NH3 is on
c
        if (lbidinh3 .and. knh3.eq.nspec+1) then
          write(iout,'(/,2a)') ' You must have NH3 to use the',
     &                         ' Bi-Di NH3 drydep scheme.'
          goto 910
        endif
c
c-----Check O3 if stratospheric ozone profile scheme is on
c
        if (lstrato3 .and. ko3.eq.nspec+1) then
          write(iout,'(/,2a)') ' You must have O3 to use the',
     &                         ' stratospheric ozone profile scheme.'
          goto 910
        endif
c
c-----Check species list for consistency with CF option
c
        if (naero.GT.0 .and. aeropt.NE.'INERT') then
c
c-----Check mandatory gas species for PM chemistry
c
          if (lvbs) then
            if (ksulf.eq.nspec+1 .or. khno3.eq.nspec+1 .or.
     &           knh3.eq.nspec+1 .or.  kso2.eq.nspec+1 .or.
     &          kh2o2.eq.nspec+1 .or.   ko3.eq.nspec+1 .or.
     &          kn2o5.eq.nspec+1 .or.
     &          kivoa.eq.nspec+1 .or. kivob.eq.nspec+1 .or.
     &          kivod.eq.nspec+1 .or. kivog.eq.nspec+1 .or.
     &          kvap1.eq.nspec+1 .or. kvap2.eq.nspec+1 .or.
     &          kvap3.eq.nspec+1 .or. kvap4.eq.nspec+1 .or.
     &          kvas1.eq.nspec+1 .or. kvas2.eq.nspec+1 .or.
     &          kvas3.eq.nspec+1 .or. kvas4.eq.nspec+1 .or.
     &          kvcp1.eq.nspec+1 .or. kvcp2.eq.nspec+1 .or.
     &          kvcp3.eq.nspec+1 .or. kvcp4.eq.nspec+1 .or.
     &          kvfp1.eq.nspec+1 .or. kvfp2.eq.nspec+1 .or.
     &          kvfp3.eq.nspec+1 .or. kvfp4.eq.nspec+1 .or.
     &          kvbs1.eq.nspec+1 .or. kvbs2.eq.nspec+1 .or.
     &          kvbs3.eq.nspec+1 .or. kvbs4.eq.nspec+1) then
              write(iout,'(/,2a)') ' You must have all of the',
     &                ' following gas species to use the VBS option:'
              write(iout,'(a)')' SULF, HNO3, NH3, SO2, H2O2, O3, N2O5,'
              write(iout,'(a)')' IVOA, IVOB, IVOD, IVOG,'
              write(iout,'(a)')' VAP1, VAP2, VAP3, VAP4,'
              write(iout,'(a)')' VAS1, VAS2, VAS3, VAS4,'
              write(iout,'(a)')' VCP1, VCP2, VCP3, VCP4,'
              write(iout,'(a)')' VFP1, VFP2, VFP3, VFP4,'
              write(iout,'(a)')' VBS1, VBS2, VBS3, VBS4.'
              goto 910
            endif
          else
            if (ksulf.eq.nspec+1  .or. khno3.eq.nspec+1 .or.
     &           knh3.eq.nspec+1  .or.  kso2.eq.nspec+1 .or.
     &          kh2o2.eq.nspec+1  .or.   ko3.eq.nspec+1 .or.
     &          kn2o5.eq.nspec+1  .or. kivoa.eq.nspec+1 .or.
     &           kacg2.eq.nspec+1 .or.  kacg1.eq.nspec+1 .or.
     &           kbcg2.eq.nspec+1 .or.  kbcg1.eq.nspec+1 .or.
     &           ksvoa.eq.nspec+1) then
              write(iout,'(/,2a)') ' You must have all of the',
     &               ' following gas species to use the CF option:'
              write(iout,'(a)')' SULF, HNO3, NH3, SO2, H2O2, O3, N2O5,'
              write(iout,'(a)')' IVOA, SVOA, ACG2, ACG1, BCG2, BCG1.'
              goto 910
            endif
          endif

          if (aeropt.EQ.'CF') then
            if (lvbs) then
c
c-----Check mandatory PM species for VBS mechanism
c
              if (kpso4.eq.nspec+1 .or. kpno3.eq.nspec+1 .or.
     &            kpnh4.eq.nspec+1 .or. kph2o.eq.nspec+1 .or.
     &            kpap0.eq.nspec+1 .or.
     &            kpap1.eq.nspec+1 .or. kpap2.eq.nspec+1 .or.
     &            kpap3.eq.nspec+1 .or. kpap4.eq.nspec+1 .or.
     &            kpas0.eq.nspec+1 .or.
     &            kpas1.eq.nspec+1 .or. kpas2.eq.nspec+1 .or.
     &            kpas3.eq.nspec+1 .or. kpas4.eq.nspec+1 .or.
     &            kpcp0.eq.nspec+1 .or.
     &            kpcp1.eq.nspec+1 .or. kpcp2.eq.nspec+1 .or.
     &            kpcp3.eq.nspec+1 .or. kpcp4.eq.nspec+1 .or.
     &            kpfp0.eq.nspec+1 .or.
     &            kpfp1.eq.nspec+1 .or. kpfp2.eq.nspec+1 .or.
     &            kpfp3.eq.nspec+1 .or. kpfp4.eq.nspec+1 .or.
     &            kpbs0.eq.nspec+1 .or.
     &            kpbs1.eq.nspec+1 .or. kpbs2.eq.nspec+1 .or.
     &            kpbs3.eq.nspec+1 .or. kpbs4.eq.nspec+1) then
                write(iout,'(/,2a)') ' You must have all of the',
     &                 ' following PM species to use the VBS option:'
                write(iout,'(a)')' PSO4, PNO3, PNH4, PH2O,'
                write(iout,'(a)')' PAP0, PAP1, PAP2, PAP3, PAP4,'
                write(iout,'(a)')' PAS0, PAS1, PAS2, PAS3, PAS4,'
                write(iout,'(a)')' PCP0, PCP1, PCP2, PCP3, PCP4,'
                write(iout,'(a)')' PFP0, PFP1, PFP2, PFP3, PFP4,'
                write(iout,'(a)')' PBS0, PBS1, PBS2, PBS3, PBS4.'
                goto 910
              endif
            else
c
c-----Check mandatory PM species for CF mechanism
c
              if (kpso4.eq.nspec+1 .or. kpno3.eq.nspec+1 .or.
     &            kpnh4.eq.nspec+1 .or. kph2o.eq.nspec+1 .or.
     &            kaoa2.eq.nspec+1 .or. kaoa1.eq.nspec+1 .or.
     &            kboa2.eq.nspec+1 .or. kboa1.eq.nspec+1 .or.
     &            kaoa0.eq.nspec+1 .or. kboa0.eq.nspec+1) then
                write(iout,'(/,2a)') ' You must have all of the',
     &                     ' following PM species to use the CF option:'
                write(iout,'(a)')' PSO4, PNO3, PNH4, PH2O,'
                write(iout,'(a)')' AOA2, AOA1, BOA2, BOA1, AOA0, BOA0.'
                goto 910
              endif
            endif
c
c-----Check for a consistent set of sea salt species, or none
c
            if (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                               khcl.eq.nspec+1) then
              continue
            elseif (kna.lt.nspec+1 .and. kpcl.lt.nspec+1 .and.
     &                                   khcl.lt.nspec+1) then
              continue
            elseif (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                                   khcl.lt.nspec+1) then
              continue
            else
              write(iout,'(/,2a)') ' You must have all or none of the',
     &               ' sea salt species to use the CF option:'
              write(iout,'(a)')' NA, PCL, HCL.'
              write(iout,'(a)')' Or, just HCL.'
              goto 910
            endif
          endif
        else
c
c-----Check for gas species for non-aerosol chemistry mechanisms
c
          if (kivoa.lt.nspec+1 .or. kivob.lt.nspec+1 .or.
     &        kivod.lt.nspec+1 .or. kivog.lt.nspec+1 .or.
     &        kacg2.lt.nspec+1 .or. kacg1.lt.nspec+1 .or.
     &        kbcg2.lt.nspec+1 .or. kbcg1.lt.nspec+1  .or.
     &        kvap1.lt.nspec+1 .or. kvap2.lt.nspec+1 .or.
     &        kvap3.lt.nspec+1 .or. kvap4.lt.nspec+1 .or.
     &        kvas1.lt.nspec+1 .or. kvas2.lt.nspec+1 .or.
     &        kvas3.lt.nspec+1 .or. kvas4.lt.nspec+1 .or.
     &        kvcp1.lt.nspec+1 .or. kvcp2.lt.nspec+1 .or.
     &        kvcp3.lt.nspec+1 .or. kvcp4.lt.nspec+1 .or.
     &        kvfp1.lt.nspec+1 .or. kvfp2.lt.nspec+1 .or.
     &        kvfp3.lt.nspec+1 .or. kvfp4.lt.nspec+1 .or.
     &        kvbs1.lt.nspec+1 .or. kvbs2.lt.nspec+1 .or.
     &        kvbs3.lt.nspec+1 .or. kvbs4.lt.nspec+1 .or.
     &        ksvoa.lt.nspec+1) then
            write(iout,'(/,a)') ' You must have none of the species'
            write(iout,'(a)')   ' IVOA, IVOB, IVOD, IVOG,'
            write(iout,'(a)')   ' ACG2, ACG1, BCG2, BCG1,'
            write(iout,'(a)')   ' VAP1, VAP2, VAP3, VAP4,'
            write(iout,'(a)')   ' VAS1, VAS2, VAS3, VAS4,'
            write(iout,'(a)')   ' VCP1, VCP2, VCP3, VCP4,'
            write(iout,'(a)')   ' VFP1, VFP2, VFP3, VFP4,'
            write(iout,'(a)')   ' VBS1, VBS2, VBS3, VBS4,'
            write(iout,'(a)')   ' SVOA'
            write(iout,'(a)')   ' with non-aerosol chemistry mechanisms.'
            goto 910
          endif
        endif
c
c-----Check for a consistent set of Mercury species, or none
c
        if (khg0.eq.nspec+1 .and. khg2.eq.nspec+1 .and.
     &      khgiip.eq.nspec+1 .and. khgiipc.eq.nspec+1) then
          continue
        elseif (khg0.lt.nspec+1 .and. khg2.lt.nspec+1 .and.
     &      khgiip.lt.nspec+1 .and. khgiipc.lt.nspec+1) then
          continue
        else
          write(iout,'(/,a)') ' You must have all or none of the '
          write(iout,'(a)')   ' Mercury gas-phase and adsorbed species '
          write(iout,'(a)')   ' HG0, HG2, HGIIP, HGIIPC'
          goto 910
        endif
        if (khg0.lt.nspec+1 .and. naero.lt.2 ) then
          write(iout,'(/,2a)') ' You must model PM10 to use the',
     &                         ' Mercury gas-phase chemistry.'
          write(iout,'(2a)')   ' HG0 and HG2 are included but',
     &          ' the number of aerosol species is less than 2'
          goto 910
        endif
        if (khg0.lt.nspec+1 .and. aeropt.ne.'CF') then
          write(iout,'(/,2a)') ' Mercury is currently available',
     &                         ' only with the CF option.'
          goto 910
        endif
        if (khg0.lt.nspec+1 .and. khcl.eq.nspec+1) then
          write(iout,'(/,2a)') ' Mercury chemistry requires HCl'
          goto 910
        endif
        if (khg0.lt.nspec+1 .and. lcig(1)) then
          write(iout,'(/,2a)') ' Mercury chemistry cannot be run',
     &                         ' with the sub-grid cloud model'
          write(iout,'(/,a)') ' Turn off sub-grid cloud option'
          goto 910
        endif
      endif
c
c-----Set up section diameters and check parameters for CF routine
c
      if (lchem .AND. aeropt.EQ.'CF') then
        ierr = 0
        call aeroset(dsec_i,ierr)
        if ( ierr .ne. 0 ) goto 910
      endif
c
c-----Reaction records
c
      if (lchem .AND. nreact.gt.0) then
        ncms = 0
        nppm = 0
        read(ichem,'(a)') record
        read(ichem,'(a)') record
        write(idiag,'(a)') 'The reaction rate parameters are'
        do i = 1,nreact
          read(ichem,'(a)') record
          read(record,*) num,rxntyp(i)
          if (abs(rxntyp(i)).lt.1 .or. abs(rxntyp(i)).gt.7) goto 902
          read(record,*,err=900) num,rxntyp(i),rxnord(i),
     &                           (rxnpar(i,j),j=1,npar(abs(rxntyp(i))))
c
          if (abs(rxntyp(i)).eq.5) then
            iref = nint(rxnpar(i,1))
            if (i.le.iref) goto 901
          endif
          if (abs(rxntyp(i)).eq.2 .or. abs(rxntyp(i)).eq.3) then
            if (rxnpar(i,1).gt.1.0D0) then
              nppm = nppm +1
            else
              ncms = ncms +1
            endif
          endif
c
          write(idiag,'(3i3,1p12d12.4)') num,rxntyp(i),rxnord(i),
     &                           (rxnpar(i,j),j=1,npar(abs(rxntyp(i))))
        enddo
c
c-----Decide what units the rate constants are provided in
c
        if (nppm.eq.0 .and. ncms.eq.0) then
          write(iout,'(/,a)') ' Unable to determine rate constant units'
          write(iout,*)       ' ncms, nppm = ', ncms, nppm
          goto 910
        elseif (nppm.gt.ncms) then
          write(idiag,'(/,a)') 'Rate constants input in ppm units'
          write(idiag,*)       'ncms, nppm = ', ncms, nppm
          iunits = 1
        else
          write(idiag,'(/,a)') 'Rate constants input in cms units'
          write(idiag,*)       'ncms, nppm = ', ncms, nppm
          iunits = 2
        endif
c
c-----Populate rate constant lookup table
c
        call exptbl(iunits,rxntyp,rxnord,rxnpar)
c
c-----Provide diagnostic info for checking rate expressions
c
        write(idiag,'(/,a,/,/,a)') 
     &        'Diagnostic info for checking rate expressions',
     &        'Rates at three temps and pressures in cm3 molec-1 s-1'
        do i=1,3
          call ktherm(tdum(i), pdum(i))
          do j=1,nreact
            kdum(j,i)=rk(j)/
     &           ( ( (298.0/tdum(i)) * (pdum(i)/1013.0) *
     &                   2.462e13 )**(rxnord(j)-1) ) / 3600.  ! cm3 molec-1 s-1
          enddo
        enddo
        write(idiag,'(a,3F12.1)') 'Temp= ', tdum
        write(idiag,'(a,3F12.1)') 'Pres= ', pdum
        write(idiag,'(a)') 'Rxn Type'
        write(idiag,'(2i3,1p3e12.4)')
     &       (j,rxntyp(j),kdum(j,1),
     &        kdum(j,2),kdum(j,3),j=1,nreact)
c
c-----Set pointers for photolysis reactions to receive T,P adjustments
c
        if (idmech.eq.5) then
          jno2rxn    = jno2s07
          jo3rxn     = jo3s07
          jhcho1rxn  = jhcho1s07
          jhcho2rxn  = jhcho2s07
          jch3chorxn = jch3chos07
        elseif (idmech.eq.7) then
          jno2rxn    = jno2cb7
          jo3rxn     = jo3cb7
          jhcho1rxn  = jhcho1cb7
          jhcho2rxn  = jhcho2cb7
          jch3chorxn = jch3chocb7
        elseif (idmech.eq.1 .OR. idmech.eq.3 .OR. idmech.eq.4) then
          jno2rxn    = jno2cb6
          jo3rxn     = jo3cb6
          jhcho1rxn  = jhcho1cb6
          jhcho2rxn  = jhcho2cb6
          jch3chorxn = jch3chocb6
        elseif (idmech.eq.9) then
          jno2rxn    = jno2rc2
          jo3rxn     = jo3rc2
          jhcho1rxn  = jhcho1rc2
          jhcho2rxn  = jhcho2rc2
          jch3chorxn = jch3chorc2
        else
          write(iout,'(/,a,i3)')
     &       'Must set T & P photolysis rxn indices for mech #', idmech
          goto 910
        endif
c
c-----Set pointers to heterogeneous rxns
c
        if (idmech.eq.1) then
          ihetrxn(1) = ihet01_cb6
          ihetrxn(2) = -999
          ihetrxn(3) = ihet03_cb6r4
          ihetrxn(4) = ihet04_cb6r4
          ihetrxn(5) = ihet05_cb6r4
          ihetrxn(6) = ihet06_cb6r5
          ihetrxn(7) = ihet07_cb6r5
        elseif (idmech.eq.3) then
          ihetrxn(1) = ihet01_cb6
          ihetrxn(2) = ihet02_cb6r5h
          ihetrxn(3) = ihet03_cb6r4
          ihetrxn(4) = -999
          ihetrxn(5) = -999
          ihetrxn(6) = ihet06_cb6r5h
          ihetrxn(7) = ihet07_cb6r5h
        elseif (idmech.eq.4) then
          ihetrxn(1) = ihet01_cb6
          ihetrxn(2) = -999
          ihetrxn(3) = ihet03_cb6r4
          ihetrxn(4) = ihet04_cb6r4
          ihetrxn(5) = ihet05_cb6r4
          ihetrxn(6) = ihet06_cb6r4
          ihetrxn(7) = ihet07_cb6r4
        elseif (idmech.eq.5) then
          ihetrxn(1) = ihet01_s07
          ihetrxn(2) = -999
          ihetrxn(3) = -999
          ihetrxn(4) = -999
          ihetrxn(5) = -999
          ihetrxn(6) = ihet06_s07
          ihetrxn(7) = ihet07_s07
        elseif (idmech.eq.7) then
          ihetrxn(1) = ihet01_cb7
          ihetrxn(2) = -999
          ihetrxn(3) = ihet03_cb7
          ihetrxn(4) = -999
          ihetrxn(5) = ihet05_cb7
          ihetrxn(6) = ihet06_cb7
          ihetrxn(7) = ihet07_cb7
        elseif (idmech.eq.9) then
          ihetrxn(1) = ihet01_rc2
          ihetrxn(2) = -999
          ihetrxn(3) = ihet03_rc2
          ihetrxn(4) = ihet04_rc2
          ihetrxn(5) = ihet05_rc2
          ihetrxn(6) = ihet06_rc2
          ihetrxn(7) = ihet07_rc2
        endif
        if (ihetrxn(1) .LE. 0) then
          write(iout,'(/,a,i3)')
     &       'Must set N2O5 hydrolysis rxn index for mech #', idmech
          goto 910
        endif
c
c-----Set reaction pointers for pig chemistry
c
        if (ipigflg .NE. 0) then
           if (naero.GT.0) then
             if (aeropt.EQ.'CF' .AND. ipigflg .EQ. IRONPIG) then
               write(iout,'(/,2a)')'IRON PiG cannot be run with',
     &                            ' CF PM option'
               write(iout,'(/,2a)')'Either turn off PiG or change',
     &                            ' to GREASD'
               goto 910
             endif
           endif
           if(idmech.eq.5) then
             ipigrxn = ipigs07
             lrmap   = lrmaps07
           elseif(idmech.eq.1 .OR. idmech.eq.3 .OR. idmech.eq.4) then
             ipigrxn = ipigcb6
             lrmap   = lrmapcb6r4
           elseif(idmech.eq.7) then
             ipigrxn = ipigcb7
             lrmap   = lrmapcb7
           elseif(idmech.eq.9) then
             ipigrxn = ipigrc2
             lrmap   = lrmaprc2
           elseif (idmech.eq.10) then
             write(iout,'(/,a)')'PiG cannot be run with Mechanism 10'
             goto 910
           else
             write(iout,'(/,a,i3)') 
     &       'Dont know how to set pig rxns for mech #', idmech
             goto 910
           endif

c ----- map species indexes for new GREASD PiG chemistry solver

           lsmap( 1) = kO
           lsmap( 2) = kO1D
           lsmap( 3) = kOH
           lsmap( 4) = kHO2
           lsmap( 5) = kNO3
           lsmap( 6) = kN2O5
           lsmap( 7) = kNO
           lsmap( 8) = kNO2
           lsmap( 9) = kO3
           lsmap(10) = kHONO
           lsmap(11) = kHNO3
           lsmap(12) = kCO
           if (idmech.eq.5 .or. idmech.eq.9) then
             lsmap(13) = kHCHO
           else
             lsmap(13) = kFORM
           endif
           lsmap(14) = kSO2
           lsmap(15) = kSULF
        endif
      endif
c
c======================== DDM Begin =======================
c
c---  Set species pointers via equivalences in ddmchm.inc
c     LMAP(non-model spc) = NGAS + 1 while KMAP(non-model spc) = NSPEC + 1,
c     which allows shorter local arrays in the chemistry solver routines.
c
      if (lchem .and. idmech.NE.10) then
        do i = 1, nspnam
          if (kmap(i).gt.ngas) then
            lmap(i) = ngas + 1
          else
            lmap(i) = kmap(i)
          endif
        enddo
      endif
c
c======================== DDM End   =======================
c
c  --- set the flag for determining if a species is gaseous ---
c
      do i=1,nrad
        lgas(i) = .FALSE.
      enddo
      do i=nrad+1,ngas
        lgas(i) = .TRUE.
      enddo
      do i=ngas+1,nspec
        lgas(i) = .FALSE.
      enddo
c
c-----Check for and process surface model records
c
      if (lsrfmod) then
        read(ichem,'(a)',err=904,end=904) record
        if (record(1:13).ne.'Surface Model') goto 904
c
c-----Species list and associated parameters
c
        read(ichem,'(a)',err=905,end=906) record
        do n = 1,nsmspc
          read(ichem,'(a)',err=905,end=906) record
          read(record,'(5x,a10,6e10.0)',err=905) smspc(n),smssrb(n),
     &                                  smlch(n),smvsrb(n),smpen(n),
     &                                  smisrb(n),smmlt(n)
          idsmsp(n) = 0
          do l = 1,nspec
            if (smspc(n).eq.spname(l)) idsmsp(n) = l
          enddo
          if (idsmsp(n) .eq. 0) then
            write(iout,'(/,2a)') 'Surface Model species cannot be ',
     &                           'found in master species list: '
            write(iout,'(a)') smspc(n)
            goto 910
          endif
        enddo
        write(idiag,'(/,a)') 'SURFACE MODEL PARAMETERS'
        write(idiag,'(2a)')'     Species     SoilSorb SoilLeach   ',
     &                     'VegSorb    VegPen'
        do n = 1,nsmspc
          write(idiag,'(i3,2x,a10,6(1pe10.2))') idsmsp(n),smspc(n),
     &                         smssrb(n),smlch(n),smvsrb(n),smpen(n),
     &                         smisrb(n),smmlt(n)
        enddo
c
c-----Reaction list and associated parameters
c
        read(ichem,'(a)',err=905,end=906)  record
        do n = 1,nsmrxn
          read(ichem,'(a)',err=905,end=906) record
          read(record,'(5x,2a10,6e10.0)',err=905) smpre(n),smprd(n),
     &                                         smskrat(n),smsjrat(n),
     &                                         smvkrat(n),smvjrat(n),
     &                                         smikrat(n),smijrat(n)
          idsmpre(n) = 0
          idsmprd(n) = 0
          do l = 1,nsmspc
            if (smpre(n).eq.smspc(l)) idsmpre(n) = l
            if (smprd(n).eq.smspc(l)) idsmprd(n) = l
          enddo
          if (idsmpre(n) .eq. 0) then
            write(iout,'(/,2a)')'Reaction precursor species cannot be ',
     &                          'found in Surface Model species list: '
            write(iout,'(a)') smpre(n)
            goto 910
          endif
          if (idsmprd(n) .eq. 0) then
            write(iout,'(/,2a)')'Reaction product species cannot be ',
     &                          'found in Surface Model species list: '
            write(iout,'(a)') smprd(n)
            if (smprd(n).eq.'          ') then
              write(iout,'(a)') 'Reaction product will be ignored'
            else
              goto 910
            endif
          endif
        enddo
        write(idiag,'(2a)')'     Precursor Product        Krate     ',
     &                     'Jrate'
        do n = 1,nsmrxn
          write(idiag,'(5x,2a10,6(1pe10.2))') smpre(n),smprd(n),
     &                                     smskrat(n),smsjrat(n),
     &                                     smvkrat(n),smvjrat(n),
     &                                     smikrat(n),smijrat(n)
        enddo
      endif
c
      write(idiag,*)
      call flush(idiag)
c
      return
c
 900  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(4(a,i4))')
     &  'reaction number', num, ' of type', rxntyp(i),
     &  ' should have', npar(abs(rxntyp(i))), ' parameters'
      goto 910
c
 901  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)')  'For reaction type 5'
      write(iout,'(a,i4)') 'reference reaction # ', iref
      write(iout,'(a,i4)') 'must come before this reaction', i
      goto 910
c
 902  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)') 'Reaction type out of bounds (1-7)'
      write(iout,'(a,i4)') 'for reaction ', num
      write(iout,'(a,i4)') 'value was', rxntyp(i)
      goto 910
c
 903  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(/,a,/)') record
      write(iout,'(2a)') 'If you include aerosol species you must',
     &                   ' specify: number, coupling time step, '
      write(iout,'(a)')  'number of sizes, and size intervals'
      goto 910
c
 904  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file'
      write(iout,'(a)') 'Expecting record label: Surface Model.'
      write(iout,'(a)') 'You have selected the Surface Model treatment:'
      write(iout,'(3a)') 'Either provide a chemistry paramters file ',
     &                   'with Surface Model records or turn off ',
     &                   'Surface Model.'
      goto 910
c
 905  write(iout,'(//,2a)') 'ERROR: Reading CHEMPARAM Surface Model',
     &                      ' record'
      write(iout,'(/,a,/)') record
      write(iout,'(a)') 'Review your Surface Model input data'
      goto 910
c
 906  write(iout,'(//,2a)') 'ERROR: Reading CHEMPARAM Surface Model',
     &                     ' record'
      write(iout,'(a)') 'Unexpected end of file'
      write(iout,'(a)') 'Review your Surface Model input data'
      goto 910
c
 907  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file'
      write(iout,'(2a)')   'Species: ',spname(l)
      write(iout,'(2a)')   'Exceeds the 8-character maximum.',
     &                     ' Shorten the name.'
      goto 910
c
 910  write(*,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(idiag,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(iout,'(//,a)') 'ERROR in READCHM.'
      call camxerr()
c
      end
