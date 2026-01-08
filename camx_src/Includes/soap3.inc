**** SOAP3
c
c----CAMx v7.32 250801
c
c
c----------------------------------------------------------------------
c
c    Include file for SOAP3 parameters and variables.
c
c      Copyright 1996 - 2025
c     Ramboll Environ
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     01/03/24   -gwilson-  original development
c
c-----------------------------------------------------------------------
c    Parameters for SOAP model species;
c-----------------------------------------------------------------------
c
c   NUM_SOAP_SPECIES  I  number of species in SOAP scheme
c
      integer NUM_SOAP_SPECIES
c
      parameter( NUM_SOAP_SPECIES = 11 )
c
c-----------------------------------------------------------------------
c    Variables for species names:
c-----------------------------------------------------------------------
c
c   soap2_spec_names  C  species names for SOAP2 inputs
c   soap3_spec_names  C  species names for SOAP3 inputs
c
      character*10 soap2_spec_names(NUM_SOAP_SPECIES)
      character*10 soap3_spec_names(NUM_SOAP_SPECIES)
c
      common /soap_spec/ soap2_spec_names, soap3_spec_names
c
c-----------------------------------------------------------------------
c    Parameters for SOAP emissions species:
c-----------------------------------------------------------------------
c
c  NUM_SOAP_EMISS    I number of emissions species for SOAP
c  IDX_SOAP_IVOC     I index of IVOC in SOAP emissions names array
c  IDX_SOAP_IVOD     I index of IVOD in SOAP emissions names array
c  IDX_SOAP_IVOA     I index of IVOA in SOAP emissions names array
c  IDX_SOAP_IVOB     I index of IVOB in SOAP emissions names array
c  IDX_SOAP_POA_OP   I index of POA_OP in SOAP emissions names array
c  IDX_SOAP_POA_GV   I index of POA_GV in SOAP emissions names array
c  IDX_SOAP_POA_DV   I index of POA_DV in SOAP emissions names array
c  IDX_SOAP_POA_MC   I index of POA_MC in SOAP emissions names array
c  IDX_SOAP_POA_IC   I index of POA_IC in SOAP emissions names array
c  IDX_SOAP_POA_BB   I index of POA_BB in SOAP emissions names array
c
      integer NUM_SOAP_EMISS
      integer IDX_SOAP_IVOC
      integer IDX_SOAP_IVOD
      integer IDX_SOAP_IVOA
      integer IDX_SOAP_IVOB
      integer IDX_SOAP_POA_OP
      integer IDX_SOAP_POA_GV
      integer IDX_SOAP_POA_DV
      integer IDX_SOAP_POA_MC
      integer IDX_SOAP_POA_IC
      integer IDX_SOAP_POA_BB
      integer IDX_SOAP_POA_AV
c
      parameter( NUM_SOAP_EMISS  = 11 )
      parameter( IDX_SOAP_IVOC   =  1 )
      parameter( IDX_SOAP_IVOD   =  2 )
      parameter( IDX_SOAP_IVOA   =  3 )
      parameter( IDX_SOAP_IVOB   =  4 )
      parameter( IDX_SOAP_POA_OP =  5 )
      parameter( IDX_SOAP_POA_GV =  6 )
      parameter( IDX_SOAP_POA_DV =  7 )
      parameter( IDX_SOAP_POA_MC =  8 )
      parameter( IDX_SOAP_POA_IC =  9 )
      parameter( IDX_SOAP_POA_AV = 10 )
      parameter( IDX_SOAP_POA_BB = 11 )
c
c-----------------------------------------------------------------------
c    Variables for species names:
c-----------------------------------------------------------------------
c
c   soap_emiss_names  C  species names for SOAP emissions species
c
      character*10 soap_emiss_names(NUM_SOAP_EMISS)
c
      common /soap_emiss/ soap_emiss_names
c
c-----------------------------------------------------------------------
c    Variables for evaporation factors:
c-----------------------------------------------------------------------
c
c   soap3_evap_factors R coefficients of equations for evap functions
c
      real soap3_evap_factors(3,NUM_SOAP_EMISS)
      common /soap_evap/ soap3_evap_factors
c
      real, parameter ::
     &      soap3_evap_factors_15dbs(3,IDX_SOAP_POA_OP:IDX_SOAP_POA_BB) = reshape(
     &           (/  3.50E-05,    !   IDX_SOAP_POA_OP
     &              -2.98E-02,    !   IDX_SOAP_POA_OP
     &               6.21E+00,    !   IDX_SOAP_POA_OP
     &              -4.09E-05,    !   IDX_SOAP_POA_GV
     &               1.75E-02,    !   IDX_SOAP_POA_GV
     &              -8.47E-01,    !   IDX_SOAP_POA_GV
     &              -8.97E-05,    !   IDX_SOAP_POA_DV
     &               4.42E-02,    !   IDX_SOAP_POA_DV
     &              -4.46E+00,    !   IDX_SOAP_POA_DV
     &              -2.60E-05,    !   IDX_SOAP_POA_MC
     &               1.13E-02,    !   IDX_SOAP_POA_MC
     &              -2.09E-01,    !   IDX_SOAP_POA_MC
     &              -6.53E-05,    !   IDX_SOAP_POA_IC
     &               3.09E-02,    !   IDX_SOAP_POA_IC
     &              -2.65E+00,    !   IDX_SOAP_POA_IC
     &              -1.60E-05,    !   IDX_SOAP_POA_AV
     &               3.11E-03,    !   IDX_SOAP_POA_AV
     &               1.21E+00,    !   IDX_SOAP_POA_AV
     &              -1.89E-05,    !   IDX_SOAP_POA_BB
     &               3.99E-04,    !   IDX_SOAP_POA_BB
     &               2.11E+00 /), !   IDX_SOAP_POA_BB
     &           (/ 3, IDX_SOAP_POA_BB-IDX_SOAP_POA_OP+1 /) )
c
      real, parameter :: 
     &      soap3_evap_factors_complx(3,IDX_SOAP_POA_OP:IDX_SOAP_POA_BB) = reshape(
     &           (/  3.71E-05,    !   IDX_SOAP_POA_OP
     &              -2.88E-02,    !   IDX_SOAP_POA_OP
     &               5.86E+00,    !   IDX_SOAP_POA_OP
     &              -3.78E-05,    !   IDX_SOAP_POA_GV
     &               1.61E-02,    !   IDX_SOAP_POA_GV
     &              -6.86E-01,    !   IDX_SOAP_POA_GV
     &              -7.85E-05,    !   IDX_SOAP_POA_DV
     &               3.79E-02,    !   IDX_SOAP_POA_DV
     &              -3.57E+00,    !   IDX_SOAP_POA_DV
     &              -2.50E-05,    !   IDX_SOAP_POA_MC
     &               1.09E-02,    !   IDX_SOAP_POA_MC
     &              -1.72E-01,    !   IDX_SOAP_POA_MC
     &              -5.82E-05,    !   IDX_SOAP_POA_IC
     &               2.70E-02,    !   IDX_SOAP_POA_IC
     &              -2.13E+00,    !   IDX_SOAP_POA_IC
     &              -1.60E-05,    !   IDX_SOAP_POA_AV
     &               3.11E-03,    !   IDX_SOAP_POA_AV
     &               1.21E+00,    !   IDX_SOAP_POA_AV
     &               2.41E-05,    !   IDX_SOAP_POA_BB
     &              -2.05E-02,    !   IDX_SOAP_POA_BB
     &               4.60E+00 /), !   IDX_SOAP_POA_BB
     &           (/ 3, IDX_SOAP_POA_BB-IDX_SOAP_POA_OP+1 /) )
c
c-----------------------------------------------------------------------
c    Variables for extrac classes for SOAP3 emissions:
c-----------------------------------------------------------------------
c
c   soap3_ntrcls I number of classes to use, including SOAP3 emissions classes
c   soap3_clsnam C name of each tracer class for SOAP3 emissions species
c   soap3_iptcls I pointer into tracer arrays for SOAP3 emissions class
c   soap3_nptcls I number of tracer species in each class
c   soap3_iemcls I (redundant) pointer into tracer arrays for SOAP3 emissions class
c
      integer     soap3_ntrcls
      character*3 soap3_clsnam(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      integer     soap3_iptcls(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      integer     soap3_nptcls(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      integer     soap3_iemcls(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
c
      common /soap3_ptr_chr/ soap3_clsnam
      common /soap3_ptr_dat/ soap3_ntrcls, soap3_iptcls, soap3_nptcls, soap3_iemcls
