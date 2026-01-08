      function chk_csa_cell(icl_mod,jcl_mod,kcl_mod,
     &                     ncol_mod,nrow_mod,nlay_mod,ipa_cel)
      implicit none
      integer chk_csa_cell
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c         Calls all of the routines that perform the averaging for 
c         the various CSA concentation fields. This will do all of
c         CSA dubdomains.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c       icl_mod   I-cell in th current modeling grid
c       jcl_mod   J-cell in th current modeling grid
c       kcl_mod   layer in th current modeling grid
c       ncol_mod  number of columns in modling grid
c       nrow_mod  number of columns in modling grid
c       nlay_mod  number of columns in modling grid
c       ipa_cel   array of indexes into CSA grid
c     Output:  
c        Return value is index into CSA arrays (-9 = not in a grid)
c
c
c     Copyright 1996 - 2025
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument delarations:
c-----------------------------------------------------------------------
c
      integer icl_mod
      integer jcl_mod
      integer kcl_mod
      integer ncol_mod
      integer nrow_mod
      integer nlay_mod
      integer ipa_cel(ncol_mod,nrow_mod,nlay_mod)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer idx_csa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      chk_csa_cell = -9
      idx_csa = ipa_cel(icl_mod,jcl_mod,kcl_mod)
      if( idx_csa .GT. 0 ) chk_csa_cell = idx_csa 
      return
      end
