c*** CSA_RESET
c
      subroutine csa_reset(igrid_mod,ncol_mod,nrow_mod,nlay_mod,
     &               nspc_mod,ngas_mod,nrad_mod,nhet,ipa_cel)
      use grid
      use camxfld
      use procan
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c   Description:
c      Reloads the regular model concentratons into the CSA gridded
c      array. This will happen each hour.
c
c     Copyright 1996 - 2025
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c         igrid_mod  number of columns in modeling grid
c         ncol_mod   number of columns in modeling grid
c         nrow_mod   number of rows in modeling grid
c         nlay_mod   number of layers in modeling grid
c         nspc_mod   total number of modelled species
c         ngas_mod   number of modeled gas species
c         nrad_mod   number of modeled radicals
c         nhet       number of heterogeneous rxns
c         conc_mod   gridded array of model concentrations for this grid
c         ipa_cel    index of model cell into CSA grid
c     Outputs:
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
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
      integer igrid_mod
      integer ncol_mod
      integer nrow_mod
      integer nlay_mod
      integer nspc_mod
      integer ngas_mod
      integer nrad_mod
      integer nhet
      integer ipa_cel(ncol_mod,nrow_mod,nlay_mod)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer icl_mod, jcl_mod, kcl_mod, ipa_idx, icsa_grid
      integer icsa_icl, icsa_jcl, icsa_kcl, ncol_csa, nrow_csa, nlay_csa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- loop over modelling cells ---
c
      do kcl_mod=1,nlay_mod
          do jcl_mod=1,nrow_mod
              do icl_mod=1,ncol_mod
c
c  --- skip if this cell is not in a CSA grid ---
c
                 ipa_idx = ipa_cel(icl_mod,jcl_mod,kcl_mod)
                 if( ipa_idx .LE. 0 ) cycle
c
c  --- get CSA grid values for this cell ---
c
                 icsa_grid = ipadom(ipa_idx)
                 icsa_icl = ipax(ipa_idx) - i_sw(icsa_grid) + 1
                 icsa_jcl = ipay(ipa_idx) - j_sw(icsa_grid) + 1
                 icsa_kcl = ipaz(ipa_idx) - b_lay(icsa_grid) + 1
                 ncol_csa = i_ne(icsa_grid) - i_sw(icsa_grid) + 1
                 nrow_csa = j_ne(icsa_grid) - j_sw(icsa_grid) + 1
                 nlay_csa = t_lay(icsa_grid) - b_lay(icsa_grid) + 1
c
c  --- call routine to load the model concs into CSA array ---
c
                 call loadcsa_conc(ncol_csa,nrow_csa,
     &             nlay_csa,ntotsp,csaconc(iptrcsa(icsa_grid)),
     &             icsa_icl,icsa_jcl,icsa_kcl,icl_mod,jcl_mod,
     &             kcl_mod,nddmsp,nspc_mod,ngas_mod,nrad_mod,nhet,
     &             ncol_mod,nrow_mod,nlay_mod,conc(iptr4d(igrid_mod)))
              enddo
          enddo
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
