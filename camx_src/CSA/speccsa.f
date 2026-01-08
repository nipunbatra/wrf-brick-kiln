c**** SPECCSA
c
      subroutine speccsa( )
      use filunit
      use chmstry
      use tracer
      use grid
      use procan
c
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets up the species names and pointers into the species
c   for all of the CSA species.  Pointers will be set up for both the
c   concentration array and the emissions array.
c
c     Copyright 1996 - 2025
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/23/99   --gwilson--    Original development (taken from SPECDDM)
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'camx_aero.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 caffec, cinflu, srcnam
      integer      mvec4d, ispc, iaffec, icsa, mvecedge, i, l, idxcsa
      logical      lout, lsns
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- echo the table of species to the output file ---
c
      write(idiag,*) ' '
      write(idiag,*) ' Affected   Influencing   Source',
     &                   '                       Long            Short' 
      write(idiag,*) ' Species      Species      Type ',
     &                   '      Group   Region   Name            Name' 
      write(idiag,'(1X,79A)') ('-',i=1,79)
c
c  --- calculate the number of CSA species needed and allocate
c      some arrays ---
c
      nddmsp = ntermddm + nrateddm
      ntotsp = nspec * nddmsp

      if ( nddmsp .GT. MXFDDM ) goto 7000
      if ( ntotsp .GT. MXTRSP ) goto 7001
      
      lsns = .NOT. lmpi
      mvecedge = MAX(ncol(1),nrow(1))
      call alloc_ddm(lsns,ngrid,ncol,nrow,nlay,nlayers_ems,nspec,
     &                                    nddmsp,nhddm,mvecedge,iout)
      call alloc_csa_conc(npa_cels,npadom,ntotsp,
     &            i_sw,j_sw,b_lay,i_ne,j_ne,t_lay,iptrcsa)
c
c  --- loop over the modeled species and setup the names ---
c
      do ispc=1,nspec
c
c  --- set the pointer into the array for this family ---
c
          iptddm(ispc) = (ispc-1)*nddmsp + 1
          iaffec = ispc
          caffec = spname(ispc)
c
c  --- set the flag for determining if this species should
c      be output to average file ---
c
          lout = .FALSE.
          do i=1,navspc
            if( lavmap(i) .EQ. ispc ) lout = .TRUE.
          enddo
c
c  --- intialize the index into the array ---
c
          idxcsa = iptddm(ispc) - 1
c
c  --- loop over all of the rate constant sensitivity groups ---
c
          do icsa = 1,nrateddm
             cinflu = rateddm(icsa)
             srcnam = 'RATE'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxcsa,icsa,0,lout,0.0)
          enddo
c
c  --- loop over all of the rate term sensitivity groups ---
c
          do icsa = 1,ntermddm
             cinflu = termddm(icsa)
             srcnam = 'TERM'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxcsa,icsa,0,lout,0.0)
          enddo
c
c  --- get then next affect species ---
c
      enddo
      ntotsp = idxcsa
      ipttim = ntotsp + 1
      nsaspc = ntotsp
      write(idiag,'(1X,79A)') ('-',i=1,79)
      write(idiag,*) 
c
c  --- set the flag for gaseous species ---
c
      lsagas = .FALSE.
      do i=iptddm(nrad+1),iptddm(ngas)+nddmsp-1
         lsagas(i) = .TRUE.
      enddo
c
c  --- initialize all of the tracers concs to zero to start off ---
c
      csaconc = 0.
      csaavrg = 0.
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in SPECCSA:'
      write(iout,'(/,1X,A,I5)')
     &            'ERROR: Number of CSA families exceeds max: ',MXFDDM
      write(iout,'(1X,A,I5)')
     &        'You need room for at least this many families: ',nddmsp
      write(iout,'(1X,A)') 'Increase parameter MXFDDM in camx.prm'
      write(iout,'(1X,2A,I5)') 'Also increase parameter MXTRSP',
     &                     ' in camx.prm if it is less than ',ntotsp
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SPECCSA:'
      write(iout,'(/,1X,A,I5)')
     &                 'Number of tracer species exceeds max: ',MXTRSP
      write(iout,'(1X,A,I5)')
     &         'You need room for at least this many tracers: ',ntotsp
      write(iout,'(1X,A)') 'Increase parameter MXTRSP in camx.prm'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
