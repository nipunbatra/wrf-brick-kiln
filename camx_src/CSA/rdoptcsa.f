c**** RDOPTCSA
c
      subroutine rdoptcsa()
      use filunit
      use grid
      use chmstry
      use procan
      use tracer
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     This routine loads all of the user options and flags for the
c     CSA algorithm.  
c
c     Copyright 1996 - 2025
c     Ramboll
c
c     Argument description:
c           none
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     04/26/19   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'namelist.inc'
c
c-----------------------------------------------------------------------
c    Local Variables:
c-----------------------------------------------------------------------
c
      character*200 strtmp
      character*10  tmpnam
      integer i, j, k1, k2, k3, n, igrd, nox, noy
      integer numrxn,numreac,numprod,idxrxn,idxspc
      logical lerror
c
      integer, allocatable :: reactmp(:,:), prodtmp(:,:)
      real,    allocatable :: coeftmp(:,:)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. lcsa ) goto 9999
c
c  --- check that switches are compatible ---
c
      if( naero .GT. 0 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTCSA:'
         write(iout,*) 'Cannot use CSA with Aerosol chemistry.'
         write(iout,*)
     &   '  Turn off CSA or use a different chemical mechanism.'
         call camxerr()
      endif
c
c   --- number of rate constant sensitivity groups ---
c
      nrateddm = Number_of_Rate_Const_Groups
      allocate( rateddm(nrateddm) )
      allocate( iprate(0:nreact, nrateddm ) )
      iprate = 0
      do i = 1, nrateddm
        strtmp = ADJUSTL(Rate_Const_Groups(i))
        j = INDEX(strtmp,':')
        if (j.eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
          write(iout,'(1X,A,I3,A)')
     &        'Delimiter (:) is not found in Rate_Const_Groups(',i,').'
          call camxerr()
        endif
        rateddm(i) = strtmp(:j-1)
        do
          strtmp = strtmp(j+1:)
          j = INDEX(strtmp,',')
          if (j.eq.0) then
            j = LEN_TRIM(strtmp) + 1
            if (j.eq.1) EXIT
          endif
          iprate(0,i) = iprate(0,i) + 1
          if (iprate(0,i).gt.nreact) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &            'Too many rxns (>nreact) in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
          if ( VERIFY( TRIM( ADJUSTL( strtmp(:j-1) ) ), '1234567890' )
     &                                                    .ne. 0 ) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &                 'Invalid rxn index in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
          read(strtmp(:j-1),'(I5)') iprate(iprate(0,i),i)
          if ( iprate(iprate(0,i),i).lt.1 .or.
     &         iprate(iprate(0,i),i).gt.nreact ) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &            'Out-of-range rxn index in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
        enddo
        if (iprate(0,i).eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
          write(iout,'(1X,A,I3,A)')
     &                   'No rxn is found in Rate_Const_Groups(',i,').'
          call camxerr()
        endif
      enddo
c
c   --- number of rate term sensitivity groups ---
c
      ntermddm = Number_of_Rate_Term_Groups
      allocate( termddm(ntermddm) )
      termddm = ' '
      allocate( ipterm(0:ngas, nreact, ntermddm ) )
      ipterm = 0
      allocate( wfterm(  ngas, nreact, ntermddm ) )
      wfterm = 1.0
      allocate( reactmp(ngas,nreact) )
      allocate( prodtmp(ngas,nreact) )
      allocate( coeftmp(ngas,nreact) )
      do i = 1, ntermddm
        numrxn = 0
        reactmp = 0
        prodtmp = 0
        coeftmp = 1.0
        do n = 1, MXNAM
          strtmp = ADJUSTL(Rate_Term_Groups(i,n))
          if ( LEN_TRIM(strtmp).eq.0 ) CYCLE
          j = INDEX(strtmp,':')
          if (j.eq.0) then
            write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
            write(iout,'(1X,A,I3,A)')
     &          'Delimiter (:) is not found in Rate_Term_Groups(',i,').'
            call camxerr()
          endif
          if ( LEN_TRIM(termddm(i)).eq.0 ) then
            termddm(i) = strtmp(:j-1)
          else
            if ( termddm(i).ne.strtmp(:j-1) ) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A)')
     &                   'Inconsistent name of Rate_Term_Groups(',i,').'
              call camxerr()
            endif
          endif
          do
            strtmp = strtmp(j+1:)
            j = INDEX(strtmp,',')
            if (j.eq.0) then
              j = LEN_TRIM(strtmp) + 1
              if (j.eq.1) EXIT
            endif
            k1 = INDEX(strtmp(:j-1),'[R')
            k2 = INDEX(strtmp(:j-1),'[P')
            k3 = INDEX(strtmp(:j-1),']')
            if (k1+k2.eq.0 .or. k3.eq.0 .or. k3.lt.k1+k2+2) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A)')
     &                     'Invalid keyword in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            if ( VERIFY( TRIM( ADJUSTL( strtmp(:k1+k2-1) ) ),
     &                                     '1234567890' ) .ne. 0 ) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A)')
     &                   'Invalid rxn index in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            read(strtmp(:k1+k2-1),*) idxrxn
            if ( idxrxn.lt.1 .or. idxrxn.gt.nreact ) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A)')
     &              'Out-of-range rxn index in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            numrxn = numrxn + 1
            ipterm(0,idxrxn,i) = 1
            tmpnam = strtmp(k3+1:j-1)
            do idxspc = 1, ngas
              if (tmpnam.eq.spname(idxspc)) EXIT
            enddo
            if (idxspc.gt.ngas) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A,I3,A)')
     &          'Invalid spc name in rxn(',idxrxn,
     &                                   ') of Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            if (k1.gt.0) then ! Reactant
              reactmp(idxspc,idxrxn) = 1
            else              ! Product
              prodtmp(idxspc,idxrxn) = 1
            endif
            if (strtmp(k1+k2+2:k1+k2+2).eq.'#') then
              if ( k2.ne.0 ) then ! Check if a weighting factor is assigned to a product; may be allowed later
                write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
                write(iout,'(1X,A)')
     &                        'Products must not have weighting factors'
                call camxerr()
              endif
              if ( VERIFY( TRIM( ADJUSTL( strtmp(k1+k2+3:k3-1) ) ),
     &                                '1234567890.+-Ee' ) .ne. 0 ) then
                write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
                write(iout,'(1X,A,I3,A)')
     &                 'Invalid coefficient in Rate_Term_Groups(',i,').'
                call camxerr()
              endif
              read(strtmp(k1+k2+3:k3-1),*) coeftmp(idxspc,idxrxn)
            endif
          enddo
        enddo ! loop n
        if (numrxn.eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
          write(iout,'(1X,A,I3,A)')
     &                     'No rxn is found in Rate_Term_Groups(',i,').'
          call camxerr()
        endif
        do idxrxn = 1, nreact
          numreac = 0
          do idxspc = 1, ngas
            if ( reactmp(idxspc,idxrxn).eq.0 ) CYCLE
            numreac = numreac + 1
            ipterm(1+numreac,idxrxn,i) = idxspc
            wfterm(1+numreac,idxrxn,i) = coeftmp(idxspc,idxrxn)
          enddo
          ipterm(1,idxrxn,i) = numreac
          numprod = 0
          do idxspc = 1, ngas
            if ( prodtmp(idxspc,idxrxn).eq.0 ) CYCLE
            numprod = numprod + 1
            if ( 2+numreac+numprod.gt.ngas ) then
              write(iout,'(//,A)') 'ERROR in RDOPTCSA:'
              write(iout,'(1X,A,I3,A,I3,A)')
     &          'Too many spc selected in rxn(',idxrxn,
     &                                  ') of Rate_Const_Groups(',i,').'
              call camxerr()
            endif
            ipterm(2+numreac+numprod,idxrxn,i) = idxspc
          enddo
          ipterm(2+numreac,idxrxn,i) = numprod
        enddo ! loop idxrxn
      enddo ! loop i
      deallocate( reactmp )
      deallocate( prodtmp )
      deallocate( coeftmp )
c
c   ---- Number of Process Analysis domains
c
      npadom = Number_of_CSA_Domains
c
c   ---- Allocate the arrays that are dimensioned by grids ---
c
      if( npadom .LE. 0 ) goto 7000
      call alloc_procan(ngrid,ncol,nrow,nlay,.TRUE.)
c
c   ---- CSA domain definitions ---
c
      do n = 1,npadom
c
c   --- the grid number ---
c
         ipagrd(n) = Within_CAMx_Grid(n)
         if( ipagrd(n) .LE. 0 .OR. ipagrd(n) .GT. ngrid ) goto 7001
c
c   --- the cell indexes in the X and Y directions
c
         i_sw(n) = CSA_Beg_I_Index(n)
         i_ne(n) = CSA_End_I_Index(n)
         j_sw(n) = CSA_Beg_J_Index(n)
         j_ne(n) = CSA_End_J_Index(n)
         b_lay(n) = CSA_Beg_K_Index(n)
         t_lay(n) = CSA_End_K_Index(n)
c
c   --- check for valid cell indexes ---
c
         lerror = .FALSE.
         igrd = ipagrd(n)
         if( igrd .GT. 1 ) then
            nox = (inst2(igrd) - inst1(igrd) + 1 ) * meshold(igrd) + 2
            noy = (jnst2(igrd) - jnst1(igrd) + 1 ) * meshold(igrd) + 2
         else
            nox = ncol(1)
            noy = nrow(1)
         endif
         if( i_sw(n) .LE. 0 .OR. i_sw(n) .GT. i_ne(n) ) lerror = .TRUE. 
         if( i_ne(n) .GT. nox ) lerror = .TRUE. 
         if( j_sw(n) .LE. 0 .OR. j_sw(n) .GT. j_ne(n) ) lerror = .TRUE. 
         if( j_ne(n) .GT. noy ) lerror = .TRUE. 
         if( b_lay(n) .LE. 0 .OR. b_lay(n) .GT. t_lay(n) ) 
     &                                                  lerror = .TRUE.
         if( t_lay(n) .GT. nlay(igrd) ) lerror = .TRUE.
         if( lerror ) goto 7002
      enddo
c
c   --- dummy array for tracer wetdep field ---
c
      allocate( ptwetfld(1) )
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
      write(iout,'(//,a)') 'ERROR in RDOPTCSA:'
      write(iout,'(1X,2A,I3)') 'Invalid value for number of ',
     &                   'CSA domains: ',npadom
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDOPTCSA:'
      write(iout,'(1X,2A,I3)') 'Invalid value for grid number for ',
     &                         'CSA domain: ',ipagrd(n)
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDOPTCSA:'
      write(iout,'(1X,A,I3,A)') 'Cell indexes for CSA sub-domain: ',n,
     &                                                 ' are invalid.'
      write(iout,'(10X,A,I3,10X,A)') 'Sub-Domain #',n,'Grid Definition'
      write(iout,'(14X,2A,18X,A,I3)') 'Beg ',' End','Grid #',igrd
      write(iout,'(A10,1X,2I5,20X,I5)') 'Rows    :',i_sw(n),i_ne(n),nox
      write(iout,'(A10,1X,2I5,20X,I5)') 'Columns :',j_sw(n),j_ne(n),noy
      write(iout,'(A10,1X,2I5,20X,I5)') 'Layers  :',
     &                                    b_lay(n),t_lay(n),nlay(igrd)
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
