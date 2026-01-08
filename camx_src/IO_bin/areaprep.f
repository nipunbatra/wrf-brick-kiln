      subroutine areaprep(igrid,idxfile,begtim,begdate,endtim,enddate,iounit,
     &                    iout_iounit,idiag_iounit,dxmod,dymod)
      use grid
      use chmstry
      use filunit
c 
c----CAMx v7.32 250801
c 
c     AREAPREP reads the header of binary area source emissions file,
c     and maps the area source species list to the internal CAMx species list
c                           
c      Copyright 1996 - 2025
c     Ramboll
c           
c     Modifications: 
c        1/20/99   Grid cell size from file should be meters for all cartesian
c                  projections (UTM, LCP, PSP)
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c        11/06/01  Input dates are now Julian
c        02/09/02  Added code to handle end of year dates
c        12/07/14  Revised for VBS emissions
c        01/08/16  Updated for Revised VBS emissions
c        05/13/16  Added checks for I2/HOI with in-line Ix emissions
c        11/22/17  Added checks for PFE/PMN/PK/PCA/PMG emissions
c 
c     Input arguments: 
c        igrid               grid index
c        idxfile             file index
c        begtim              model start time (HHMM) 
c        begdate             model start date (YYJJJ) 
c        endtim              model end time (HHMM) 
c        enddate             model end date (YYJJJ)
c        iounit              area emissions file unit
c        iout_iounit         output message file unit
c        idiag_iounit        output diagnostic file unit
c        dxmod               model grid size in x-direction (deg or km) 
c        dymod               model grid size in y-direction (deg or km)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        STARTUP 
c 
      include 'camx.prm'
      include 'flags.inc'
      include 'vbs.inc'
      include 'soap.inc'
      include 'soap3.inc'
c
      character*4 ifile(10),note(60)
      integer begdate,enddate
      character*10 arfil,infil,arspc 
c
      character*4 arspec(10,MXSPEC)
c
      integer, parameter :: num_elem = 5
      integer :: idx_elem(num_elem)
      integer    idxar, idx_svoa, idx_hoa
      logical, allocatable, dimension(:) :: lfound
      logical lfound_poa(IDX_SOAP_POA_OP:IDX_SOAP_POA_BB)
      logical lfound_svoa, lfound_hoa
c
      integer isoa
c
      data arfil /'EMISSIONS '/
c
c-----Entry point
c
c-----Read 1st AREA header record and check inputs 
c             
      idx_svoa = -9
      idx_hoa = -9
      rewind(iounit)
      read(iounit) ifile,note,nseg,narspc(igrid,idxfile),idat1,tim1,idat2,tim2
      allocate( lfound(narspc(igrid,idxfile)) )
c             
      if( INT(tim2) .EQ. 24 ) then
          idat2 = idat2 + 1
          tim2 = 0.
          if( MOD(idat2,1000) .GT. 365 ) then 
             if( MOD(INT(idat2/1000),4) .EQ. 0 ) then 
                if( MOD(idat2,1000) .EQ. 367 )
     &                    idat2 = (INT(idat2/1000)+1)*1000 + 1
             else
                idat2 = (INT(idat2/1000)+1)*1000 + 1 
             endif
          endif
      endif
      write(infil,'(10a1)') (ifile(n),n=1,10) 
      if (infil.ne.arfil) then 
        write(iout_iounit,'(//,a)') 'ERROR in AREAPREP:'
        write(iout_iounit,*)'AREA input file is not labelled EMISSIONS' 
        call camxerr()
      endif   
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if (idat1.gt.begdate) then 
        write(iout_iounit,'(//,a)') 'WARNING in AREAPREP:'
        write(iout_iounit,*)'AREA start date > simulation start date' 
        write(iout_iounit,*)'AREA file: ',idat1
        write(iout_iounit,*)'Sim start: ',begdate 
        if (.not.le1day) then
          write(iout_iounit,*)'CAMx expecting day-specific emissions: Stopping'
          call camxerr()
        endif
      elseif (idat1.eq.begdate .and. tim1.gt.begtim) then 
        write(iout_iounit,'(//,a)') 'ERROR in AREAPREP:'
        write(iout_iounit,*)'AREA start time > simulation start time' 
        write(iout_iounit,*)'AREA file: ',tim1
        write(iout_iounit,*)'Sim start: ',begtim 
        call camxerr()
      elseif (idat2.lt.enddate) then 
        write(iout_iounit,'(//,a)') 'WARNING in AREAPREP:'
        write(iout_iounit,*)'AREA end date < simulation end date' 
        write(iout_iounit,*)'AREA file: ',idat2
        write(iout_iounit,*)'  Sim end: ',enddate 
        if (.not.le1day) then
          write(iout_iounit,*)'CAMx expecting day-specific emissions: Stopping'
          call camxerr()
        endif
      elseif (idat2.eq.enddate .and. tim2.lt.endtim) then 
        write(iout_iounit,'(//,a)') 'WARNING in AREAPREP:'
        write(iout_iounit,*)'AREA end time < simulation end time' 
        write(iout_iounit,*)'AREA file: ',tim2
        write(iout_iounit,*)'  Sim end: ',endtim 
        call camxerr()
      endif 
c 
c-----Read 2nd AREA header record and check inputs 
c 
      read(iounit) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz 
      if (.NOT.llatlon) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if (abs(dx-dxmod).gt.0.001 .or. abs(dy-dymod).gt.0.001) then
        write(iout_iounit,'(//,a)') 'WARNING in AREAPREP:'
        write(iout_iounit,*)'AREA cell size not equal to model cell size'
        write(iout_iounit,'(a,2f10.4)')'AREA file: ',dx,dy
        write(iout_iounit,'(a,2f10.4)')'    model: ',dxmod,dymod
        write(iout_iounit,*)
      endif
      if (nx.ne.ncol(igrid) .or. ny.ne.nrow(igrid)) then
           if( nx .EQ. ncol(igrid)-2 .AND. ny .EQ. nrow(igrid)-2 ) then
              buffer_offset_iarem(igrid,idxfile) = 1
           else 
              write(iout_iounit,'(//,a)') 'ERROR in AREAPREP:'
              write(iout_iounit,*)'AREA grid size not equal to model grid size '
              write(iout_iounit,*)'   grid #: ',igrid
              write(iout_iounit,*)'AREA file: ',nx,ny,nz
              write(iout_iounit,*)'    model: ',ncol(igrid),nrow(igrid),nlay(igrid)
              write(iout_iounit,*)
              call camxerr()
           endif
      else
           buffer_offset_iarem(igrid,idxfile) = 0
      endif 
c 
c-----Read 3rd & 4th AREA header 
c 
      read(iounit) 
      read(iounit) ((arspec(n,l),n=1,10),l=1,narspc(igrid,idxfile)) 
c
c----Find this species in master list ---
c
      do l=1,narspc(igrid,idxfile)
        write(arspc,'(10a1)') (arspec(n,l),n=1,10) 
        if( arspc .EQ. 'HNO2      ') arspc = 'HONO      '
        if( arspc .EQ. 'HCHO      ' .AND. kHCHO .EQ. nspec+1)
     &                                       arspc = 'FORM      '
        do isoa=1,NUM_SOAP_SPECIES
             if( arspc(1:4) .EQ. 'POA ' ) cycle
             if( arspc .EQ. soap2_spec_names(isoa) )
     &                            arspc = soap3_spec_names(isoa)
        enddo
        if( arspc(1:4) .EQ. 'HOA ') arspc = 'POA       '
        if( arspc(1:4) .EQ. 'POA ' .AND. luse_soap3 ) arspc = 'POA_AV    '
        idxems_files(igrid,idxfile,l) = 0
c
c  --- find species in list ---
c
        lfound(l) = .FALSE.
        do ispc=1,nspec
           if(arspc .EQ. spname(ispc)) lfound(l) = .TRUE.
           if(spname(ispc)(1:4) .EQ. 'HOA ') idx_hoa = ispc
           if(spname(ispc)(1:5) .EQ. 'SVOA ') idx_svoa = ispc
        enddo
        do ispc=1,NUM_SOAP_EMISS
           if( arspc .eq. soap_emiss_names(ispc) ) lfound(l) = .TRUE.
       enddo
        if( .NOT. lfound(l) ) cycle
        do k=1,nemspc
          if(arspc .EQ. emspcname(k) ) idxems_files(igrid,idxfile,l) = k
        enddo
        if( idxems_files(igrid,idxfile,l) .EQ. 0 ) then
           nemspc = nemspc + 1
           emspcname(nemspc) = arspc
           idxems_files(igrid,idxfile,l) = nemspc
           idxar = nemspc
        endif
      enddo
c 
c-----Map AREA species to model species 
c 
      do 15 lar = 1,narspc(igrid,idxfile)
        if( .NOT. lfound(lar) ) cycle
        write(arspc,'(10a1)') (arspec(n,lar),n=1,10) 
        if (arspc.eq.'HNO2      ') arspc = 'HONO      '
        if (arspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        arspc = 'FORM      '
        do isoa=1,NUM_SOAP_SPECIES
             if( arspc(1:4) .EQ. 'POA ' ) cycle
             if( arspc .EQ. soap2_spec_names(isoa) )
     &                            arspc = soap3_spec_names(isoa)
        enddo
        if( arspc(1:4) .EQ. 'HOA ') arspc = 'POA       '
        if( arspc(1:4) .EQ. 'POA ' .AND. luse_soap3 ) arspc = 'POA_AV    '
c
c   --- update names for SOAP3 ---
c
        do isoa=1,NUM_SOAP_SPECIES
           if( arspc(1:4) .EQ. 'POA ' ) cycle
           if( arspc .EQ. soap2_spec_names(isoa) ) then
               arspc = soap3_spec_names(isoa)                   
           endif
        enddo
        do 20 l=1,nspec 
          
          if(arspc.eq.spname(l)) then 
            lemmap(l) = idxems_files(igrid,idxfile,lar)
            write(idiag_iounit,'(2(a,i5,2x,a))')
     &                   'Area source species ',lar,arspc, 
     &                   ' mapped to model species ',l,spname(l) 
            if (arspc.eq.'I2        ' .OR. arspc.eq.'HOI       ') then
              if (lixemis) then
                write(iout_iounit,'(//,A)') 'ERROR in AREAPREP:'
                write(iout_iounit,'(A)') 'In-line Ix emissions are invoked,'
                write(iout_iounit,'(2A)')'but I2 and/or HOI found',
     &                            ' in an input gridded emissions file.'
                write(iout_iounit,'(A)') 'Either turn off in-line Ix emissions'
                write(iout_iounit,'(2A)')'or remove I2 and HOI emissions from',
     &                            ' the gridded emissions file.'
                call camxerr()
              endif
            endif
            goto 15 
          endif 
c
c  --- check for SOAP3 species ---
c
          if( luse_soap3 ) then
            if( spname(l)(1:4) .NE. 'HOA ' ) cycle
            do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
              if( arspc .EQ. soap_emiss_names(isoa) ) then
                  lemmap(l) = idxems_files(igrid,idxfile,lar)
                  write(idiag_iounit,'(2(a,i5,2x,a))')
     &                   'Area source species ',lar,arspc, 
     &                   ' mapped to model species ',l,spname(l) 
              endif
            enddo
          endif
 20     continue 
        write(idiag_iounit,*)'AREA species: ',arspc,' not modeled'
 15   continue
c
c --- make sure that SVOA is in the list of emissions species ---
c
      if( luse_soap3 .AND. idx_svoa .GT. 0 ) then
          lfound_svoa = .FALSE.
          do i=1,nemspc
            if( emspcname(i)(1:5) .EQ. 'SVOA ' ) then
                lfound_svoa = .TRUE.
                lemmap(idx_svoa) = i
            endif
          enddo
          if( .NOT. lfound_svoa ) then
             nemspc = nemspc + 1
             emspcname(nemspc) = 'SVOA      ' 
             lemmap(idx_svoa) = nemspc
          endif
      endif
      if( luse_soap3 .AND. idx_hoa .GT. 0 ) then
          lfound_hoa = .FALSE.
          do i=1,nemspc
            if( emspcname(i)(1:4) .EQ. 'POA ' ) then
                lfound_hoa = .TRUE.
                lemmap(idx_hoa) = i
            endif
          enddo
          if( .NOT. lfound_hoa ) then
             nemspc = nemspc + 1
             emspcname(nemspc) = 'POA       ' 
             lemmap(idx_hoa) = nemspc
          endif
      endif
      if( luse_soap3 .AND. idx_hoa .GT. 0 ) then
          do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
             do i=1,nemspc
               if( soap_emiss_names(isoa) .EQ. emspcname(i) ) then
                   lfound_poa(isoa) = .TRUE.
               endif
             enddo
             if( .NOT. lfound_poa(isoa) ) then
                nemspc = nemspc + 1
                emspcname(nemspc) = soap_emiss_names(isoa)
             endif
          enddo
      endif
c
      write(idiag_iounit,*)
      deallocate( lfound )
c
 9999 continue
c             
      return
      end
