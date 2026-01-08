      program wrfcamx
!
!-----WRFCAMx reads WRF meteorological output files and 
!     interpolates/aggregates the data to a CAMx grid.
!
!     THIS VERSION FOR WRF/ARW v3.2+ to CAMx v7.0+
!
!     WRF data on an Arakawa C grid:
!                                              
!                 v(i,j+1)                    j
!                ----.----                    |
!               |         |                   |
!               | T,q,p,Kv|                   |
!         u(i,j).    x    .u(i+1,j)           x-------i
!               |  (i,j)  |
!               |         |
!                ----.----
!                  v(i,j)
!
!     CAMx data on an Arakawa C grid:
!                                              
!                  v(i,j)                     j
!                ----.----                    |
!               |         |                   |
!               | T,q,p,Kv|                   | 
!       u(i-1,j).    x    .u(i,j)             x-------i
!               |  (i,j)  |
!               |         |
!                ----.---- 
!                 v(i,j-1)
!
!     NOTE: this version contains 2 domain mapping routines:
!         - WRF_WINDOW windows WRF data to a CAMx grid with the same projection
!           as WRF.
!           THIS CHOICE LEADS TO STAGGERED WINDS
!         - WRF_INTERP horizontally interpolates WRF data to a CAMx grid with a
!           different projection from WRF and rotates winds to the new 
!           projection.
!           THIS CHOICE LEADS TO ALL VARIABLES AT CELL CENTERS
!         Choice of interpolation is selected using the "PROJECT" variable
!
!         The following mappings are allowed:
!           WRF LCP -> CAMx LCP (identical or different)
!                   -> CAMx UTM
!                   -> CAMx LATLON
!           WRF PS  -> CAMx PS (identical)
!           WRF MRC -> CAMx MRC (identical or different)
!                   -> CAMx UTM
!                   -> CAMx LATLON
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Copyright (C) 2006-2022  Ramboll
!c
!c This program is free software; you can redistribute it and/or
!c modify it under the terms of the GNU General Public License
!c as published by the Free Software Foundation; either version 2
!c of the License, or (at your option) any later version.
!c
!c This program is distributed in the hope that it will be useful,
!c but WITHOUT ANY WARRANTY; without even the implied warranty of
!c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!c GNU General Public License for more details.
!c
!c To obtain a copy of the GNU General Public License
!c go to the Free Software Foundation at http://www.fsf.org.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      USE wrf_netcdf
      USE wrf_fields
      USE camx_fields
      USE ncf_iomod

      implicit none
      include 'netcdf.inc'
!
!-----WRF netCDF time variables
!
      character*80 varnam
      character*80, allocatable, dimension(:) :: times
      integer dimids(10)
      integer ivtype
      integer cdfid, rcode, id_var
      integer ndims, natts, dims(4)
      integer id_time, n_times, itimes
      integer istart_t(2), iend_t(2)
      integer wrfyr,wrfmo,wrfdy,wrfhr,wrfmn,wrfsc
      integer ierr
!
!-----CAMx parameters and variable attributes for output
!
      integer mx3dmet,mx2dmet,mx3dcld,mx2dcld,mxkv,nlu,mxsrf
      parameter (mx3dmet =  6)
      parameter (mx2dmet = 14)
      parameter (mx3dcld =  9)
      parameter (mx2dcld =  7)
      parameter (mxkv    =  1)
      parameter (nlu     = 26)
      parameter (mxsrf   = 27)

      integer itzon,izone,iproj,istag
      integer n2dsrf,n3dmet,n2dmet,n3dcld,n2dcld,nkv
      integer o3dmet,o2dmet,o2dsrf
      integer kcmaq,kmyj,kysu
      integer x3dmet(mx3dmet),x2dmet(mx2dmet),x3dcld(mx3dcld),x2dcld(mx2dcld)

      real dxcamx,dycamx,x0camx,y0camx,olonin,olatin,tlat1in,tlat2in

      character*60 name_srf(mxsrf),units_srf(mxsrf),lnam_srf(mxsrf), &
                   desc_srf(mxsrf),coord_srf(mxsrf)
      character*60 name_3dmet(mx3dmet),units_3dmet(mx3dmet), &
                   lnam_3dmet(mx3dmet),desc_3dmet(mx3dmet),coord_3dmet(mx3dmet)
      character*60 name_2dmet(mx2dmet),units_2dmet(mx2dmet), &
                   lnam_2dmet(mx2dmet),desc_2dmet(mx2dmet),coord_2dmet(mx2dmet)
      character*60 name_3dcld(mx3dcld),units_3dcld(mx3dcld), &
                   lnam_3dcld(mx3dcld),desc_3dcld(mx3dcld),coord_3dcld(mx3dcld)
      character*60 name_2dcld(mx2dcld),units_2dcld(mx2dcld), &
                   lnam_2dcld(mx2dcld),desc_2dcld(mx2dcld),coord_2dcld(mx2dcld)
      character*60 name_kv(mxkv),units_kv(mxkv),lnam_kv(mxkv), &
                   desc_kv(mxkv),coord_kv(mxkv)
!
!-----NetCDF output arrays
!
      real, allocatable, dimension(:,:,:,:) :: out3d
      real, allocatable, dimension(:,:,:)   :: out2d
!
!-----Misc variables
!
      integer addday,subday
      integer sdate,edate,stime,etime,ioffset,joffset,k,nfiles, &
              wdate,wtime,cdate,ctime,cdatep,ctimep,i,j,kk
      integer startdat
      integer snodate,snotime,nxs,nys,isno,osno,snowagin
      integer dtout,kup,nfrc,n,cmn,cmnp,nout_step
      integer kpbl_cmaq,kpbl_ysu,kpbl_myj
      integer ktop
      integer, allocatable, dimension(:)  ::  kz1,kz2,kzin

      real xcamx,ycamx
      real lwat,tkemin,gamma,deltax,th,qvc,press0,wind,ustar,eli,wstar,kvmin, &
           risfc,sumtkz,sumtk,xinf,q,z,dzz,eee,al1,rwat,swat,gwat,delz,volrat, &
           sumf,sumff,osum
      real wsum,psum,esum,dsum,flx
      real qhi,qlo
      real dz1,dz2,zr,critk
      real, allocatable, dimension(:) ::  zm,zzc,thetav,uwind,vwind,  &
                                          tkep,elp,rkc,ttt,ppp,qqq,ccc

      character*200 fname,note,kvfil
      character*60 progname
      character*10 kvmeth,scmeth,project
      character*6 lucat
      logical lfirstw,lfirstc,lnew,lstagw,lnest,ldiag,lwindow,lstrat,debug,lice
!
      data progname /'WRFCAMx v5.2 01-10-22'/
!
      data lfirstw /.true./
      data lfirstc /.true./
      data lstagw /.true./
      data lnest  /.false./
      data debug /.false./
      data tkemin /1.e-6/
      data gamma /0.286/
!
!-----Initialize variables
!
      ioffset = -1
      joffset = -1
      izone = 0 
      olonin = 0.
      olatin = 0.
      tlat1in = 0.
      tlat2in = 0.
      kcmaq  = 0
      kmyj   = 0
      kysu   = 0
!
!-----Prepare output files for chunking/compression
!
      ncf_compress = .false.
#ifdef CHUNK
      ncf_compress = .true.
      call ncf_set_cache()
#endif
!
!-----Get user-supplied inputs
!
      read(*,'(20x,a)') note
      write(*,'(a)') note
      read(*,'(20x,l)') lnest
      write(*,'(a,l10)') '     This is a CAMx Nested Grid:',lnest
      if (lnest) then
        write(*,'(a,l10)') ' Adding Nested Grid Buffer Cells'
      else
        write(*,'(a,l10)') '     No Nested Grid Buffer Cells'
      endif
      read(*,'(20x,l)') ldiag
      write(*,'(a,l10)') '       Output Diagnostic Fields:',ldiag
      read(*,'(20x,l)') lice
      write(*,'(a,l10)') '     Adjust water LU for seaice:',lice
      read(*,'(20x,a10)') kvmeth
      if (kvmeth.ne.'MYJ' .and. kvmeth.ne.'CMAQ' .and. &
          kvmeth.ne.'YSU' .and. kvmeth.ne.'ALL') then
        write(*,*) 'Unknown Kv Method keyword: ',kvmeth
        stop
      endif
      write(*,'(a,1x,a10)') '                      Kv Method:',kvmeth
      read(*,'(20x,f10.0)') kvmin
      write(*,'(a,f10.2)')  '                     Minimum Kv:',kvmin
      read(*,'(20x,a10)') project
      write(*,'(a,1x,a10)') '                     Projection:',project
      read(*,'(20x,a10)') scmeth
      if (scmeth.ne.'KF' .and. scmeth.ne.'DIAG' .and. scmeth.ne.'NONE') then
        write(*,*) 'Unknown Sub-grid Cloud Method keyword: ',scmeth
        stop
      endif
      write(*,'(a,1x,a10)') '   Diagnose sub-grid convection:',scmeth
      read(*,'(20x,l)') lstrat
      write(*,'(a,l10)')  '   Diagnose sub-grid stratiform:',lstrat

      read(*,'(20x,a)') fname
      read(fname,*) sdate,edate
      stime = mod(sdate,100)
      etime = mod(edate,100)
      sdate = sdate/100
      edate = edate/100
      startdat = sdate
      call juldate(sdate)
      call juldate(edate)
      if (stime.gt.23) then
        stime = stime - 24
        sdate = addday(sdate)
      endif
      if (etime.gt.23) then
        etime = etime - 24
        edate = addday(edate)
      endif
      stime = 100*stime
      etime = 100*etime
      write(*,'(a,i6.5,i5.4)') '   Start date/time (YYJJJ HHMM):', &
                                sdate,stime
      write(*,'(a,i6.5,i5.4)') '     End date/time (YYJJJ HHMM):', &
                                edate,etime
      cdate = sdate
      ctime = stime

      read(*,'(20x,a)') fname
      read(fname,*) dtout
      write(*,'(a,i10)')       '          WRF output freq (min):',dtout
      if (dtout.lt.60 .and. amod(60.,float(dtout)).ne.0.) then
        write(*,*)'Met output frequency does not divide an hour evenly'
        write(*,*)'Program stopping'
        stop
      endif
      read(*,'(20x,i10)') itzon
      write(*,'(a,i10)')       '                      Time zone:',itzon
      read(*,'(20x,a)') fname
      read(fname,*) nxc,nyc,nzc
      write(*,'(a,3i10)')'                 CAMx grid size:',nxc,nyc,nzc

      if (lnest) then
        nxc = nxc + 2
        nyc = nyc + 2
      endif
      allocate(kz1(nzc))
      allocate(kz2(nzc))
      allocate(kzin(nzc))
      allocate(zm(nzc))
      allocate(zzc(nzc))
      allocate(thetav(nzc))
      allocate(uwind(nzc))
      allocate(vwind(nzc))
      allocate(tkep(nzc))
      allocate(elp(nzc))
      allocate(rkc(nzc))
      allocate(ttt(nzc))
      allocate(ppp(nzc))
      allocate(qqq(nzc))
      allocate(ccc(nzc))
      allocate(snowage(nxc,nyc))
      allocate(camxlon(nxc,nyc))
      allocate(camxlat(nxc,nyc))

      read(*,'(20x,a)') fname
      call parsedxy(fname,dxcamx,dycamx)
      write(*,'(a,2f10.3)')    '              CAMx grid spacing:',dxcamx,dycamx
      if (project.eq.'UTM       ') then
        iproj = 1
        read(*,'(20x,a)') fname
        read(fname,*) x0camx,y0camx,izone
        write(*,'(a,2f10.3)')  '    CAMx UTM Origin (SW corner):',x0camx,y0camx
        write(*,'(a,i10,/)')   '                       UTM zone:',izone
      elseif (project.eq.'LATLON    ') then
        iproj = 0
        read(*,'(20x,a)') fname
        read(fname,*) x0camx,y0camx
        write(*,'(a,2f10.3,/)')'CAMx Lat/Lon Origin (SW corner):',x0camx,y0camx
      elseif (project.eq.'LAMBERT   ') then
        iproj = 2
        read(*,'(20x,a)') fname
        read(fname,*) x0camx,y0camx,olonin,olatin,tlat1in,tlat2in
        write(*,'(a,2f10.3)')  '    CAMx LCP Origin (SW corner):',x0camx,y0camx
        write(*,'(a,4f10.3,/)') &
                '    CAMx LCP Projection Params :',olonin,olatin,tlat1in,tlat2in
      elseif (project.eq.'POLAR     ') then
        iproj = 4
        read(*,'(20x,a)') fname
        read(fname,*) x0camx,y0camx,olonin,olatin,tlat1in
        write(*,'(a,2f10.3)')  '    CAMx PS Origin (SW corner):',x0camx,y0camx
        write(*,'(a,3f10.3,/)') &
                         '    CAMx PS Projection Params :',olonin,olatin,tlat1in
      elseif (project.eq.'MERCATOR  ') then
        iproj = 5
        read(*,'(20x,a)') fname
        read(fname,*) x0camx,y0camx,olonin,olatin,tlat1in
        write(*,'(a,2f10.3)')  '    CAMx MRC Origin (SW corner):',x0camx,y0camx
        write(*,'(a,3f10.3,/)') &
                        '    CAMx MRC Projection Params :',olonin,olatin,tlat1in
      else
        write(*,'(a,a)') 'Unknown Projection Keyword: ',project
        stop
      endif
      if (project.ne.'LATLON' .and. dxcamx.ne.dycamx) then
        write(*,*)'CAMx DX must equal DY for cartesian projections:'
        write(*,*)'LAMBERT,UTM,POLAR,MERCATOR'
        write(*,*)'CAMx dx, dy: ',dxcamx,dycamx
        write(*,*)'Program stopping'
        stop
      endif
      if (lnest) then
        x0camx = x0camx - dxcamx
        y0camx = y0camx - dycamx
      endif

      read(*,'(20x,a)') fname
      read(fname,*) (kzin(k),k=1,nzc)

      read(*,'(20x,a)') fname
      call ncf_createfile(fname,'CAMx LU file',o2dsrf)
      write(*,*)'Opened CAMx LU file: ',fname

      read(*,'(20x,a)') fname
      call ncf_createfile(fname,'CAMx Met file',o3dmet)
      write(*,*)'Opened CAMx 3D Met file: ',fname

      read(*,'(20x,a)') fname
      call ncf_createfile(fname,'CAMx Met file',o2dmet)
      write(*,*)'Opened CAMx 2D Met file: ',fname

      read(*,'(20x,a)') fname
      if (kvmeth.eq.'CMAQ' .or. kvmeth.eq.'ALL') then
        kvfil = trim(fname)//'.CMAQ'
        call ncf_createfile(kvfil,'CAMx Kv file (CMAQ)',kcmaq)
        write(*,*)'Opened CAMx Kv file (CMAQ): ',kvfil
      endif
      if (kvmeth.eq.'MYJ' .or. kvmeth.eq.'ALL') then
        kvfil = trim(fname)//'.MYJ'
        call ncf_createfile(kvfil,'CAMx Kv file (MYJ)',kmyj)
        write(*,*)'Opened CAMx Kv file (MYJ): ',kvfil
      endif
      if (kvmeth.eq.'YSU' .or. kvmeth.eq.'ALL') then
        kvfil = trim(fname)//'.YSU'
        call ncf_createfile(kvfil,'CAMx Kv file (YSU)',kysu)
        write(*,*)'Opened CAMx Kv file (YSU): ',kvfil
      endif

      read(*,'(20x,a)') fname
      if (fname.ne.' ') then
        isno = 18
        open(unit=isno,file=fname,form='unformatted')
        write(*,*)'Opened input intermediate snow age file: ',fname
      else
        isno = 0
        write(*,*)'Skipping input intermediate snow age file'
        write(*,*)'Initial snow age will be set to user-defined value'
      endif
      read(*,'(20x,a)') fname
      if (isno.eq.0) then
        read(fname,*) snowagin
        write(*,*)'Initial snow age (hours):',snowagin
      endif
      read(*,'(20x,a)') fname
      if (fname.ne.' ') then
        osno = 19
        open(unit=osno,file=fname,form='unformatted')
        write(*,*)'Opened output intermediate snow age file: ',fname
      else
        osno = 0
        write(*,*)'Intermediate snow age file will not be output'
      endif
!
!-----Determine WRF layer mapping in CAMx layers
!
      do k = 1,nzc
        kz2(k) = kzin(k)
        if (k.eq.1) then
          kz1(k) = 1
        else
          kz1(k) = kz2(k-1) + 1
        endif
        write(*,'(a,i2,a,i2,a,i3)')'CAMx layer: ',k, &
              ' contains WRF layers ',kz1(k),' to',kz2(k)
      enddo
!
!-----Determine lat/lon of CAMx grid cell midpoints
!
      do j = 1,nyc
        ycamx = y0camx + dycamx*(j - 0.5)
        do i = 1,nxc
          xcamx = x0camx + dxcamx*(i - 0.5)
          if (project.eq.'LATLON    ') then
            camxlon(i,j) = xcamx
            camxlat(i,j) = ycamx
          elseif (project.eq.'LAMBERT   ') then
            call lcpgeo(1,olonin,olatin,tlat1in,tlat2in,xcamx,ycamx, &
                        camxlon(i,j),camxlat(i,j))
          elseif (project.eq.'UTM       ') then
            call utmgeo(1,izone,xcamx,ycamx,camxlon(i,j),camxlat(i,j))
          elseif (project.eq.'POLAR     ') then
            call pspgeo(1,olonin,olatin,tlat1in,xcamx,ycamx,camxlon(i,j), &
                        camxlat(i,j))
          elseif (project.eq.'MERCATOR  ') then
            call mrcgeo(1,olonin,olatin,tlat1in,xcamx,ycamx,camxlon(i,j), &
                        camxlat(i,j))
          endif
        enddo
      enddo
!
!-----Read intermediate snow age file
!
      snowage = 0.
      if (isno.gt.0) then
        write(*,*)'Reading input snow age file'
        read(isno) snodate,snotime,nxs,nys
        if (nxs.ne.nxc .or. nys.ne.nyc) then
          write(*,*)'Grid dimensions do not match'
          write(*,*)'CAMx grid: ',nxc,nyc
          write(*,*)'Snow grid: ',nxs,nys
          write(*,*)'STOPPING'
          stop
        endif
        if (snodate.ne.sdate .or. snotime.ne.stime) then
          write(*,*)'Snow age date/time does not match CAMx start date/time'
          write(*,*)'CAMx: ',sdate,stime
          write(*,*)'Snow: ',snodate,snotime
          write(*,*)'STOPPING'
          stop
        endif
        read(isno) ((snowage(i,j),i=1,nxs),j=1,nys)
      else
        do j = 1,nyc
          do i = 1,nxc
            snowage(i,j) = float(snowagin)
          enddo
        enddo
      endif 
!
!-----Loop over number of input WRF files
!
      nfiles = 0
      lnew = .false.
 200  read(*,'(20x,a)',end=201) fname
      rcode = nf_open(fname, NF_NOWRITE, cdfid)
      if (rcode .ne. 0) then
        write(*,*) 'Error opening WRF netcdf file ',fname
        write(*,*)'Program stopping'
        stop
      endif
      write(*,*)
      write(*,*) 'Opened input WRF file: ',fname
      write(*,*) 'cdfid = ', cdfid
      nfiles = nfiles + 1
      if (nfiles.gt.1) lnew = .true.
!
!-----Get the times on this WRF file
!
      rcode = nf_inq_varid(cdfid, 'Times', id_var)
      if (rcode .ne. 0) then
        write(*,*)'Cannot find any times on file'
        write(*,*)'This looks like a static file'
        write(*,*)'Program stopping'
        stop
      else
        id_time = ncvid(cdfid, 'Times', rcode)
      endif
      rcode = nf_inq_var(cdfid, id_time, varnam, ivtype, ndims, dimids, natts)
      do i = 1,ndims
        rcode = nf_inq_dimlen(cdfid, dimids(i), dims(i))
      enddo
      n_times = dims(2)
      if (allocated(times)) deallocate(times)
      allocate(times(n_times))

      write(*,'(/,a,t45,i6.5,i5.4)')  &
              'Looking for CAMx date/time (YYJJJ HHMM):',cdate,ctime
!
!-----Loop over times in this WRF output file to find the desired CAMx date/time
!
      do 100 itimes = 1,n_times
        istart_t(1) = 1
        iend_t(1)   = dims(1)
        istart_t(2) = itimes
        iend_t(2)   = 1
        rcode = NF_GET_VARA_TEXT(cdfid,id_time,istart_t,iend_t,times(itimes))
!
!-----Convert date from YYYYMMDDHH to a julian day and convert hour to
!     user-specified time zone
!
        read(times(itimes)(1:4),'(i4)') wrfyr
        read(times(itimes)(6:7),'(i2)') wrfmo
        read(times(itimes)(9:10),'(i2)') wrfdy
        read(times(itimes)(12:13),'(i2)') wrfhr
        read(times(itimes)(15:16),'(i2)') wrfmn
        read(times(itimes)(18:19),'(i2)') wrfsc
!
        wrfmn = wrfmn + nint(wrfsc/60.)
        nfrc = 60/dtout
        do n = 0,nfrc
          if (abs(wrfmn - n*dtout).lt.2) then
            wrfmn = n*dtout
            goto 20
          endif
        enddo
        write(*,'(a,a)')'WRF output frequency does not synchronize to', &
                        ' user-specified output inverval'
        write(*,*)'Program stopping'
        stop
!
 20     wdate = mod(wrfyr,100)*10000 + wrfmo*100 + wrfdy
        call juldate(wdate)
        if (wrfmn.gt.59) then
          wrfmn = wrfmn - 60
          wrfhr = wrfhr + 1
          if (wrfhr.gt.23) then
            wrfhr = wrfhr - 24
            wdate = addday(wdate)
          endif
        endif
        wrfhr = wrfhr - itzon
        if (wrfhr.lt.0) then
          wrfhr = wrfhr + 24
          wdate = subday(wdate)
        endif
        if (wrfhr.gt.23) then
          wrfhr = wrfhr - 24
          wdate = addday(wdate)
        endif
        wtime = 100*wrfhr + wrfmn
!
!-----Is this the date/time we are looking for?
!
        if (wdate.ne.cdate .or. wtime.ne.ctime) then
!
!-----If the second or higher WRF file was just opened, it has only one
!     timestamp, and WRF date/time doesn't match CAMx date/time, then assume
!     the user is using this WRF file to reset rain rates from a different WRF
!     run that overlaps the previous WRF file in time.
!
          if (lnew .and. n_times.eq.1) then
            call readwrf(lfirstw,lnew,cdfid,itimes,project,olatin,olonin, &
                         tlat1in,tlat2in,x0camx,y0camx,dxcamx,dycamx, &
                         deltax,kvmeth,kcmaq,kmyj,kysu, &
                         scmeth,lstrat,ioffset,joffset,izone,dtout,lucat, &
                         ldiag,lwindow,lice,debug)
          endif
          goto 100
        endif
!
!-----WRF date/time = CAMx date/time: set special date/time variables for
!     precip and clouds
!
        write(*,'(a,t45,i6.5,i5.4,/)') &
              '       Found WRF date/time (YYJJJ HHMM):',wdate,wtime
        cdatep = cdate
        ctimep = ctime/100
        cmnp = mod(ctime,100)
        if (dtout.lt.60) then
          cmnp = cmnp - dtout
          if (cmnp.lt.0) then
            cmnp = cmnp + 60
            ctimep = ctimep - 1
          endif
        else
          ctimep = ctimep - dtout/60
        endif
        if (ctimep.lt.0) then
          ctimep = ctimep + 24
          cdatep = subday(cdatep)
        endif
        ctimep = 100*ctimep + cmnp
!
!-----If the second or higher WRF file was just opened and the timestamp
!     to read is not the first on the file, this is probably from a
!     different WRF run that overlaps the previous WRF file in time.
!     Read WRF data for the previous timestamp that we actually want so
!     that the rain rates are set correctly.
!
        if (lnew .and. itimes.gt.1) then
          call readwrf(lfirstw,lnew,cdfid,itimes-1,project,olatin,olonin, &
                       tlat1in,tlat2in,x0camx,y0camx,dxcamx,dycamx, &
                       deltax,kvmeth,kcmaq,kmyj,kysu, &
                       scmeth,lstrat,ioffset,joffset,izone,dtout,lucat, &
                       ldiag,lwindow,lice,debug)
        endif
        lnew = .false.
!
!-----Read the WRF file for this timestamp
!
        call readwrf(lfirstw,lnew,cdfid,itimes,project,olatin,olonin, &
                     tlat1in,tlat2in,x0camx,y0camx,dxcamx,dycamx, &
                     deltax,kvmeth,kcmaq,kmyj,kysu, &
                     scmeth,lstrat,ioffset,joffset,izone,dtout,lucat, &
                     ldiag,lwindow,lice,debug)
!
!-----Horizontally interpolate, vertically aggregate, WRF data to CAMx grid
!
!     This is where the user should substitute new routines that handle
!     interpolation to a specific grid type/projection
!
        if (lwindow) then
          lstagw = .true.
          call wrf_window(kz1,kz2,ioffset,joffset,nlu,kcmaq, &
                          kmyj,kysu,ldiag)
        else
          lstagw = .false.
          call wrf_interp(project,kz1,kz2,nlu,deltax,kcmaq, &
                          kmyj,kysu,ldiag)
        endif
!
        do j = 1,nyc
          do i = 1,nxc
!
!-----Check for (and fill) humidity holes (Qa < 1 ppm)
!
            do k = 1,nzc
              if (qac(i,j,k).le.1.) then
                qhi = 0.
                qlo = 0.
                if (k.gt.1) qlo = qac(i,j,k-1)
                do kup = k+1,nzc
                  if (qac(i,j,kup).gt.1.) then
                    qhi = qac(i,j,kup)
                    goto 30
                  endif
                enddo
 30             if (qhi.le.1. .and. qlo.le.1.) then
                  write(*,'(a,i3,a,i3,a)') 'Entire column (',i,',',j &
                                    ,') has unreasonably low humidty!'
                  write(*,'(a)') 'Stopping! Check WRF Data'
                  stop
                elseif (qhi.gt.1. .and. qlo.le.1.) then
                  do kk = 1,kup-1
                    qac(i,j,kk) = qhi
                  enddo
                elseif (qhi.le.1. .and. qlo.gt.1) then
                  do kk = k,nzc
                    qac(i,j,kk) = qlo
                  enddo
                else
                  do kk = k,kup-1
                    qac(i,j,kk) = (qhi + qlo)*0.5
                  enddo
                endif
                goto 31
              endif
            enddo
 31         continue
!
!-----Diagnose Kv profiles: MYJ, YSU, CMAQ options
!
            do k = 1,nzc
              zzc(k) = zhc(i,j,k)
              zm(k) = zzc(k)/2.
              if (k.gt.1) zm(k) = (zzc(k) + zzc(k-1))/2.

              th = tac(i,j,k)*(1000./pac(i,j,k))**gamma
              qvc = qac(i,j,k)*18./28.8/1.e6
              thetav(k) = th*(1. + 0.61*qvc)
              qqq(k) = qac(i,j,k)

              uwind(k) = uac(i,j,k)
              if (i.gt.1) uwind(k) = (uac(i,j,k) + uac(i-1,j,k))/2.
              vwind(k) = vac(i,j,k)
              if (j.gt.1) vwind(k) = (vac(i,j,k) + vac(i,j-1,k))/2.

              if (kmyj .ne. 0) then
                tkep(k) = tkc(i,j,k)
                elp(k) = elc(i,j,k)
              endif
            enddo
            press0 = pac(i,j,1) - zm(1)*(pac(i,j,2) - pac(i,j,1))/ &
                                        (zm(2) - zm(1))
            wind = sqrt(uwind(1)**2 + vwind(1)**2)

            pblc(i,j) = max(pblc(i,j),zhc(i,j,1))
            call micromet(tac(i,j,1),tsfc(i,j),pac(i,j,1),press0, &
                          zm(1),wind,z0c(i,j),pblc(i,j),ustar,eli,wstar, &
                          risfc)

            if (kysu.ne.0) then
              call kv_ysu(nzc,pblc(i,j),ustar,eli,wstar,zm,zzc, &
                          thetav,uwind,vwind,qqq,kvmin,risfc,rkc)
              do k = 1,nzc
                rkv_ysu(i,j,k) = rkc(k)
              enddo
            endif

            if (kcmaq.ne.0) then
              call kv_cmaq(nzc,pblc(i,j),ustar,eli,wstar,zm,zzc, &
                           thetav,uwind,vwind,kvmin,rkc)
              do k = 1,nzc
                rkv_cmaq(i,j,k) = rkc(k)
              enddo
            endif

            if (kmyj.ne.0) then
              call kv_tkemyj(nzc,thetav,uwind,vwind,zm,zzc, &
                            tkep,elp,kvmin,rkc)
              do k = 1,nzc
                rkv_myj(i,j,k) = rkc(k)
              enddo
            endif
!
!-----Prepare cloud fields
!
            kff(i,j) = 0.
            kft(i,j) = 0.
            do k = 1,nzc
              cwc(i,j,k) = 0.
              pwr(i,j,k) = 0.
              pws(i,j,k) = 0.
              pwg(i,j,k) = 0.
              cod(i,j,k) = 0.
              kfw(i,j,k) = 0.
              kfp(i,j,k) = 0.
              kfe(i,j,k) = 0.
              kfd(i,j,k) = 0.
            enddo
!
!-----Calculate mean CWC and PWC for each CAMx layer
!
            do k = 1,nzc
              lwat = 0.
              rwat = 0.
              swat = 0.
              gwat = 0.
              osum = 0.
              wsum = 0.
              psum = 0.
              esum = 0.
              dsum = 0.
              delz = zhc(i,j,k)
              if (k.gt.1) delz = zhc(i,j,k) - zhc(i,j,k-1)
              do kk = kz1(k),kz2(k)
                dzz = ztmp(i,j,kk)
                if (kk.gt.1) dzz = ztmp(i,j,kk) - ztmp(i,j,kk-1)
                lwat = lwat + cwtmp(i,j,kk)*dzz
                rwat = rwat + prtmp(i,j,kk)*dzz
                swat = swat + pstmp(i,j,kk)*dzz
                gwat = gwat + pgtmp(i,j,kk)*dzz
                osum = osum + odtmp(i,j,kk)
                wsum = wsum + kwtmp(i,j,kk)*dzz
                psum = psum + kptmp(i,j,kk)*dzz
                esum = esum + ketmp(i,j,kk)
                dsum = dsum + kdtmp(i,j,kk)
              enddo
              if (lwat.gt.0.) cwc(i,j,k) = lwat/delz
              if (rwat.gt.0.) pwr(i,j,k) = rwat/delz
              if (swat.gt.0.) pws(i,j,k) = swat/delz
              if (gwat.gt.0.) pwg(i,j,k) = gwat/delz
              if (osum.gt.0.) cod(i,j,k) = osum
              if (wsum.gt.0.) kfw(i,j,k) = wsum/delz
              if (psum.gt.0.) kfp(i,j,k) = psum/delz
              if (esum.gt.0.) kfe(i,j,k) = esum
              if (dsum.gt.0.) kfd(i,j,k) = dsum
            enddo
!
!-----Run through column and check K-F consistency
!
            if (scmeth.eq.'KF') then
              do k = 1,nzc
                if (kfw(i,j,k).gt.0.01) goto 40
              enddo
              kff(i,j) = 0.
              kft(i,j) = 0.
              do k = 1,nzc
                kfw(i,j,k) = 0.
                kfp(i,j,k) = 0.
                kfe(i,j,k) = 0.
                kfd(i,j,k) = 0.
              enddo
              goto 42
 40           continue
              kff(i,j) = kftmp(i,j)
              kft(i,j) = tctmp(i,j)
              flx = kfe(i,j,1) - kfd(i,j,1)
              ktop = 0
              do k = 2,nzc
                if (kfw(i,j,k).le.0.01) then
                  if (ktop.ne.0) goto 41
                else
                  ktop = k
                endif
                flx = flx + kfe(i,j,k) - kfd(i,j,k)
              enddo
 41           continue
              if (flx.gt.0.) then
                kfd(i,j,ktop) = kfd(i,j,ktop) + flx
              else
                kfe(i,j,ktop) = kfe(i,j,ktop) - flx
              endif
              do k = ktop+1,nzc
                kfw(i,j,k) = 0.
                kfp(i,j,k) = 0.
                kfe(i,j,k) = 0.
                kfd(i,j,k) = 0.
              enddo
            endif
 42         continue
!
!-----Update snow age for cells covered in snow
!
            if (snow(i,j).eq.0.) then
              snowage(i,j) = 0. 
            else
              if (.not.lfirstc) then
                snowage(i,j) = snowage(i,j) + dtout/60.
                if ((snow(i,j) - snowo(i,j)).ge.0.001) snowage(i,j) = 0.
              endif
            endif
            snowo(i,j) = snow(i,j)
!
          enddo
        enddo
!
!-----Calculate some remaining diagnostic fields (precip rate, total cloud OD,
!     PBL depth)
!
        do k = 1,nzc
          do j = 1,nyc
            do i = 1,nxc
              pwtr(i,j,k) = pwr(i,j,k) + pws(i,j,k) + pwg(i,j,k)
            enddo
          enddo
        enddo
        do j = 1,nyc
          do i = 1,nxc
            volrat = pwtr(i,j,1)/1.e6
            if (scmeth.eq.'KF') volrat = volrat + kfp(i,j,1)*kff(i,j)/1.e6
            rain(i,j) = (volrat/1.0e-7)**1.27
            codsum(i,j) = 0.
            do k = 1,nzc
              codsum(i,j) = codsum(i,j) + cod(i,j,k)
            enddo

            kpbl_cmaq = 1
            kpbl_myj  = 1
            kpbl_ysu  = 1
            do k = 2,nzc
              dz1 = zhc(i,j,k-1)
              if (k.gt.2) dz1 = zhc(i,j,k-1) - zhc(i,j,k-2)
              dz2 = zhc(i,j,k) - zhc(i,j,k-1)
              zr = dz2/dz1
              critk = 0.03*dz1*dz1*(1. + zr)/200.
              if (kcmaq.ne.0 .and. &
                  kpbl_cmaq.eq.k-1 .and. rkv_cmaq(i,j,k-1).gt.critk) &
                  kpbl_cmaq = k
              if (kmyj.ne.0 .and. &
                  kpbl_myj .eq.k-1 .and. rkv_myj(i,j,k-1).gt.critk ) &
                  kpbl_myj  = k
              if (kysu.ne.0 .and. &
                  kpbl_ysu .eq.k-1 .and. rkv_ysu(i,j,k-1).gt.critk ) &
                  kpbl_ysu  = k
            enddo
            if (kcmaq.ne.0) pbl_cmaq(i,j) = zhc(i,j,kpbl_cmaq)
            if (kmyj.ne.0)  pbl_myj(i,j)  = zhc(i,j,kpbl_myj)
            if (kysu.ne.0)  pbl_ysu(i,j)  = zhc(i,j,kpbl_ysu)
          enddo
        enddo
!
!-----Write to CAMx files
!
        if (lfirstc) then
          nout_step = 1
!
!-----Set CAMx file attributes and write header info
!
          call ncf_set_tstep(sdate,stime,edate,etime,dtout)
          call ncf_set_global(sdate,stime,startdat,dtout,progname,lnest, &
                              lstagw,itzon,izone,iproj,nxc,nyc,nzc,dxcamx, &
                              dycamx,x0camx,y0camx,olonin,olatin,tlat1in, &
                              tlat2in,note)
          call ncf_set_varatt_grid()
          call ncf_set_varatt_2dsrf(mxsrf,n2dsrf,name_srf,units_srf,lnam_srf, &
                                    desc_srf,coord_srf)
          call ncf_set_varatt_3dmet(mx3dmet,mx2dmet,n3dmet,n2dmet,ldiag,&
                                    kcmaq,kmyj,kysu,x3dmet,x2dmet, &
                                    name_3dmet,units_3dmet,lnam_3dmet, &
                                    desc_3dmet,coord_3dmet, &
                                    name_2dmet,units_2dmet,lnam_2dmet, &
                                    desc_2dmet,coord_2dmet)
          call ncf_set_varatt_3dcld(mx3dcld,mx2dcld,n3dcld,n2dcld,ldiag,&
                                    scmeth,x3dcld,x2dcld, &
                                    name_3dcld,units_3dcld,lnam_3dcld, &
                                    desc_3dcld,coord_3dcld, &
                                    name_2dcld,units_2dcld,lnam_2dcld, &
                                    desc_2dcld,coord_2dcld)
          call ncf_set_varatt_kv(mxkv,nkv,name_kv,units_kv,lnam_kv,desc_kv, &
                                 coord_kv)

!         ---2-D surface file---

          call ncf_wrt_dim('CAMx LU file',o2dsrf,nxc,nyc,1,n2dsrf)
          call ncf_wrt_global('CAMx LU file',o2dsrf,.true.,.false., &
                              'LANDUSE',n2dsrf,name_srf)
          call ncf_wrt_vars_grid('CAMx LU file',o2dsrf)
          call ncf_wrt_vars('CAMx LU file',o2dsrf,nxc,nyc,1,n2dsrf,name_srf, &
                            units_srf,lnam_srf,desc_srf,coord_srf,4)
          call ncf_enddef_file('CAMx LU file',o2dsrf)
          call ncf_wrt_data_grid('CAMx LU file',o2dsrf,nxc,nyc,nzc,x0camx, &
                                 y0camx,dxcamx,dycamx,camxlat,camxlon)

!         ---3-D met file---

          call ncf_wrt_dim('CAMx 3D Met file',o3dmet,nxc,nyc,nzc,n3dmet+n3dcld)
          call ncf_wrt_global2('CAMx 3D Met file',o3dmet,.false.,.false., &
                               '3D METEOROLOGY',n3dmet,n3dcld,name_3dmet, &
                               name_3dcld)
          call ncf_wrt_vars_grid('CAMx 3D Met file',o3dmet)
          call ncf_wrt_vars('CAMx 3D Met file',o3dmet,nxc,nyc,nzc,n3dmet, &
                            name_3dmet,units_3dmet,lnam_3dmet,desc_3dmet, &
                            coord_3dmet,4)
          call ncf_wrt_vars('CAMx 3D Met file',o3dmet,nxc,nyc,nzc,n3dcld, &
                            name_3dcld,units_3dcld,lnam_3dcld,desc_3dcld, &
                            coord_3dcld,4)
          call ncf_enddef_file('CAMx 3D Met file',o3dmet)
          call ncf_wrt_data_grid('CAMx 3D Met file',o3dmet,nxc,nyc,nzc,x0camx, &
                                 y0camx,dxcamx,dycamx,camxlat,camxlon)

!         ---2-D met file---

          call ncf_wrt_dim('CAMx 2D Met file',o2dmet,nxc,nyc,1,n2dmet+n2dcld)
          call ncf_wrt_global2('CAMx 2D Met file',o2dmet,.false.,.true., &
                               '2D METEOROLOGY',n2dmet,n2dcld,name_2dmet, &
                               name_2dcld)
          call ncf_wrt_vars_grid('CAMx 2D Met file',o2dmet)
          call ncf_wrt_vars('CAMx 2D Met file',o2dmet,nxc,nyc,1,n2dmet, &
                            name_2dmet,units_2dmet,lnam_2dmet,desc_2dmet, &
                            coord_2dmet,4)
          call ncf_wrt_vars('CAMx 2D Met file',o2dmet,nxc,nyc,1,n2dcld, &
                            name_2dcld,units_2dcld,lnam_2dcld,desc_2dcld, &
                            coord_2dcld,4)
          call ncf_enddef_file('CAMx 2D Met file',o2dmet)
          call ncf_wrt_data_grid('CAMx 2D Met file',o2dmet,nxc,nyc,1,x0camx, &
                                 y0camx,dxcamx,dycamx,camxlat,camxlon)

!         ---Various Kv files depending on which are selected---

          if (kcmaq.ne.0) then
            call ncf_wrt_dim('CAMx Kv file (CMAQ)',kcmaq,nxc,nyc,nzc,nkv)
            call ncf_wrt_global('CAMx Kv file (CMAQ)',kcmaq,.false.,.false., &
                                'VERTICAL DIFFUSIVITY (CMAQ)',nkv,name_kv)
            call ncf_wrt_vars_grid('CAMx Kv file (CMAQ)',kcmaq)
            call ncf_wrt_vars('CAMx Kv file (CMAQ)',kcmaq,nxc,nyc,nzc,nkv, &
                              name_kv,units_kv,lnam_kv,desc_kv,coord_kv,4)
            call ncf_enddef_file('CAMx Kv file (CMAQ)',kcmaq)
            call ncf_wrt_data_grid('CAMx Kv file (CMAQ)',kcmaq,nxc,nyc,nzc, &
                                   x0camx,y0camx,dxcamx,dycamx,camxlat,camxlon)
          endif
          if (kmyj.ne.0) then
            call ncf_wrt_dim('CAMx Kv file (MYJ)',kmyj,nxc,nyc,nzc,nkv)
            call ncf_wrt_global('CAMx Kv file (MYJ)',kmyj,.false.,.false., &
                                'VERTICAL DIFFUSIVITY (MYJ)',nkv,name_kv)
            call ncf_wrt_vars_grid('CAMx Kv file (MYJ)',kmyj)
            call ncf_wrt_vars('CAMx Kv file (MYJ)',kmyj,nxc,nyc,nzc,nkv, &
                              name_kv,units_kv,lnam_kv,desc_kv,coord_kv,4)
            call ncf_enddef_file('CAMx Kv file (MYJ)',kmyj)
            call ncf_wrt_data_grid('CAMx Kv file (MYJ)',kmyj,nxc,nyc,nzc, &
                                   x0camx,y0camx,dxcamx,dycamx,camxlat,camxlon)
          endif
          if (kysu.ne.0) then
            call ncf_wrt_dim('CAMx Kv file (YSU)',kysu,nxc,nyc,nzc,nkv)
            call ncf_wrt_global('CAMx Kv file (YSU)',kysu,.false.,.false., &
                                'VERTICAL DIFFUSIVITY (YSU)',nkv,name_kv)
            call ncf_wrt_vars_grid('CAMx Kv file (YSU)',kysu)
            call ncf_wrt_vars('CAMx Kv file (YSU)',kysu,nxc,nyc,nzc,nkv, &
                              name_kv,units_kv,lnam_kv,desc_kv,coord_kv,4)
            call ncf_enddef_file('CAMx Kv file (YSU)',kysu)
            call ncf_wrt_data_grid('CAMx Kv file (YSU)',kysu,nxc,nyc,nzc, &
                                   x0camx,y0camx,dxcamx,dycamx,camxlat,camxlon)
          endif
!
!-----Landuse and topography: renormalize landuse fractions
!
          do j = 1,nyc
            do i = 1,nxc
              if (lucat.eq.'USGS  ' .or. lucat.eq.'NLCD  ' .or. &
                  lucat.eq.'IGBP  ') then
                sumf = 0.
                sumff = 0.
                do k = 1,nlu
                  if (lucx(i,j,k).lt.0.) lucx(i,j,k) = 0.
                  sumf = sumf + lucx(i,j,k)
                enddo
                do k = 1,nlu
                  lucx(i,j,k) = lucx(i,j,k)/sumf
                  sumff = sumff + lucx(i,j,k)
                enddo
                if (sumff.lt.0.99 .or. sumff.gt.1.01) then
                  write(*,*)'Land cover total not = 1!',i,j,sumff
                  stop
                endif
              else
                do k = 1,nlu
                  lucx(i,j,k) = 0.
                enddo
              endif
            enddo
          enddo

          call ncf_wrt_data_tstep('CAMx LU file',o2dsrf,nout_step,n2dsrf)
          call ncf_wrt_data('CAMx LU file',o2dsrf,1,nxc,nyc,1,nout_step, &
                            name_srf(1),topcx)
          call ncf_wrt_data('CAMx LU file',o2dsrf,n2dsrf-1,nxc,nyc,1, &
                            nout_step,name_srf(2),lucx)
          ierr = nf_close(o2dsrf)

          call stats('Landuse fraction    ',lucx,nxc,nyc,nlu)
          write(*,*)
          call stats('Topography (m)      ',topcx,nxc,nyc,1)
          write(*,*)

          lfirstc = .false.
        endif
!
!-----Time-varying met fields
!
!-----3D Met File (MAKE SURE MET_3DNAME INDEX MATCHES VARIABLES TO SEND)
!
        call ncf_wrt_data_tstep('CAMx 3D Met file',o3dmet,nout_step, &
                                n3dmet+n3dcld)
        allocate( out3d(nxc,nyc,nzc,n3dmet) )
        out3d(:,:,:,x3dmet(i3dz)) = zhc(:,:,:)
        out3d(:,:,:,x3dmet(i3dp)) = pac(:,:,:)
        out3d(:,:,:,x3dmet(i3dt)) = tac(:,:,:)
        out3d(:,:,:,x3dmet(i3dq)) = qac(:,:,:)
        out3d(:,:,:,x3dmet(i3du)) = uac(:,:,:)
        out3d(:,:,:,x3dmet(i3dv)) = vac(:,:,:)
        call ncf_wrt_data('CAMx 3D Met file',o3dmet,n3dmet,nxc,nyc,nzc, &
                          nout_step,name_3dmet,out3d)
        deallocate( out3d )
!
!-----2D Met File (MAKE SURE MET_2DNAME INDEX MATCHES VARIABLES TO SEND)
!
        call ncf_wrt_data_tstep('CAMx 2D Met file',o2dmet,nout_step, &
                                n2dmet+n2dcld)
        allocate( out2d(nxc,nyc,n2dmet) )
        out2d(:,:,x2dmet(i2dts)) = tsfc(:,:)
        out2d(:,:,x2dmet(i2dsd)) = snow(:,:)
        out2d(:,:,x2dmet(i2dsa)) = snowage(:,:)
        if (ldiag) then
          out2d(:,:,x2dmet(i2du10)) = u10(:,:)
          out2d(:,:,x2dmet(i2dv10)) = v10(:,:)
          out2d(:,:,x2dmet(i2dt2))  = t2(:,:)
          out2d(:,:,x2dmet(i2dsw))  = swsfc(:,:)
          out2d(:,:,x2dmet(i2dsm))  = soilm(:,:)
          out2d(:,:,x2dmet(i2dst))  = soilt(:,:)
          out2d(:,:,x2dmet(i2dty))  = styp(:,:)
          out2d(:,:,x2dmet(i2dwrf)) = pblc(:,:)
          if (kcmaq.ne.0) &
            out2d(:,:,x2dmet(i2dcmaq)) = pbl_cmaq(:,:)
          if (kmyj.ne.0) &
            out2d(:,:,x2dmet(i2dmyj)) = pbl_myj(:,:)
          if (kysu.ne.0) &
            out2d(:,:,x2dmet(i2dysu)) = pbl_ysu(:,:)
        endif
        call ncf_wrt_data('CAMx 2D Met file',o2dmet,n2dmet,nxc,nyc,1, &
                          nout_step,name_2dmet,out2d)
        deallocate( out2d )
 
        call stats('Height (m)          ',zhc,nxc,nyc,nzc)
        write(*,*)
        call stats('Pressure (mb)       ',pac,nxc,nyc,nzc)
        write(*,*)
        if (ldiag) call stats('U at 10 m (m/s)     ',u10,nxc,nyc,1)
        call stats('U-comp (m/s)        ',uac,nxc,nyc,nzc)
        write(*,*)
        if (ldiag) call stats('V at 10 m (m/s)     ',v10,nxc,nyc,1)
        call stats('V-comp (m/s)        ',vac,nxc,nyc,nzc)
        write(*,*)
        call stats('Sfc Temp (K)        ',tsfc,nxc,nyc,1)
        if (ldiag) call stats('Temp at 2 m (K)     ',t2,nxc,nyc,1)
        call stats('Temperature (K)     ',tac,nxc,nyc,nzc)
        write(*,*)
        call stats('Humidity (ppm)      ',qac,nxc,nyc,nzc)
        write(*,*)
        call stats('Snow Depth (cm)     ',1000.*snow,nxc,nyc,1)
        call stats('Snow Age (hr)       ',snowage,nxc,nyc,1)
        write(*,*)
        if (ldiag) then
          call stats('Sfc shortwave (W/m2)',swsfc,nxc,nyc,1)
          call stats('Soil moist (m3/m3)  ',soilm,nxc,nyc,1)
          call stats('Soil temperature (K)',soilt,nxc,nyc,1)
          call stats('Soil type           ',styp,nxc,nyc,1)
          write(*,*)
        endif
        call stats('WRF PBL Depth (m)   ',pblc,nxc,nyc,1)
        write(*,*)
!
!-----Kv file
!
        if (kcmaq.ne.0) then
          call ncf_wrt_data_tstep('CAMx Kv file (CMAQ)',kcmaq,nout_step,nkv)
          call ncf_wrt_data('CAMx Kv file (CMAQ)',kcmaq,1,nxc,nyc,nzc, &
                            nout_step,name_kv,rkv_cmaq)
          call stats('Kv CMAQ (m2/s)      ',rkv_cmaq,nxc,nyc,nzc)
          call stats('CMAQ PBL Depth (m)  ',pbl_cmaq,nxc,nyc,1)
          write(*,*)
        endif
        if (kmyj.ne.0) then
          call ncf_wrt_data_tstep('CAMx Kv file (MYJ)',kmyj,nout_step,nkv)
          call ncf_wrt_data('CAMx Kv file (MYJ)',kmyj,1,nxc,nyc,nzc, &
                            nout_step,name_kv,rkv_myj)
          call stats('Kv MYJ (m2/s)       ',rkv_myj,nxc,nyc,nzc)
          call stats('MYJ PBL Depth (m)   ',pbl_myj,nxc,nyc,1)
          write(*,*)
        endif
        if (kysu.ne.0) then
          call ncf_wrt_data_tstep('CAMx Kv file (YSU)',kysu,nout_step,nkv)
          call ncf_wrt_data('CAMx Kv file (YSU)',kysu,1,nxc,nyc,nzc, &
                            nout_step,name_kv,rkv_ysu)
          call stats('Kv YSU (m2/s)       ',rkv_ysu,nxc,nyc,nzc)
          call stats('YSU PBL Depth (m)   ',pbl_ysu,nxc,nyc,1)
          write(*,*)
        endif
!
!-----Cloud/rain fields
!
        if (cdatep.lt.sdate .or. &
           (cdatep.eq.sdate .and. ctimep.lt.stime)) goto 101
        write(*,'(a,t30,i6.5,i5.4,/)') &
              'Cld/rn date/time (YYJJJ HHMM):',cdatep,ctimep

        allocate( out3d(nxc,nyc,nzc,n3dcld) )
        out3d(:,:,:,x3dcld(i3dcwtr)) = cwc(:,:,:)
        out3d(:,:,:,x3dcld(i3drwtr)) = pwr(:,:,:)
        out3d(:,:,:,x3dcld(i3dswtr)) = pws(:,:,:)
        out3d(:,:,:,x3dcld(i3dgwtr)) = pwg(:,:,:)
        out3d(:,:,:,x3dcld(i3dcod))  = cod(:,:,:)
        if (scmeth.eq.'KF') then
          out3d(:,:,:,x3dcld(i3dkf_cwtr)) = kfw(:,:,:)
          out3d(:,:,:,x3dcld(i3dkf_pwtr)) = kfp(:,:,:)
          out3d(:,:,:,x3dcld(i3dkf_ent))  = kfe(:,:,:)
          out3d(:,:,:,x3dcld(i3dkf_det))  = kfd(:,:,:)
        endif
        call ncf_wrt_data('CAMx 3D Met file',o3dmet,n3dcld,nxc,nyc,nzc, &
                          nout_step-1,name_3dcld,out3d)
        deallocate( out3d )

        allocate( out2d(nxc,nyc,n2dcld) )
        if (ldiag) then
          out2d(:,:,x2dcld(i2dprat)) = rain(:,:)
          out2d(:,:,x2dcld(i2dtcod)) = codsum(:,:)
          out2d(:,:,x2dcld(i2dct))   = ctop(:,:)
          out2d(:,:,x2dcld(i2dcp))   = cape(:,:)
          out2d(:,:,x2dcld(i2dcf))   = cldfrc(:,:)
        endif
        if (scmeth.eq.'KF') then
          out2d(:,:,x2dcld(i2dkf_frc))  = kff(:,:)
          out2d(:,:,x2dcld(i2dkf_tscl)) = kft(:,:)
        endif
        call ncf_wrt_data('CAMx 2D Met file',o2dmet,n2dcld,nxc,nyc,1, &
                          nout_step-1,name_2dcld,out2d)
        deallocate( out2d )

        call stats('Cloud Water (g/m3)  ',cwc,nxc,nyc,nzc)
        write(*,*)
        call stats('Precip Water(g/m3)  ',pwtr,nxc,nyc,nzc)
        write(*,*)
        call stats('Cloud Optical Depth ',cod,nxc,nyc,nzc)
        if (ldiag) call stats('Cloud Total OD      ',codsum,nxc,nyc,1)
        write(*,*)
        if (ldiag) then
          call stats('Precip Rate (mm/hr) ',rain,nxc,nyc,1)
          call stats('Conv cloud top (km) ',ctop,nxc,nyc,1)
          call stats('CAPE (J/kg)         ',cape,nxc,nyc,1)
          call stats('Cloud Frac (dec)    ',cldfrc,nxc,nyc,1)
          write(*,*)
        endif

        if (scmeth.eq.'KF') then
          call stats('KF Cloud Fraction   ',kff,nxc,nyc,1)
          call stats('KF Cloud Tscale  (s)',kft,nxc,nyc,1)
          write(*,*)
          call stats('KF Cloud Water(g/m3)',kfw,nxc,nyc,nzc)
          write(*,*)
          call stats('KF Precp Water(g/m3)',kfp,nxc,nyc,nzc)
          write(*,*)
          call stats('KF Entrain  (g/m2/s)',kfe,nxc,nyc,nzc)
          write(*,*)
          call stats('KF Detrain  (g/m2/s)',kfd,nxc,nyc,nzc)
          write(*,*)
        endif
!
!-----Check if this is the end of the time period to process
!
 101    continue
        snotime = ctime
        snodate = cdate

        cmn = mod(ctime,100)
        ctime = ctime/100
        if (dtout.lt.60) then
          cmn = cmn + dtout
          if (cmn.gt.59) then
            cmn = cmn - 60
            ctime = ctime + 1
          endif
        else
          ctime = ctime + dtout/60
        endif
        if (ctime.gt.23) then
          ctime = ctime - 24
          cdate = addday(cdate)
        endif
        ctime = 100*ctime + cmn
!
!-----This is the end of time; write intermediate snow file and final
!     cloud/rain fields (all zero as this is a dummy time for clouds)
!
        if ((cdate.eq.edate .and. ctime.gt.etime) .or. cdate.gt.edate) then
          if (osno.gt.0) then
            write(osno) snodate,snotime,nxc,nyc
            write(osno) ((snowage(i,j),i=1,nxc),j=1,nyc)
          endif

          allocate( out3d(nxc,nyc,nzc,n3dcld) )
          out3d(:,:,:,x3dcld(i3dcwtr)) = 0.
          out3d(:,:,:,x3dcld(i3drwtr)) = 0.
          out3d(:,:,:,x3dcld(i3dswtr)) = 0.
          out3d(:,:,:,x3dcld(i3dgwtr)) = 0.
          out3d(:,:,:,x3dcld(i3dcod))  = 0.
          if (scmeth.eq.'KF') then
            out3d(:,:,:,x3dcld(i3dkf_cwtr)) = 0.
            out3d(:,:,:,x3dcld(i3dkf_pwtr)) = 0.
            out3d(:,:,:,x3dcld(i3dkf_ent))  = 0.
            out3d(:,:,:,x3dcld(i3dkf_det))  = 0.
          endif
          call ncf_wrt_data('CAMx 3D Met file',o3dmet,n3dcld,nxc,nyc,nzc, &
                            nout_step,name_3dcld,out3d)
          deallocate( out3d )

          allocate( out2d(nxc,nyc,n2dcld) )
          if (ldiag) then
            out2d(:,:,x2dcld(i2dprat)) = 0.
            out2d(:,:,x2dcld(i2dtcod)) = 0.
            out2d(:,:,x2dcld(i2dct))   = 0.
            out2d(:,:,x2dcld(i2dcp))   = 0.
            out2d(:,:,x2dcld(i2dcf))   = 0.
          endif
          if (scmeth.eq.'KF') then
            out2d(:,:,x2dcld(i2dkf_frc))  = 0.
            out2d(:,:,x2dcld(i2dkf_tscl)) = 0.
          endif
          call ncf_wrt_data('CAMx 2D Met file',o2dmet,n2dcld,nxc,nyc,1, &
                            nout_step,name_2dcld,out2d)
          deallocate( out2d )

          write(*,'(/,a,/)')'Reached CAMx end date/hour; RUN FINISHED'
          goto 999
        endif

        nout_step = nout_step + 1
        write(*,'(/,a,t45,i6.5,i5.4)')  &
              'Looking for CAMx date/time (YYJJJ HHMM):',cdate,ctime

 100  continue
      rcode = nf_close(cdfid)
      goto 200
 201  write(*,'(/,a,a,/)')'End of WRF data, CAMx date/time not found;', &
                          ' RUN INCOMPLETE'
 999  continue
      ierr = nf_close(o3dmet)
      ierr = nf_close(o2dmet)
      if (kcmaq.ne.0) ierr = nf_close(kcmaq)
      if (kmyj.ne.0)  ierr = nf_close(kmyj)
      if (kysu.ne.0)  ierr = nf_close(kysu)
!
      stop
      end
!
!-----Date functions
!
      integer function addday(idate)
      implicit none
      integer idate,iyr,idy
      iyr = idate/1000
      idy = idate - iyr*1000
      if ((mod(iyr,4).eq.0 .and. idy.eq.366) .or. &
          (mod(iyr,4).ne.0 .and. idy.eq.365)) then
        iyr = iyr + 1
        if (iyr.gt.99) iyr = 0
        addday = iyr*1000 + 1
      else
        addday = idate + 1
      endif
      end
!
      integer function subday(idate)
      implicit none
      integer idate,iyr,idy
      iyr = idate/1000
      idy = idate - iyr*1000
      if (idy.eq.1) then
        iyr = iyr - 1
        if (iyr.lt.0) iyr = 99
        if (mod(iyr,4).eq.0) then
          idy = 366
        else
          idy = 365
        endif
        subday = iyr*1000 + idy
      else
        subday = idate - 1
      endif
      end
