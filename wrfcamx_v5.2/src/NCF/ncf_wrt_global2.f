      subroutine ncf_wrt_global2(action,iounit,l2dsrf,l2dmet,ncf_name,nvars1,
     &                           nvars2,varname1,varname2)
      implicit none
c
c-----This routine writes the Global attributes to the NetCDF file
c
c      Argument description:
c       Inputs:
c           action    C  action string
c           iounit    I  NetCDF file ID of file
c           l2dsrf    L  2D surface file flag
c           l2dsrf    L  2D met file flag
c           ncf_name  C  file label
c           nvars1    I  number of species in list 1
c           nvars2    I  number of species in list 2
c           varname1  C  name of each variable in list 1
c           varname2  C  name of each variable in list 2
c       Outputs:
c
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       nvars1,nvars2
      character*(*) varname1(nvars1),varname2(nvars2)
      character*(*) ncf_name
      logical       l2dsrf
      logical       l2dmet
c
      character*(40*16) string
      integer           ierr, i
      integer           istrln
c
c-----Entry point:
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'FTYPE', NF_INT,
     &                                                    1, ncf_ftype)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CDATE', NF_INT,
     &                                                    1, ncf_cdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CTIME', NF_INT,
     &                                                    1, ncf_ctime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'WDATE', NF_INT,
     &                                                    1, ncf_wdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'WTIME', NF_INT,
     &                                                    1, ncf_wtime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'SDATE', NF_INT,
     &                                                    1, ncf_sdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'SDATEC', NF_INT,
     &                                                    1, ncf_sdatec)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'STIME', NF_INT,
     &                                                    1, ncf_stime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'TSTEP', NF_INT,
     &                                                    1, ncf_tstep)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NTHIK', NF_INT,
     &                                                 1, ncf_nthik)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NCOLS', NF_INT,
     &                                                 1, ncf_ncols)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NROWS', NF_INT,
     &                                                 1, ncf_nrows)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NCOLS_BUF', NF_INT,
     &                                                 1, ncf_ncolsbuf)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NROWS_BUF', NF_INT,
     &                                                 1, ncf_nrowsbuf)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      if (l2dsrf .or. l2dmet) then
        ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NLAYS', NF_INT,
     &                                                 1, 1)
      else
        ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NLAYS', NF_INT,
     &                                                 1, ncf_nlays)
      endif
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NVARS', NF_INT,
     &                                               1, nvars1+nvars2)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GDTYP', NF_INT,
     &                                                   1, ncf_gdtyp)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_ALP', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_alp) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
       ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_BET', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_bet) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_GAM', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_gam) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XCENT', NF_DOUBLE,
     &                                           1, DBLE(ncf_xcent) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YCENT', NF_DOUBLE,
     &                                           1, DBLE(ncf_ycent) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XORIG', NF_DOUBLE, 
     &                                           1, DBLE(ncf_xorig) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YORIG', NF_DOUBLE,
     &                                           1, DBLE(ncf_yorig) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XORIG_BUF', NF_DOUBLE, 
     &                                           1, DBLE(ncf_xorigbuf) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YORIG_BUF', NF_DOUBLE,
     &                                           1, DBLE(ncf_yorigbuf) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XCELL', NF_DOUBLE,
     &                                           1, DBLE(ncf_xcell) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YCELL', NF_DOUBLE,
     &                                           1, DBLE(ncf_ycell) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'VGTYP', NF_INT,
     &                                                    1, ncf_vgtyp)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_real(iounit, NF_GLOBAL, 'VGTOP', NF_FLOAT, 
     &                                                    1, ncf_vgtop)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'VGLVLS', NF_INT,
     &                                                    1, ncf_vglvls)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'GDNAM',
     &                       istrln(ncf_gdnam), trim(ncf_gdnam) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'UPNAM', 
     &                       istrln(ncf_upnam), trim(ncf_upnam) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      string = varname1(1)
      do i=2,nvars1
        string = string(1:(i-1)*16) // trim(varname1(i))
      enddo
      do i=1,nvars2
        string = string(1:(nvars1+i-1)*16) // trim(varname2(i))
      enddo
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'VAR-LIST',
     &             istrln(string)+9, trim(string)//'         ' )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'FILEDESC', 
     &                       istrln(ncf_name), trim(ncf_name) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'HISTORY', 
     &                       istrln(ncf_history), trim(ncf_history) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'IUTM', NF_INT,
     &                                                   1, ncf_iutm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'ISTAG', NF_INT,
     &                                                   1, ncf_istag)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CPROJ', NF_INT,
     &                                                   1, ncf_cproj)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      if (l2dsrf) then
        ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NSTEPS', NF_INT,
     &                                                 1, 1)
      else
        ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NSTEPS', NF_INT,
     &                                                 1, ncf_nsteps)
      endif
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'CAMx_NAME',
     &                       istrln(ncf_name), trim(ncf_name) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'NOTE',
     &                       istrln(ncf_note), trim(ncf_note) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'ITZON', NF_INT,
     &                                                    1, ncf_itzon)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'UPDSC', 
     &                       istrln(ncf_upnam), trim(ncf_upnam) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_GLOBAL:'
      write(*,'(A)') trim(action)
      write(*,'(A)') 'Cannot write global atttributes to file.'
      stop
c
 9999 continue
      return
      end
