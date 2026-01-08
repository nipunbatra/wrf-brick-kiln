      subroutine ncf_wrt_data_tstep(action,iounit,ncf_cur_tstep,nvars)
      use ncf_iomod
      implicit none
c
c-----This routine writes the data for the timestep variables to the
c     NetCDF file
c
c      Argument description:
c       Inputs:
c           action      C name of file to open
c           iounit      I NetCDF file ID of file
c         ncf_cur_tstep I time step index
c           nvars       I number of variables in the file
c       Outputs:
c
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       ncf_cur_tstep
      integer       nvars
c
      character*200 this_var
      integer,      allocatable, dimension(:,:) :: iarray_2d
      integer       this_varid, ivar, istep, ierr
      integer       data_start(3), data_count(3)
c
c-----Entry point:
c
       data_start(1) = 1
       data_count(1) = 2
       data_start(2) = 1
       data_count(2) = nvars
       data_start(3) = ncf_cur_tstep
       data_count(3) = 1
c
c  --- variable for TFLAG ---
c
      allocate( iarray_2d(2,nvars) )
      do ivar=1,nvars
        iarray_2d(1,ivar) = ncf_tflag(1,ncf_cur_tstep)
        iarray_2d(2,ivar) = ncf_tflag(2,ncf_cur_tstep)
      enddo
      this_var = "TFLAG"
      this_varid = 0
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_vara_int(iounit,this_varid,data_start,data_count,iarray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( iarray_2d )
c
c  --- variable for ETFLAG ---
c
      allocate( iarray_2d(2,nvars) )
      do ivar=1,nvars
        iarray_2d(1,ivar) = ncf_etflag(1,ncf_cur_tstep)
        iarray_2d(2,ivar) = ncf_etflag(2,ncf_cur_tstep)
      enddo
      this_var = "ETFLAG"
      this_varid = 0
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_vara_int(iounit,this_varid,data_start,data_count,iarray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( iarray_2d )

      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA_TSTEP:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot find variable id for: ',trim(this_var)
      stop
c
 7001 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA_TSTEP:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot write data for the variable: ',trim(this_var)
      stop
c
 9999 continue
      return
      end
