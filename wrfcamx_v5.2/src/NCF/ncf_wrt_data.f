      subroutine ncf_wrt_data(action,iounit,num_vars,num_cols,num_rows,num_lays,
     &                        ncf_cur_tstep,var_name,varfld)
      implicit none
c
c-----This routine writes the data for the time-varying variables to the
c     NetCDF file
c
c      Argument description:
c       Inputs:
c           action        C name of file to open
c           iounit        I NetCDF file ID of file
c           num_vars      I number of variables in the array
c           num_cols      I number of columns in this grid
c           num_rows      I number of rows in this grid
c           num_lays      I number of layers in the array
c           ncf_cur_tstep I time step index
c           var_name      C variable name
c           varfld        R gridded data array
c       Outputs:
c
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       num_vars
      integer       num_cols
      integer       num_rows
      integer       num_lays
      integer       ncf_cur_tstep
      character*(*) var_name(num_vars)
      real          varfld(num_cols,num_rows,num_lays,num_vars)      
c
      character*60  this_var
      integer       data_start(4), data_count(4), i, j, k
      integer       n, ierr, this_varid
      real,         allocatable, dimension(:,:,:) :: array_3d
c
c-----Entry point:
c
c  --- set the position in the NetCDF variable to write ---
c
      data_start(1) = 1
      data_count(1) = num_cols
      data_start(2) = 1
      data_count(2) = num_rows
      data_start(3) = 1
      data_count(3) = num_lays
      data_start(4) = ncf_cur_tstep
      data_count(4) = 1
c
c  --- allocate the array that will be used to write the data ---
c
      allocate( array_3d(num_cols,num_rows,num_lays) )
c
c  --- get name for this variable, and get it's variable ID
c      in this file ---
c
      do n = 1,num_vars
        this_var = var_name(n)
        ierr = nf_inq_varid(iounit,trim(this_var),this_varid)
        if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- load the data into the local array to write ---
c
        do k=1,num_lays
           do j=1,num_rows
              do i=1,num_cols
                 array_3d(i,j,k) = varfld(i,j,k,n)
              enddo
           enddo
        enddo
c
        ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
        if( ierr .NE. NF_NOERR ) goto 7001
      enddo
c
      deallocate( array_3d )
c
      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot find variable id for: ',trim(this_var)
      stop
c
 7001 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot write data for the variable: ',trim(this_var)
      stop
c
 9999 continue
      return
      end
