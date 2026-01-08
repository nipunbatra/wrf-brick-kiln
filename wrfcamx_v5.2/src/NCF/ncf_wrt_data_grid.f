      subroutine ncf_wrt_data_grid(action,iounit,grid_ncol,grid_nrow,
     &                             grid_nlay,grid_xorig,grid_yorig,grid_deltax,
     &                             grid_deltay,grid_lat,grid_lon)
      implicit none
c
c-----Description:
c
c   This routine writes the data for the grid definiton variables to the
c   NetCDF file
c
c      Argument description:
c       Inputs:
c           action      C name of file to open
c           iounit      I NetCDF file ID of file
c           grid_ncol   I number of grid cells in X direction
c           grid_nrow   I number of grid cells in Y direction
c           grid_nlay   I number of grid layers
c           grid_xorig  R X coordinate of grid origin      
c           grid_yorig  R Y coordinate of grid origin      
c           grid_deltax R cell width in X direction
c           grid_deltay R cell width in Y direction
c           grid_lat    R latitude value of cell center
c           grid_lon    R longitude value of cell center
c       Outputs:
c
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       grid_ncol
      integer       grid_nrow
      integer       grid_nlay
      real          grid_xorig
      real          grid_yorig
      real          grid_deltax
      real          grid_deltay
      real          grid_lat(grid_ncol,grid_nrow)
      real          grid_lon(grid_ncol,grid_nrow)
c
      character*200 this_var
      real*8,       allocatable, dimension(:)   :: darray_1d
      real*8,       allocatable, dimension(:,:) :: darray_2d
      integer,      allocatable, dimension(:)   :: iarray_1d
      integer       this_varid, icol, irow, ilay, ierr
c
c-----Entry point:
c
c  --- variable for X coordinates ---
c
      this_var = "X"
      allocate( darray_1d(grid_ncol) )
      do icol=1,grid_ncol
        darray_1d(icol) = DBLE(grid_xorig + (REAL(icol)-0.5)*grid_deltax)
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( darray_1d )
c
c  --- variable for Y coordinates ---
c
      this_var = "Y"
      allocate( darray_1d(grid_nrow) )
      do irow=1,grid_nrow
        darray_1d(irow) = DBLE(grid_yorig + (REAL(irow)-0.5)*grid_deltay)
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( darray_1d )
c
c  --- variable for Layer indices ---
c
      this_var = "layer"
      allocate( iarray_1d(grid_nlay) )
      do ilay=1,grid_nlay
        iarray_1d(ilay) = ilay
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_int(iounit,this_varid,iarray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( iarray_1d )
c
c  --- variable for longitude coordinates ---
c
      allocate( darray_2d(grid_ncol,grid_nrow) )
c
      this_var = "longitude"
      do icol=1,grid_ncol
        do irow=1,grid_nrow
           darray_2d(icol,irow) = DBLE(grid_lon(icol,irow))
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  --- variable for latitude coordinates ---
c
      this_var = "latitude"
      do icol=1,grid_ncol
        do irow=1,grid_nrow
           darray_2d(icol,irow) = DBLE(grid_lat(icol,irow))
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( darray_2d )
c
      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA_GRID:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot find variable id for: ',trim(this_var)
      stop
c
 7001 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DATA_GRID:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot write data for the variable: ',trim(this_var)
      stop
c
 9999 continue
      return
      end
