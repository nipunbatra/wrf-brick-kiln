      subroutine ncf_wrt_dim(action,iounit,numcols,numrows,nlays,nvars)
      implicit none
c
c-----This routine writes the dimensions to the NetCDF file
c
c      Argument description:
c       Inputs:
c           action  C  name of file to open
c           iounit  I  NetCDF file ID of file
c           numcols I number of cols in this file
c           numrows I number of cols in this file
c           nlays   I  number of layers in this file
c           nvars   I  number of species in the file
c       Outputs:
c
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       numcols
      integer       numrows
      integer       nlays
      integer       nvars
c
      integer ierr
c
c-----Entry point:
c
      ncf_date_time = 2
      ncf_lay = nlays
      ncf_col = numcols
      ncf_row = numrows
      ncf_var = nvars
c
      ierr = nf_def_dim(iounit, "TSTEP", NF_UNLIMITED, ncf_tstep_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "DATE-TIME", ncf_date_time, ncf_date_time_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "LAY", ncf_lay, ncf_lay_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "VAR", ncf_var, ncf_var_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "ROW", ncf_row, ncf_row_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "COL", ncf_col, ncf_col_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_DIM:'
      write(*,'(A)') trim(action)
      write(*,'(A)') 'Cannot write dimensions to file.'
      stop
c
 9999 continue
      return
      end
