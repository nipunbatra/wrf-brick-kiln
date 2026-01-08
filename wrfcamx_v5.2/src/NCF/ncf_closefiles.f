      subroutine ncf_closefiles(iounit)
c     use ncf_iomod
c
c-----This routine closes a NetCDF file.
c
c       Inputs:
c           iounit I  NetCDF file ID of opened file
c       Outputs:
c
c     include 'netcdf.inc'
c
      integer iounit
c
      integer ierr
c
c-----Entry point:
c
      ierr = nf_close(iounit)
c
      return
      end
