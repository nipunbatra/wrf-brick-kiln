      subroutine ncf_enddef_file(action,iounit)
      implicit none
c
c-----This routine ends the definition mode of a NetCDF file
c
c      Argument description:
c       Inputs:
c           action  C  name of file to open
c           iounit  I  NetCDF file ID of file
c       Outputs:
c
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
c
      integer ierr
c
c-----Entry point:
c
      ierr = nf_enddef(iounit)
      if( ierr .NE. NF_NOERR ) then
        write(*,'(//,a)') 'ERROR in NCF_ENDDEF_FILE:'
        write(*,'(A)') trim(action)
        write(*,'(A)') 'Cannot terminate definition mode.'
        stop
      endif
c
      return
      end
