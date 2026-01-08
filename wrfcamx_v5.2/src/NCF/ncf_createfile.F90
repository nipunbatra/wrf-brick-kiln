      subroutine ncf_createfile(fname,action,iounit)
      implicit none
!
!-----This routine opens a NetCDF file.
!
!      Argument description:
!       Inputs:
!           fname  C  name of file to open
!           action C  description of exactly what this call is doing
!       Outputs:
!           iounit I  NetCDF file ID of opened file
!
      include 'netcdf.inc'
!
      character*(*) fname
      character*(*) action
      integer       iounit
!
      integer ierr, ncf_format
!
!-----Entry point:
!
      ncf_format = NF_CLOBBER
      ierr = nf_create(fname, ncf_format, iounit)      
      if( ierr .NE. nf_noerr ) then
         write(*,'(//,a)') 'ERROR in NCF_CREATEFILE:'
         write(*,'(A)') trim(action)
         write(*,'(2A)') 'Could not open file: ', trim(fname)
         stop
      endif
!
      return
      end
