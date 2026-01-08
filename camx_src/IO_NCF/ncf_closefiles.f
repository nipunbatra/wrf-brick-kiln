c**** NCF_CLOSEFILES
c
      subroutine ncf_closefiles()
      use grid
      use pigsty
      use filunit
      use tracer
      use procan
      use ncf_iomod
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine closes a NetCDF file.
c
c      Copyright 1996 - 2025
c     Ramboll
c      Argument description:
c       Inputs:
c       Outputs:
c           iounit I  NetCDF file ID of opened file
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer iounit
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer igrd, ierr, icsa_grid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c
c  --- CSA is always NetCDF ---
c    
      if( lcsa ) then
        do icsa_grid=1,npadom
           ierr = nf_close(iow_csa(icsa_grid))
        enddo
      endif
c
c  --- if rest of model is not NetCDF output just leave ---
c
      if( .NOT. lcdfout ) goto  9999
c
      if( allocated(iavg) ) then
         do igrd=1,ngrid
           ierr = nf_close(iavg(igrd))
           ierr = nf_close(idep(igrd))
         enddo
      endif
c
      if( allocated(isample) ) then
         do igrd=1,nsample
           ierr = nf_close(isample(igrd))
         enddo
      endif
c
      if( allocated(iowsfc) ) then
         do igrd=1,ngrid
           ierr = nf_close(iowsfc(igrd))
         enddo
      endif
c
      if( allocated(iowptdep) ) then
         do igrd=1,ngrid
           ierr = nf_close(iowptdep(igrd))
         enddo
      endif
c
      if( allocated(ismout) ) then
         do igrd=1,ngrid
           ierr = nf_close(ismout(igrd))
         enddo
      endif
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
