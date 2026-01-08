      Module ncf_iomod
      include 'ncf_iodat.inc'
c
c-----This allocates the dynamic memory arrays in the NCF_IOMOD.INC
c     include file.
c
      Contains
c
       subroutine ncf_alloc_tstep(num_tsteps)
       implicit none
c
c-----This routine allocates the array to store the NetCDF time step
c     varaibles
c
c     Input:
c       num_steps I number of time steps in this simulation
c     Output:  
c
       integer :: num_tsteps
c
c----Entry point:
c
       if( .NOT. allocated(ncf_tflag) ) allocate( ncf_tflag(2,num_tsteps) )
       if( .NOT. allocated(ncf_etflag) ) allocate( ncf_etflag(2,num_tsteps) )
c
       return
       end subroutine
c
      end Module
