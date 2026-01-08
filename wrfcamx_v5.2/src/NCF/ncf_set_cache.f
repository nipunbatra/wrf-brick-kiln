      subroutine ncf_set_cache()
      implicit none
c
c-----This routine sets the cache parameters that will control how
c     the NetCDF output files will be chunked.
c
c      Argument description:
c       Inputs:
c       Outputs:
c
      include 'ncf_iodat.inc'
c
      integer primes(166)
      logical lfound_prime
      integer i
c
      data primes /2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
     &             71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,
     &             149,151,157,163,167,173,179,181,191,193,197,199,211,
     &             223,227,229,233,239,241,251,257,263,269,271,277,281,
     &             283,293,307,311,313,317,331,337,347,349,353,359,367,
     &             373,379,383,389,397,401,409,419,421,431,433,439,443,
     &             449 457,461,463,467,479,487,491,499,503,509,521,523,
     &             541,547,557,563,569,571,577,587,593,599,601,607,613,
     &             617,619,631,641,643,647,653,659,661,673,677,683,691,
     &             701,709,719,727,733,739,743,751,757,761,769,773,787,
     &             797,809,811,821,823,827,829,839,853,857,859,863,877,
     &             881,883,887,907,911,919,929,937,941 947,953,967,971,
     &             977,983,991,997/
c
c-----Entry point:
c
      if ( NCF_CACHESIZE .LE. 0 ) then
        write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(*,'(A)') 'Invalid value provided for parameter NCF_CACHESIZE.'
        write(*,'(2A)') 'Value represents total size of the raw data ',
     &                                           'chunk cache in MegaBytes.'
        write(*,'(A)') 'It must be a postive number.'
        stop
      endif
c
      lfound_prime = .FALSE.
      do i=1,166
         if( NCF_NELEMS .EQ. primes(i) ) lfound_prime = .TRUE.
      enddo
      if( NCF_NELEMS .LE. 0 .OR. .NOT. lfound_prime ) then
        write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(*,'(A)') 'Invalid value provided for parameter NCF_NELEMS.'
        write(*,'(2A)') 'This value must be a prime number.'
        stop
      endif
c
      if( NCF_PREEMPTION .LT. 0 .OR. NCF_PREEMPTION .GT. 100 ) then
        write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(*,'(A)') 'Invalid value provided for parameter NCF_PREEMPTION.'
        write(*,'(A)') 'This value must be between 0 and 100.'
        stop
      endif
c
      if( NCF_DEFLATE_LEVEL .LT. 0 .OR.  NCF_DEFLATE_LEVEL .GT. 9 ) then
        write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(*,'(A)') 'Invalid value provided for parameter NCF_DEFLATE_LEVEL.'
        write(*,'(A)') 'This value must be between 0 and 9.'
        stop
      endif
c
      return
      end
 
