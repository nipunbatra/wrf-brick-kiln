      subroutine juldate(idate)
! 
!-----JULDATE converts date from calender (YYMMDD) format to Julian
!     (YYJJJ) format
!
      implicit none
!
      integer idate
      integer nday(12),iyear,imonth,iday,mday,n,jday
!
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
!
!-----Entry point
!
      iyear = idate/10000
      imonth = (idate - iyear*10000)/100
      iday = idate - iyear*10000 - imonth*100
!
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 n = 1,imonth-1
        mday = mday + nday(n)
 10   continue
      jday = mday + iday
      idate = iyear*1000 + jday
!
      return
      end

!-----------------------------------------------------------------------

      subroutine getime( this_date, this_time )
!
!-----This routine returns the current date/time as integer values in 
!     IOAPI format.
!
!   Arguments:
!     Outputs:
!       this_date    I   date (YYYYJJJ)
!       this_time    I   time (HHMMSS)
!
      integer this_date
      integer this_time
!
      character*8  cdate
      character*10 ctime
      character*5  czone
      integer*4    imon, iday, iyear, ihr, imin, isec
      integer*4    iarray(8)
!
!-----Entry point:
!
!   ---- call routine to get date/time as integers ---
!
      call date_and_time(cdate, ctime, czone, iarray) 
      iyear = iarray(1)
      imon = iarray(2)
      iday = iarray(3)
      ihr = iarray(5)
      imin = iarray(6)
      isec = iarray(7)
      iyear = MOD(iyear,100)
      this_date = iyear*10000 + imon*100 + iday
      call juldate(this_date)
      this_date = 2000000 + iyear*1000 + MOD(this_date,1000)
      this_time = ihr*10000 + imin*100 + isec
!
      return
      end

!-----------------------------------------------------------------------

      function istrln( string )
      integer   istrln
!
!-----This routine returns the non-blank length of a string.
!
!   Arguments:
!     Inputs:
!       string   C   string for determining length
!
      character*(*) string
!
      integer   i
!
!-----Entry point:
!
      istrln = 0
      do 10 i=LEN( string ),1,-1
         if( string(i:i) .NE. ' ' ) then
             istrln = i
             goto 9999
         endif
   10 continue
!
 9999 return
      end
