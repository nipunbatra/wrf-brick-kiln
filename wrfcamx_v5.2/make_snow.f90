program make_snow_multi
  implicit none
  integer, parameter :: nx=102, ny=102  ! CAMx grid + 2 for nested buffer
  real :: snowage(nx,ny)
  integer :: iunit
  integer :: snodate_jul, snotime  ! Internal YYJJJ format
  integer :: filedate              ! YYMMDD for filename
  integer :: file_hour             ! Hour for this file
  integer :: nhours, hour_idx
  integer :: nday(12), iyear, imonth, iday, jday
  character(len=30) :: fname

  data nday/31,28,31,30,31,30,31,31,30,31,30,31/

  snowage = 36.0

  ! Generate snow files that match CAMx start times:
  ! - 240131: needs 24032 1200 (for first day starting at noon)
  ! - 240201: needs 24033 0000 (for second day starting at midnight)
  ! - etc.

  ! File 240131 - for CAMx start 240201 12:00 (Julian 24032 1200)
  filedate = 240131
  snodate_jul = 24032
  snotime = 1200
  write(fname,'("snow.file.",I6.6)') filedate
  open(newunit=iunit, file=fname, form='unformatted', &
       access='sequential', status='replace')
  write(iunit) snodate_jul, snotime, nx, ny
  write(iunit) snowage
  close(iunit)
  print *, 'Generated: ', trim(fname), ' with Julian date ', snodate_jul, ' time ', snotime

  ! Files 240201-240210 - for CAMx start at 00:00 next day
  do iday = 1, 10
     filedate = 240200 + iday  ! 240201, 240202, ...
     snodate_jul = 24032 + iday  ! 24033, 24034, ... (next day Julian)
     snotime = 0
     write(fname,'("snow.file.",I6.6)') filedate
     open(newunit=iunit, file=fname, form='unformatted', &
          access='sequential', status='replace')
     write(iunit) snodate_jul, snotime, nx, ny
     write(iunit) snowage
     close(iunit)
     print *, 'Generated: ', trim(fname), ' with Julian date ', snodate_jul, ' time ', snotime
  end do

  print *, 'Done generating snow files'
end program make_snow_multi
