      subroutine ncf_set_tstep(begin_date,begin_time,
     &                         ending_date,ending_time,dtout)
      use ncf_iomod
      implicit none
c
c-----This routine sets the time variable attributes for the NetCDF file
c
c      Argument description:
c       Inputs:
c         begin_date  I model begin date (YYJJJ)
c         begin_time  I model begin time
c         ending_date I model end date (YYJJJ)
c         ending_time I model end time
c         dtout       I output interval (min)
c       Outputs:
c
      integer      begin_date
      integer      begin_time
      integer      ending_date
      integer      ending_time
      integer      dtout
c
      integer date_now, time_now, num_steps
      integer mn,hr
      integer addday
c
c-----Entry point:
c
c  --- initialize to start of simulation ---
c
      ncf_nsteps = 0
      date_now = begin_date
      time_now = begin_time
c
c  --- walk through time to figure out how many steps for allocation ---
c
 111  continue
      ncf_nsteps = ncf_nsteps + 1
      mn = mod(time_now,100)
      hr = time_now/100
      if (dtout.lt.60) then
        mn = mn + dtout
        if (mn.gt.59) then
          mn = mn - 60
          hr = hr + 1
        endif
      else
        hr = hr + dtout/60
      endif
      if (hr.gt.23) then
        hr = hr - 24
        date_now = addday(date_now)
      endif
      time_now = 100*hr + mn

      if ((date_now.eq.ending_date .and. time_now.lt.ending_time) .or.
     &     date_now.lt.ending_date) goto 111
c
c  --- call routine to allocate allocate the array ----
c
      ncf_nsteps = ncf_nsteps + 1
      call ncf_alloc_tstep(ncf_nsteps)
c
c  --- now walk through time again to fill the NetCDF variable ---
c
      num_steps = 0
      date_now = begin_date
      time_now = begin_time
c
 222  continue
      num_steps = num_steps + 1
      if( num_steps .GT. ncf_nsteps ) goto 9999
      if (date_now/1000.gt.80) then
        ncf_tflag(1,num_steps) = 1900000 + date_now
        ncf_etflag(1,num_steps) = 1900000 + date_now
      else
        ncf_tflag(1,num_steps) = 2000000 + date_now
        ncf_etflag(1,num_steps) = 2000000 + date_now
      endif
      ncf_tflag(2,num_steps) = time_now*100
      ncf_etflag(2,num_steps) = time_now*100
c
      mn = mod(time_now,100)
      hr = time_now/100
      if (dtout.lt.60) then
        mn = mn + dtout
        if (mn.gt.59) then
          mn = mn - 60
          hr = hr + 1
        endif
      else
        hr = hr + dtout/60
      endif
      if (hr.gt.23) then
        hr = hr - 24
        date_now = addday(date_now)
      endif
      time_now = 100*hr + mn
      goto 222
c
 9999 continue
      return
      end
