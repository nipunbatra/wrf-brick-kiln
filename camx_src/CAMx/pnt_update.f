      subroutine pnt_update()
      use filunit
      use grid
      use chmstry
      use bndary
      use ptemiss
      use tracer
c 
c----CAMx v7.32 250801
c 
c     PNT_UPDATE updates indexing of points of point sources and
c     the grid locations. This is needed because sources may have been
c     added in probing tools files.
c                           
c      Copyright 1996 - 2025
c     Ramboll
c           
c     Modifications: 
c 
c     Input arguments: 
c             
c     Output arguments: 
c             
c     Routines Called: 
c        RDPTHDR
c        GEOPSP 
c             
c     Called by: 
c        STARTUP 
c 
      include 'camx.prm'
      include 'flags.inc'

      integer ifle, icount
c
      integer indxpt(MXPTSRC)
      integer count_zero, count_outside, count_edge
c
c-----Entry point
c
c-----Coarse grid
c
      if (llatlon) then
        do n = 1,nptsrc
          xstk(n,1) = xloc(n) - xorg
          ystk(n,1) = yloc(n) - yorg
        enddo
      else
        do n = 1,nptsrc
          xstk(n,1) = xloc(n)/1000. - xorg
          ystk(n,1) = yloc(n)/1000. - yorg
        enddo
      endif
c
      count_zero = 0
      count_outside = 0
      count_edge = 0
      do 30 n=1,nptsrc
        ii = 1 + FLOOR(xstk(n,1)/delx)
        jj = 1 + FLOOR(ystk(n,1)/dely)
        isrc(n,1) = ii
        jsrc(n,1) = jj
        indxpt(n) = 1
        if (abs(dstk(n)) .LT. 0.01 .OR. abs(hstk(n) .LT. 0.01 .OR.
     &                                     vstk(n) .LT. 0.01 ) ) then
          count_zero = count_zero + 1
       endif
c
c-----Note sources outside domain
c
        if (ii.lt.1 .or. ii.gt.ncol(1) .or. jj.lt.1 .or.
     &                               jj.gt.nrow(1)) then
          count_outside = count_outside + 1 
          indxpt(n) = 0
        elseif (ii.le.1 .or. ii.ge.ncol(1) .or.
     &                 jj.le.1 .or. jj.ge.nrow(1)) then
          count_edge = count_edge + 1 
          indxpt(n) = 0
        endif
 30   continue
c
c --- Write a summary of "suspect" sources ---
c
      if( count_zero .GT. 0 ) then
          write(idiag,'(2A,I10)') 'Number of sources found with ',
     &                                    'zero stack parameters: ',count_zero
          write(idiag,'(2A,/,3A20,/,3(13X,F5.3))') 'Zero stack parameters ',
     &                 'are replaced with defaults:','Height (m)','Diameter (m)',
     &                                        '   Exit Velo (m/s)',1.0, 0.1, 0.1
      endif
      if( count_outside .GT. 0 ) then
          write(idiag,'(/,2A,I10)') 'Number of sources found ',
     &                            'that are outside the domain: ',count_outside
          write(idiag,'(A)') 'These sources will be skipped.'
      endif
      if( count_edge .GT. 0 ) then
          write(idiag,'(/,2A,I10)') 'Number of sources found that ',
     &                            'are located in an edge cell: ',count_edge
          write(idiag,'(A)') 'These sources will be skipped.'
      endif
c
c-----Map source to coarse grid
c
      nsrc = 0
      do n=1,nptsrc
        if (indxpt(n).eq.1) then
          nsrc = nsrc + 1
          idsrc(nsrc,1) = n
          isrc(nsrc,1) = isrc(n,1) 
          jsrc(nsrc,1) = jsrc(n,1) 
        endif
      enddo
      nosrc(1) = nsrc
c
c-----Fine grids
c
      do igrd = 2,ngrid
        nsrc = 0
        do n = 1,nptsrc
          if (indxpt(n).eq.1) then
            xstk(n,igrd) = xstk(n,1) - (inst1(igrd) - 1)*delx
            ystk(n,igrd) = ystk(n,1) - (jnst1(igrd) - 1)*dely
            ii = 2 + FLOOR(xstk(n,igrd)/delx*FLOAT( meshold(igrd) ) )
            jj = 2 + FLOOR(ystk(n,igrd)/dely*FLOAT( meshold(igrd) ) )
c
            if (ii.gt.1 .and. ii.lt.ncol(igrd) .and.
     &          jj.gt.1 .and. jj.lt.nrow(igrd)) then
              nsrc = nsrc + 1
              idsrc(nsrc,igrd) = n
              isrc(nsrc,igrd) = ii
              jsrc(nsrc,igrd) = jj
            endif
          endif
          nosrc(igrd) = nsrc
        enddo
      enddo
c
      write(idiag,*)
      do igrd = 1,ngrid
        write(idiag,*) 'Point sources in grid #',igrd,' = ',nosrc(igrd)
      enddo
      write(idiag,*)
c
      return
      end
