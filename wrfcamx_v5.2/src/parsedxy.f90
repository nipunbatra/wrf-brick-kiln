      subroutine parsedxy(inrec,dx,dy)

!-----Parse input record for Dx and Dy (grid cell size) to account for fractions
!     of km (e.g., expressed as 4/3).
!     Thanks to Alex Cohan, LADCO for this contribution.

      implicit none
      character*200 inrec
      integer j,i,j2,j3,ii
      real tt1,tt2
      real dx,dy

      j = 0
      i = 1
      ii = len_trim(inrec)

!-----Search for any instance of "/"

      do while (j.eq.0 .and. i.le.ii)
        if (inrec(i:i).eq.'/') j = i
        i = i + 1
      enddo

!-----No "/" found, read grid cell size directly

      if (j.eq.0) then
       read(inrec,*) dx,dy
       return
      endif

!-----Otherwise, search for intervening " " and second "/"

      j2 = 0
      do while (j2.eq.0 .and. i.le.ii)
        if (inrec(i:i) .eq. ' ') j2 = i
        i = i + 1
      enddo
      if (j2.eq.0) then
        write(*,*)'ERROR: In PARSEDXY, cannot find space separating', &
                  ' CAMx Dx and Dy entries'
        stop
      endif

      j3 = 0
      do while (j3.eq.0 .and. i.le.ii)
        if (inrec(i:i).eq.'/') j3 = i
        i = i + 1
      enddo
      if (j3.eq.0) then
        write(*,*)'ERROR: In PARSEDXY, cannot find two slashes', &
                  ' defining fractional CAMx Dx and Dy entries'
        stop
      endif

      read(inrec(1:(j-1)),*) dx
      read(inrec((j+1):(j2-1)),*) tt1
      read(inrec((j2+1):(j3-1)),*) dy
      read(inrec((j3+1):ii),*) tt2
      dx = dx/tt1
      dy = dy/tt2

      return
      end
