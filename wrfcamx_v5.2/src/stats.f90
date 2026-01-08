      subroutine stats(arrnm,arrin,nx,ny,nz)
!
!-----STATS calculates min/avg/max of input array by layer, and outputs
!     the results to standard output

      implicit none
!
      integer nx,ny,nz
      real arrin(nx,ny,nz)
      character*20 arrnm
!
      integer i,j,k,nsum
      real rmin,rmax,sum

      do k = 1,nz
        rmin = 1.e10
        rmax = -1.e10
        sum = 0.
        nsum = 0
        do j = 1,ny
          do i = 1,nx
            if (arrnm(1:1).eq.'U' .or. arrnm(1:1).eq.'V') then
              sum = sum + abs(arrin(i,j,k))
              rmax = amax1(rmax,abs(arrin(i,j,k)))
              rmin = amin1(rmin,abs(arrin(i,j,k)))
            elseif (arrnm(1:5).eq.'Cloud'  .or.    &
                    arrnm(1:6).eq.'Precip' .or.    &
                    arrnm(1:4).eq.'Snow'   .or.    &
                    arrnm(1:2).eq.'KF'     .or.    &
                    arrnm(1:4).eq.'Conv'   .or.    &
                    arrnm(1:4).eq.'CAPE') then
              if (arrin(i,j,k).ne.0.) then
                sum = sum + arrin(i,j,k)
                nsum = nsum + 1
                rmax = amax1(rmax,arrin(i,j,k))
                rmin = amin1(rmin,arrin(i,j,k))
              endif
            else
              sum = sum + arrin(i,j,k)
              rmax = amax1(rmax,arrin(i,j,k))
              rmin = amin1(rmin,arrin(i,j,k))
            endif
          enddo
        enddo
        if (arrnm(1:5).eq.'Cloud' .or. arrnm(1:6).eq.'Precip' .or. &
            arrnm(1:4).eq.'Snow' .or. &
            arrnm(1:4).eq.'KF C' .or. arrnm(1:4).eq.'KF P' .or. &
            arrnm(1:4).eq.'Conv' .or. arrnm(1:4).eq.'CAPE') then
          if (sum.ne.0.) write(*,'(a20,i3,3f15.2)')arrnm,k,rmin,sum/nsum,rmax
        elseif (arrnm(1:6).eq.'KF Ent' .or. arrnm(1:6).eq.'KF Det') then
          if (sum.ne.0.) write(*,'(a20,i3,3f15.2)')arrnm,k,1000.*rmin, &
                                               1000.*sum/nsum,1000.*rmax
        else
          write(*,'(a20,i3,3f15.2)')arrnm,k,rmin,sum/(nx*ny),rmax
        endif
      enddo
      return
      end
