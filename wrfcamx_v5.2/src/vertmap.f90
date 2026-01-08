      subroutine vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,xin,xout)
!
!-----VERTMAP vertically aggregates WRF data on the met model grid
!     to the CAMx grid 
! 
!     NOTE: the CAMx physical height grid is a coarser set of the 
!           WRF eta coordinate system  
!
      implicit none
!
      integer nxc,nyc,nzc,nz,kz1(nzc),kz2(nzc)
      real xin(nxc,nyc,nz),xout(nxc,nyc,nzc)
      real eta(0:nz)
!
      integer i,j,k,kk
      real sum,deta
!
      do j = 1,nyc
        do i = 1,nxc
          do k = 1,nzc
            sum = 0.
            do kk = kz1(k),kz2(k) 
              deta = eta(kk-1) - eta(kk)
              sum = sum + xin(i,j,kk)*deta
            enddo
            deta = eta(kz1(k)-1) - eta(kz2(k))
            xout(i,j,k) = sum/deta
          enddo
        enddo
      enddo
!
      return
      end
