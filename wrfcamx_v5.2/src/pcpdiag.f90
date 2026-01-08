      subroutine pcpdiag(nz,lconpcp,precip,zh,ta,qc,qpr,qpg,qps)
!
!-----PCPDIAG calculates precipitation water profiles based on resolved
!     cloud inputs or from un-resolved cloud fields derived from CLDDIAG.
!
!     Input:
!        nz       number of layers
!        lconpcp  flag indicating convective precip is present
!        precip   total precip rate (mm/hr)
!        zh       Layer interface heights (m)
!        ta       Layer temperature (K)
!        qc       Layer cloud water content (g/m3)
!     
!     Output:  
!        qpr      Layer liquid precip content (g/m3)
!        qpg      Layer graupel precip content (g/m3)
!        qps      Layer snow precip content (g/m3)
!
      implicit none
!
!-----Argument variables
!
      integer nz
      real precip
      real zh(nz),ta(nz),qc(nz),qpr(nz),qpg(nz),qps(nz)
      logical lconpcp
!
!-----Local variables
!
      integer k,inuf,ktop
      real pi,dz,zmin,drpdia,drpmas,drpvel,cscav,cliq
      real pp(nz),qp(nz)
!
      pi = acos(-1.0)
      do k = 1,nz
        pp(k) = 0.
        qp(k) = 0.
        qpr(k) = 0.
        qpg(k) = 0.
        qps(k) = 0.
      enddo
      if (precip .lt. 0.2) goto 999
!
!-----Find top of precip profile while T > 233K
!
      dz = 0.
      inuf = 0
      zmin = 300. + (2700./25.)*precip
      ktop = 1
      do k = 1,nz
        if (qc(k) .gt. 0.01 .and. ta(k) .gt. 233.) then
          ktop = k
          if (k .eq. 1) then
            dz = zh(k)
          else
            dz = dz + (zh(k) - zh(k-1))
          endif
          if (dz .ge. zmin .and. inuf .eq. 0) inuf = k
        endif
      enddo
      if (inuf .eq. 0) then
        precip = 0.
        goto 999
!       dz = 0.
!       do k = 2,nz
!         qc(k) = max(qc(k),min(1.,1.e6*1.0e-7*precip**0.79))
!         dz = dz + (zh(k) - zh(k-1))
!         if (dz.ge.zmin) then
!           if (k.gt.ktop) ktop = k
!           goto 203
!         endif
!       enddo
!       ktop = nz
      endif
!
!-----Go up the column and define precip water based on presence
!     of cloud water
!
 203  pp(1) = precip
      qp(1) = 1.e6*1.0e-7*pp(1)**0.79
      do k = 2,ktop
        pp(k) = pp(k-1)
        qp(k) = qp(k-1)
      enddo
!
!-----Load precip water into resolved or convective arrays
!
      do k = 1,nz
        if (ta(k) .gt. 273.) then
          qpr(k) = qp(k)
        else
          if (lconpcp) then
            qpg(k) = qp(k) 
          else
            qps(k) = qp(k)
          endif
        endif
      enddo
!
 999  return
      end
