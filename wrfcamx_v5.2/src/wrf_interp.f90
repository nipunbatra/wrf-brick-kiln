      subroutine wrf_interp(project,kz1,kz2,nlu,deltax,kcmaq, &
                            kmyj,kysu,ldiag)
!
!-----WRF_INTERP horizontally interpolates, then vertically aggregates,
!     WRF data on the met model grid to a different CAMx grid
!
!     NOTE: CAMx grid must be wholly contained within the WRF data grid
!           limits, and the CAMx physical height grid is a coarser set
!           of the WRF eta coordinate system
!
!     NOTE: This option puts all variables from WRF Arakawa C split to
!           cell centers on the CAMx grid
!
      USE wrf_fields
      USE camx_fields

      implicit none
 
      integer nlu,kcmaq,kmyj,kysu
      integer kz1(nzc),kz2(nzc)
      real deltax
      character*10 project

      integer i,j,k,n
      real rotate,ws,wd,wdc
      real*8 rearth,pi,r2d,cfact,ax,ay,bx,by,cx,cy,distac,distbc,arg1
      real*8 dbl_one
      logical ldiag
      dbl_one = -1.0
      rearth = 6370.
      pi     = dacos(dbl_one)
      r2d    = 180./pi
      cfact  = 1.
      if (project.eq.'LATLON    ') cfact = r2d
!
!-----First interpolate Arakawa C staggered winds to WRF cell centers
!
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ucrs(i,j,k) = (ua(i+1,j,k) + ua(i,j,k))/2.
            vcrs(i,j,k) = (va(i,j+1,k) + va(i,j,k))/2.
          enddo
        enddo
      enddo
!
!-----Rotate winds from WRF projection to CAMx grid projection:
!     Take 2 points with the same WRF i value.  If the CAMx i coordinates are
!     identical, the WRF grid is lined up with the CAMx grid and no
!     angle rotation is needed. Otherwise, find the correction angle
!     using trig functions at the points A(xa,ya),B(xb,yb),C(xa,yb)
!
      do j = 1,ny
        do i = 1,nx
          ax = xdprj(i,j)/cfact
          ay = ydprj(i,j)/cfact
          bx = xdprj(i,j+1)/cfact
          by = ydprj(i,j+1)/cfact
          cx = ax
          cy = by
          if (project.eq.'LATLON    ') then
            arg1 = dcos(cx-ax)*dcos(cy)*dcos(ay) + dsin(cy)*dsin(ay)
            arg1 = dmin1(dble(1.),dmax1(dble(-1.),arg1))
            distac = rearth*dacos(arg1)
            arg1 = dcos(cx-bx)*dcos(cy)*dcos(by) + dsin(cy)*dsin(by)
            arg1 = dmin1(dble(1.),dmax1(dble(-1.),arg1))
            distbc = rearth*dacos(arg1)
          else
            distac = sqrt((cx-ax)**2 + (cy-ay)**2)
            distbc = sqrt((cx-bx)**2 + (cy-by)**2)
          endif
          rotate = r2d*datan(distbc/distac)
          if (cx .gt. bx) rotate = -rotate
!
!-----Compute wind speed and direction from the WRF u and v components
!
          do k = 1,nz
            ws = sqrt((ucrs(i,j,k)**2) + (vcrs(i,j,k)**2))
            if (ucrs(i,j,k) .gt. 0.) then
              wd = 270. - atan(vcrs(i,j,k)/ucrs(i,j,k))*r2d
            elseif (ucrs(i,j,k) .lt. 0.) then
              wd = 90. - atan(vcrs(i,j,k)/ucrs(i,j,k))*r2d
            elseif ((ucrs(i,j,k) .eq. 0.) .and. (vcrs(i,j,k) .gt. 0.)) then
              wd = 180.
            else
              wd = 0.
            endif
!
!-----Add grid rotation to wind direction and convert wind back to 
!     u and v components
!
            wdc = wd + rotate
            ucrs(i,j,k) = -ws*sin(wdc/r2d)
            vcrs(i,j,k) = -ws*cos(wdc/r2d)
          enddo
        enddo
      enddo
!
!-----Couple WRF variables to P*
!
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ucrs(i,j,k) = ucrs(i,j,k)*mut(i,j,k)
            vcrs(i,j,k) = vcrs(i,j,k)*mut(i,j,k)
            ta(i,j,k) = ta(i,j,k)*mut(i,j,k)
            pa(i,j,k) = pa(i,j,k)*mut(i,j,k)
            qa(i,j,k) = qa(i,j,k)*mut(i,j,k)
          enddo
        enddo
      enddo
!
!-----Interpolate WRF data onto the CAMx grid (all data to cell centers)
!
      do k = 1,nz
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,mut(1,1,k),mutmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,ucrs(1,1,k),utmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,vcrs(1,1,k),vtmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,ta(1,1,k),ttmp(1,1,k),xcrs,ycrs, &
                      xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,pa(1,1,k),ptmp(1,1,k),xcrs,ycrs, &
                      xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,qa(1,1,k),qtmp(1,1,k),xcrs,ycrs, &
                      xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,zh(1,1,k),ztmp(1,1,k),xcrs,ycrs, &
                      xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,qc(1,1,k),cwtmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,qpr(1,1,k),prtmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,qps(1,1,k),pstmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,qpg(1,1,k),pgtmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,tau(1,1,k),odtmp(1,1,k),xcrs, &
                      ycrs,xcwrf,ycwrf,deltax)
        if (kmyj.ne.0 .and. k.lt.nz) then
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,tke(1,1,k+1),tktmp(1,1,k), &
                        xcrs,ycrs,xcwrf,ycwrf,deltax)
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,el(1,1,k+1),eltmp(1,1,k),xcrs, &
                        ycrs,xcwrf,ycwrf,deltax)
        endif
        if ((kcmaq.ne.0 .or. kysu.ne.0) .and. k.eq.1) then
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,z0,z0c,xcrs,ycrs, &
                        xcwrf,ycwrf,deltax)
        endif
        if (k.eq.1) then
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,pbl,pblc,xcrs,ycrs, &
                        xcwrf,ycwrf,deltax)
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,tsrf,tsfc,xcrs,ycrs, &
                        xcwrf,ycwrf,deltax)
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,topo,topcx,xcrs,ycrs, &
                        xcwrf,ycwrf,deltax)
          call interp2d(nx,ny,nxc,nyc,icrs,jcrs,snocvr,snow,xcrs,ycrs, &
                        xcwrf,ycwrf,deltax)
          if (ldiag) then
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,u10in,u10,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,v10in,v10,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,t2in,t2,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,swd,swsfc,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,smin,soilm,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,stin,soilt,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,ctopw,ctop,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,capew,cape,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,stypin,styp,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
            call interp2d(nx,ny,nxc,nyc,icrs,jcrs,cfin,cldfrc,xcrs,ycrs, &
                          xcwrf,ycwrf,deltax)
          endif
        endif
      enddo
      do n = 1,nlu
        call interp2d(nx,ny,nxc,nyc,icrs,jcrs,clu(1,1,n),lucx(1,1,n), &
                      xcrs,ycrs,xcwrf,ycwrf,deltax)
      enddo
!
!-----Map momentum and thermodynamic variables onto the CAMx vertical
!     grid structure
!
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,mutmp,mutc)
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,utmp,uac)
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,vtmp,vac)
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,ttmp,tac)
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,ptmp,pac)
      call vertmap(nxc,nyc,nzc,nz,kz1,kz2,eta,qtmp,qac)
!
!-----Decouple vertically interpolated variables from P*
!
      do k = 1,nzc
        do j = 1,nyc
          do i = 1,nxc
            uac(i,j,k) = uac(i,j,k)/mutc(i,j,k)
            vac(i,j,k) = vac(i,j,k)/mutc(i,j,k)
            tac(i,j,k) = tac(i,j,k)/mutc(i,j,k)
            pac(i,j,k) = pac(i,j,k)/mutc(i,j,k)
            qac(i,j,k) = qac(i,j,k)/mutc(i,j,k)
          enddo
        enddo
      enddo
!
!-----Map layer interface heights and TKE to CAMx vertical grid
!
      do j = 1,nyc 
        do i = 1,nxc 
          do k = 1,nzc 
            zhc(i,j,k) = ztmp(i,j,kz2(k))
            if (kmyj.ne.0) then
              if (kz2(k) .eq. nz) then
                tkc(i,j,k) = 0.0
                elc(i,j,k) = 0.0
              else
                tkc(i,j,k) = tktmp(i,j,kz2(k))
                elc(i,j,k) = eltmp(i,j,kz2(k))
              endif
            endif
          enddo 
        enddo 
      enddo 

      return
      end
!
!-------------------------------------------------------------------------
!
      subroutine interp2d(nx,ny,nxc,nyc,iwrf,jwrf,valm,valc,xm,ym,xc,yc,dx)
!
!-----INTERP2D horizontally interpolates an input field to a different grid
!
      implicit none
      integer nx,ny,nxc,nyc,iwrf(nxc,nyc),jwrf(nxc,nyc)
      real valm(nx,ny),valc(nxc,nyc)
      real xm(nx),ym(ny),xc(nxc,nyc),yc(nxc,nyc),dx

      integer i,j,im,jm
      real dcdx1,dcdx2,c1,c2,dcdy

      do j = 1,nyc
        do i = 1,nxc
          im = iwrf(i,j)
          jm = jwrf(i,j)
          dcdx1 = (valm(im,jm-1) - valm(im-1,jm-1))/dx
          dcdx2 = (valm(im,jm) - valm(im-1,jm))/dx
          c1 = valm(im-1,jm-1) + dcdx1*(xc(i,j)-xm(im-1))
          c2 = valm(im-1,jm) + dcdx2*(xc(i,j)-xm(im-1))
          dcdy = (c2 - c1)/dx
          valc(i,j) = c1 + dcdy*(yc(i,j)-ym(jm-1))
        enddo
      enddo

      return
      end
