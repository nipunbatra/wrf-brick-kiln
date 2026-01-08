      subroutine wrf_window(kz1,kz2,ioffset,joffset,nlu,kcmaq, &
                            kmyj,kysu,ldiag)
!
!-----WRF_WINDOW horizontally assigns, then vertically aggregates,
!     WRF data on the met model grid to the CAMx grid
!
!     NOTE: CAMx grid is simply a windowed area of the WRF grid, and the CAMx 
!           physical height grid is a coarser set of the WRF eta coordinate
!           system. 
!
      USE wrf_fields
      USE camx_fields

      implicit none

      integer nlu,kcmaq,kmyj,kysu
      integer kz1(nzc),kz2(nzc),ioffset,joffset

      integer i,j,k,n,iiu,jju,iiv,jjv,iix,jjx
      real mutu, mutv
      logical ldiag
!
!-----Couple WRF variables to P*
!
      do k = 1,nz
        do j = 1,ny
          ua(1,j,k) = ua(1,j,k)*mut(1,j,k)
          ua(nx+1,j,k) = ua(nx+1,j,k)*mut(nx,j,k)
          do i = 2,nx
            mutu = (mut(i,j,k) + mut(i-1,j,k))/2.
            ua(i,j,k) = ua(i,j,k)*mutu
          enddo
        enddo
        do i = 1,nx
          va(i,1,k) = va(i,1,k)*mut(i,1,k)
          va(i,ny+1,k) = va(i,ny+1,k)*mut(i,ny,k)
          do j = 2,ny
            mutv = (mut(i,j,k) + mut(i,j-1,k))/2.
            va(i,j,k) = va(i,j,k)*mutv
          enddo
        enddo
        do j = 1,ny
          do i = 1,nx
            ta(i,j,k) = ta(i,j,k)*mut(i,j,k)
            pa(i,j,k) = pa(i,j,k)*mut(i,j,k)
            qa(i,j,k) = qa(i,j,k)*mut(i,j,k)
          enddo
        enddo
      enddo
!
!-----Simple case of windowing the WRF grid with 1:1 meshing
!
      do k = 1,nz
        do j = 1,nyc
          do i = 1,nxc
            iiu = ioffset + i + 1
            jju = joffset + j
            utmp(i,j,k) = ua(iiu,jju,k)
            iiv = ioffset + i
            jjv = joffset + j + 1
            vtmp(i,j,k) = va(iiv,jjv,k)
            iix = i + ioffset
            jjx = j + joffset
            mutmp(i,j,k) = mut(iix,jjx,k)
            ttmp(i,j,k) = ta(iix,jjx,k)
            ptmp(i,j,k) = pa(iix,jjx,k)
            qtmp(i,j,k) = qa(iix,jjx,k)
            cwtmp(i,j,k) = qc(iix,jjx,k)
            prtmp(i,j,k) = qpr(iix,jjx,k)
            pstmp(i,j,k) = qps(iix,jjx,k)
            pgtmp(i,j,k) = qpg(iix,jjx,k)
            odtmp(i,j,k) = tau(iix,jjx,k)
            kwtmp(i,j,k) = kw(iix,jjx,k)
            kptmp(i,j,k) = kp(iix,jjx,k)
            ketmp(i,j,k) = ke(iix,jjx,k)
            kdtmp(i,j,k) = kd(iix,jjx,k)
!
!-----      Note: TKE is read in starting at the surface.  Store TKE
!           beginning at the top of the first layer
!
            if (kmyj.ne.0) then
              if (k .lt. nz) then
                tktmp(i,j,k) = tke(iix,jjx,k+1)
                eltmp(i,j,k) = el(iix,jjx,k+1)
              else
                tktmp(i,j,k) = 0.0
                eltmp(i,j,k) = 0.0
              endif
            endif
            ztmp(i,j,k) = zh(iix,jjx,k)
            if (k.eq.1) then
              tsfc(i,j)  = tsrf(iix,jjx)
              topcx(i,j) = topo(iix,jjx)
              snow(i,j)  = snocvr(iix,jjx)
              pblc(i,j)  = pbl(iix,jjx)
              kftmp(i,j) = cfr(iix,jjx)
              tctmp(i,j) = tc(iix,jjx)
              if (ldiag) then
                u10(i,j) = u10in(iix,jjx)
                v10(i,j) = v10in(iix,jjx)
                t2(i,j)  = t2in(iix,jjx)
                swsfc(i,j) = swd(iix,jjx)
                soilm(i,j) = smin(iix,jjx)
                soilt(i,j) = stin(iix,jjx)
                ctop(i,j)  = ctopw(iix,jjx)
                cape(i,j)  = capew(iix,jjx)
                styp(i,j)  = stypin(iix,jjx)
                cldfrc(i,j) = cfin(iix,jjx)
              endif
              do n = 1,nlu
                lucx(i,j,n)  = clu(iix,jjx,n)
              enddo
              if (kcmaq.ne.0 .or. kysu.ne.0) then
                z0c(i,j) = z0(iix,jjx)
              endif
            endif
          enddo
        enddo
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
             if (i.eq.nxc) then
               uac(i,j,k) = uac(i,j,k)/mutc(i,j,k)
             else
               mutu = (mutc(i,j,k) + mutc(i+1,j,k))/2.
               uac(i,j,k) = uac(i,j,k)/mutu
             endif
             if (j.eq.nyc) then
               vac(i,j,k) = vac(i,j,k)/mutc(i,j,k)
             else
               mutv = (mutc(i,j,k) + mutc(i,j+1,k))/2.
               vac(i,j,k) = vac(i,j,k)/mutv
             endif
           enddo
         enddo
         do j = 1,nyc
           do i = 1,nxc
             tac(i,j,k) = tac(i,j,k)/mutc(i,j,k)
             pac(i,j,k) = pac(i,j,k)/mutc(i,j,k)
             qac(i,j,k) = qac(i,j,k)/mutc(i,j,k)
           enddo
         enddo
       enddo
!
!-----Map layer interface heights, TKE, and mixing lengths to CAMx
!     vertical grid
!
      do j = 1,nyc 
        do i = 1,nxc 
          do k = 1,nzc 
            zhc(i,j,k) = ztmp(i,j,kz2(k))
            if (kmyj.ne.0) then
              tkc(i,j,k) = tktmp(i,j,kz2(k))
              elc(i,j,k) = eltmp(i,j,kz2(k))
            endif
          enddo 
        enddo 
      enddo 
!
      return
      end
