  MODULE wrf_fields
!
!-----WRF fields
!
  real, allocatable, dimension(:,:,:) :: ua,va,ucrs,vcrs,ta,qa,tke,el
  real, allocatable, dimension(:,:,:) :: qc,qpr,qps,qpg,tau
  real, allocatable, dimension(:,:,:) :: zh, ph, phb
  real, allocatable, dimension(:,:,:) :: pa, pab
  real, allocatable, dimension(:,:,:) :: ql,qi
  real, allocatable, dimension(:,:,:) :: kf,kw,kp,kd,ke
  real, allocatable, dimension(:,:,:) :: ksf,kdf,kl,ki,kud,kue,kdd,kde
  real, allocatable, dimension(:,:,:) :: clu,rluf
  real, allocatable, dimension(:,:,:) :: taold,paold,qvold,qcold
  real, allocatable, dimension(:,:,:) :: smois,stemp
  real, allocatable, dimension(:,:,:) :: mut

  real, allocatable, dimension(:,:)   :: xdprj,ydprj
  real, allocatable, dimension(:,:)   :: mu,mub,pbl,rlu
  real, allocatable, dimension(:,:)   :: tsrf,rainr,rainc,z0,ylat,xlon
  real, allocatable, dimension(:,:)   :: snocvr,topo
  real, allocatable, dimension(:,:)   :: rainro,rainco
  real, allocatable, dimension(:,:)   :: cfr,timec,tc
  real, allocatable, dimension(:,:)   :: smin,stin,stypin
  real, allocatable, dimension(:,:)   :: ctopw,capew,cfin
  real, allocatable, dimension(:,:)   :: si
  integer, allocatable, dimension(:,:) :: istin

  real, allocatable, dimension(:)     :: xcrs,ycrs,ydot
  real, allocatable, dimension(:)     :: eta,c1h,c2h
  real, allocatable, dimension(:)     :: zz,dz,tt,pr,qq,qv,cfrac,rh
  real, allocatable, dimension(:)     :: qr,qs,qg

  real, allocatable, dimension(:,:)   :: u10in,v10in,t2in,swd

  integer nx,ny,nz,nsrf
  integer ihyb
!
!--------------------------------------------------------------------
!
  CONTAINS

  subroutine alloc_wrf(nx,ny,nz)
  implicit none

! (nx) -- Max WEST-EAST STAGGERED
! (ny) -- Max SOUTH-NORTH STAGGERED
! (nz) -- Max BOTTOM-TOP STAGGERED
! (nx-1) -- Max WEST-EAST UNSTAGGERED
! (ny-1) -- Max SOUTH-NORTH UNSTAGGERED
! (nz-1) -- Max BOTTOM-TOP UNSTAGGERED
! (nsrf) -- Number of soil layers

  integer nx,ny,nz

     allocate( ua(nx,ny-1,nz-1)      )
     allocate( va(nx-1,ny,nz-1)      )
     allocate( ucrs(nx-1,ny-1,nz-1)  )
     allocate( vcrs(nx-1,ny-1,nz-1)  )
     allocate( ta(nx-1,ny-1,nz-1)    )
     allocate( qa(nx-1,ny-1,nz-1)    )
     allocate( zh(nx-1,ny-1,nz-1)    )
     allocate( pa(nx-1,ny-1,nz-1)    )
     allocate( pab(nx-1,ny-1,nz-1)   )
     allocate( ph(nx-1,ny-1,nz)      )
     allocate( phb(nx-1,ny-1,nz)     )
     allocate( mut(nx-1,ny-1,nz-1)   )
     allocate( qc(nx-1,ny-1,nz-1)    )
     allocate( qpr(nx-1,ny-1,nz-1)   )
     allocate( qps(nx-1,ny-1,nz-1)   )
     allocate( qpg(nx-1,ny-1,nz-1)   )
     allocate( tau(nx-1,ny-1,nz-1)   )
     allocate( tke(nx-1,ny-1,nz)     )
     allocate( el(nx-1,ny-1,nz)      )
     allocate( ql(nx-1,ny-1,nz-1)    )
     allocate( qi(nx-1,ny-1,nz-1)    )
     allocate( kl(nx-1,ny-1,nz-1)    )
     allocate( ki(nx-1,ny-1,nz-1)    )
     allocate( kf(nx-1,ny-1,nz-1)    )
     allocate( kw(nx-1,ny-1,nz-1)    )
     allocate( kp(nx-1,ny-1,nz-1)    )
     allocate( kd(nx-1,ny-1,nz-1)    )
     allocate( ke(nx-1,ny-1,nz-1)    )
     allocate( kdf(nx-1,ny-1,nz-1)   )
     allocate( ksf(nx-1,ny-1,nz-1)   )
     allocate( kud(nx-1,ny-1,nz-1)   )
     allocate( kue(nx-1,ny-1,nz-1)   )
     allocate( kde(nx-1,ny-1,nz-1)   )
     allocate( kdd(nx-1,ny-1,nz-1)   )
     allocate( taold(nx-1,ny-1,nz-1) )
     allocate( paold(nx-1,ny-1,nz-1) )
     allocate( qvold(nx-1,ny-1,nz-1) )
     allocate( qcold(nx-1,ny-1,nz-1) )

     allocate( clu(nx-1,ny-1,26)  )
     allocate( rluf(nx-1,ny-1,50) )

     allocate( smois(nx-1,ny-1,4) )
     allocate( stemp(nx-1,ny-1,4) )
 
     allocate( xdprj(nx-1,ny)    )
     allocate( ydprj(nx-1,ny)    )
     allocate( mu(nx-1,ny-1)     )
     allocate( mub(nx-1,ny-1)    )
     allocate( rainc(nx-1,ny-1)  )
     allocate( rainr(nx-1,ny-1)  )
     allocate( tsrf(nx-1,ny-1)   )
     allocate( pbl(nx-1,ny-1)    )
     allocate( rlu(nx-1,ny-1)    )
     allocate( z0(nx-1,ny-1)     )
     allocate( ylat(nx-1,ny-1)   )
     allocate( xlon(nx-1,ny-1)   )
     allocate( snocvr(nx-1,ny-1) )
     allocate( topo(nx-1,ny-1)   )
     allocate( rainro(nx-1,ny-1) )
     allocate( rainco(nx-1,ny-1) )
     allocate( cfr(nx-1,ny-1)    )
     allocate( timec(nx-1,ny-1)  )
     allocate( tc(nx-1,ny-1)     )
     allocate( smin(nx-1,ny-1)   )
     allocate( stin(nx-1,ny-1)   )
     allocate( stypin(nx-1,ny-1) )
     allocate( ctopw(nx-1,ny-1)  )
     allocate( capew(nx-1,ny-1)  )
     allocate( cfin(nx-1,ny-1)   )
     allocate( si(nx-1,ny-1)     )
     allocate( istin(nx-1,ny-1)  )


     allocate( eta(0:nz-1) )
     allocate( c1h(nz-1) )
     allocate( c2h(nz-1) )

     allocate( ydot(ny)    )
     allocate( xcrs(nx-1)  )
     allocate( ycrs(ny-1)  )
     allocate( zz(nz-1)    )
     allocate( dz(nz-1)    )
     allocate( tt(nz-1)    )
     allocate( pr(nz-1)    )
     allocate( qq(nz-1)    )
     allocate( qv(nz-1)    )
     allocate( rh(nz-1)    )
     allocate( cfrac(nz-1) )
     allocate( qr(nz-1)    )
     allocate( qs(nz-1)    )
     allocate( qg(nz-1)    )

     allocate( u10in(nx-1,ny-1) )
     allocate( v10in(nx-1,ny-1) )
     allocate( t2in(nx-1,ny-1)  )
     allocate( swd(nx-1,ny-1)   )

  end subroutine
!--------------------------------------------------------------
END MODULE
