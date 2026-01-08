  MODULE camx_fields
!
!-----CAMx grid parameters
!
  real, allocatable, dimension(:,:)    :: xcwrf,ycwrf
  real, allocatable, dimension(:,:)    :: camxlon,camxlat
  real, allocatable, dimension(:)      :: xc,yc
  integer, allocatable, dimension(:,:) :: idot,icrs,jdot,jcrs
!
!-----Intermediate CAMx fields (horizontally interpolated)
!
  real, allocatable, dimension(:,:,:) :: utmp,vtmp,ttmp,qtmp,ptmp,tktmp,eltmp
  real, allocatable, dimension(:,:,:) :: ztmp,cwtmp,prtmp,pstmp,pgtmp,odtmp
  real, allocatable, dimension(:,:,:) :: kwtmp,kptmp,ketmp,kdtmp
  real, allocatable, dimension(:,:,:) :: mutmp
  real, allocatable, dimension(:,:,:) :: lucx
  real, allocatable, dimension(:,:)   :: rain,tsfc,pblc,z0c,snow,snowo,  &
                                         snowage,topcx
  real, allocatable, dimension(:,:)   :: pbl_cmaq,pbl_myj,pbl_ysu
  real, allocatable, dimension(:,:)   :: codsum,u10,v10,t2,swsfc
  real, allocatable, dimension(:,:)   :: kftmp,kff,tctmp,kft
  real, allocatable, dimension(:,:)   :: soilm,soilt,styp
  real, allocatable, dimension(:,:)   :: ctop,cape,cldfrc
!
!-----Final CAMx fields (vertically aggregated)
!
  real, allocatable, dimension(:,:,:) :: uac,vac,tac,qac,pac,tkc,elc
  real, allocatable, dimension(:,:,:) :: zhc,cwc,pwr,pws,pwg,cod,pwtr
  real, allocatable, dimension(:,:,:) :: rkv_cmaq,rkv_myj,rkv_ysu
  real, allocatable, dimension(:,:,:) :: kfw,kfp,kfe,kfd
  real, allocatable, dimension(:,:,:) :: mutc

  integer nxc,nyc,nzc
!
!--------------------------------------------------------------------
!
  CONTAINS

  subroutine alloc_camx(nxc,nyc,nzc,nz)
  implicit none

! (nxc) -- CAMx Max WEST-EAST
! (nyc) -- CAMx Max SOUTH-NORTH
! (nzc) -- CAMx Max BOTTOM-TOP
! (nz)  -- WRF Max BOTTOM-TOP UNSTAGGERED

  integer nxc,nyc,nzc,nz

     allocate( mutmp(nxc,nyc,nz) )
     allocate( utmp(nxc,nyc,nz)  )
     allocate( vtmp(nxc,nyc,nz)  )
     allocate( ttmp(nxc,nyc,nz)  )
     allocate( qtmp(nxc,nyc,nz)  )
     allocate( ptmp(nxc,nyc,nz)  )
     allocate( tktmp(nxc,nyc,nz) )
     allocate( eltmp(nxc,nyc,nz) )
     allocate( ztmp(nxc,nyc,nz)  )
     allocate( cwtmp(nxc,nyc,nz) )
     allocate( prtmp(nxc,nyc,nz) )
     allocate( pstmp(nxc,nyc,nz) )
     allocate( pgtmp(nxc,nyc,nz) )
     allocate( odtmp(nxc,nyc,nz) )
     allocate( kwtmp(nxc,nyc,nz) )
     allocate( kptmp(nxc,nyc,nz) )
     allocate( ketmp(nxc,nyc,nz) )
     allocate( kdtmp(nxc,nyc,nz) )
 
     allocate( mutc(nxc,nyc,nzc)     )
     allocate( uac(nxc,nyc,nzc)      )
     allocate( vac(nxc,nyc,nzc)      )
     allocate( tac(nxc,nyc,nzc)      )
     allocate( qac(nxc,nyc,nzc)      )
     allocate( pac(nxc,nyc,nzc)      )
     allocate( zhc(nxc,nyc,nzc)      )
     allocate( rkv_cmaq(nxc,nyc,nzc) )
     allocate( rkv_myj(nxc,nyc,nzc)  )
     allocate( rkv_ysu(nxc,nyc,nzc)  )
     allocate( tkc(nxc,nyc,nzc)      )
     allocate( elc(nxc,nyc,nzc)      )
     allocate( cwc(nxc,nyc,nzc)      )
     allocate( pwr(nxc,nyc,nzc)      )
     allocate( pws(nxc,nyc,nzc)      )
     allocate( pwg(nxc,nyc,nzc)      )
     allocate( cod(nxc,nyc,nzc)      )
     allocate( pwtr(nxc,nyc,nzc)     )
     allocate( kfw(nxc,nyc,nzc)      )
     allocate( kfp(nxc,nyc,nzc)      )
     allocate( kfe(nxc,nyc,nzc)      )
     allocate( kfd(nxc,nyc,nzc)      )

     allocate( lucx(nxc,nyc,26) )

     allocate( xcwrf(nxc,nyc)   )
     allocate( ycwrf(nxc,nyc)   )
     allocate( idot(nxc,nyc)    )
     allocate( icrs(nxc,nyc)    )
     allocate( jdot(nxc,nyc)    )
     allocate( jcrs(nxc,nyc)    )
     allocate( rain(nxc,nyc)    )
     allocate( tsfc(nxc,nyc)    )
     allocate( pblc(nxc,nyc)    )
     allocate( z0c(nxc,nyc)     )
     allocate( snow(nxc,nyc)    )
     allocate( snowo(nxc,nyc)   )
     allocate( topcx(nxc,nyc)   )
     allocate( kftmp(nxc,nyc)   )
     allocate( kff(nxc,nyc)     )
     allocate( tctmp(nxc,nyc)   )
     allocate( kft(nxc,nyc)     )

     allocate( xc(nxc) )
     allocate( yc(nyc) )

     allocate( pbl_cmaq(nxc,nyc) )
     allocate( pbl_myj(nxc,nyc)  )
     allocate( pbl_ysu(nxc,nyc)  )
     allocate( codsum(nxc,nyc)   )
     allocate( u10(nxc,nyc)      )
     allocate( v10(nxc,nyc)      )
     allocate( t2(nxc,nyc)       )
     allocate( swsfc(nxc,nyc)    )
     allocate( soilm(nxc,nyc)    )
     allocate( soilt(nxc,nyc)    )
     allocate( styp(nxc,nyc)     )
     allocate( ctop(nxc,nyc)     )
     allocate( cape(nxc,nyc)     )
     allocate( cldfrc(nxc,nyc)   )

  end subroutine
!--------------------------------------------------------------
END MODULE
