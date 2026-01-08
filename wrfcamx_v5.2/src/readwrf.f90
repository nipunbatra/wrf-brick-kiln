  subroutine readwrf(lfirst,lnew,cdfid,itimes,project,olatin,olonin,tlat1in, &
                     tlat2in,x0camx,y0camx,dxcamx,dycamx,deltax,kvmeth, &
                     kcmaq,kmyj,kysu,scmeth,lstrat,ioffset,joffset,izone, &
                     dtout,lucat,ldiag,lwindow,lice,debug)
!
!-----READWRF reads hourly WRF data, and converts gridded data into
!     conventions for CAMx
!
  USE wrf_netcdf
  USE wrf_fields
  USE camx_fields
  USE diag_buoyancy

  implicit none
  include 'netcdf.inc'

  integer i,j,k,ic,jc,itimes,cdfid,id_var,izone,dtout
  integer kcmaq,kmyj,kysu
  integer rcode
 
  character*20 wrfproj
  character*10 project,kvmeth,scmeth
  character*40 lucatin
  character*6  lucat

  logical lfirst,lnew,lconv,debug,pcflag,prflag,lstrat,ldiag,lwindow,lice
!
!-----For WRF & CAMx grid
!
  integer ipbl                                 ! PBL Scheme
  integer icu                                  ! Sub-grid cumulus scheme
  integer isrf
  save ipbl,icu
  integer iproj,nest,parent_id
  integer n,nxin,nyin,nzin,ioffset,joffset,lu
  integer ktop,kbot
  integer ostat
  real olatin,olonin,tlat1in,tlat2in,x0camx,y0camx
  real deltax,ylato,xlono,tlat1,tlat2
  real x0wrf,y0wrf,xfwrf,yfwrf,dxcamx,dycamx
  real dylato,dxlono,dtlat1a,dtlat1b,dtlat2a,dtlat2b,dxdif
  real xsw,ysw,xne,yne,xse,yse,xnw,ynw
  real xmin,xmax,ymin,ymax
  real clat,clon
  real xoffset,yoffset,rlon,rlat
  real rtmpr,rtmpc,tv,rho,precip
  real tauqq,taukw
  real zsum
  real dum1,dum2,dum3,dum4
  real esat,qsat
  real tausum,fcld,ctrns,solmax
!
  REAL    , PARAMETER :: g            = 9.81
  REAL    , PARAMETER :: rd           = 287.
  REAL    , PARAMETER :: cp           = 7.*rd/2.
  REAL    , PARAMETER :: rcp          = rd/cp
  REAL    , PARAMETER :: p1000mb      = 100000.
  REAL    , PARAMETER :: lv           = 2.5e6    !J/kg
  REAL    , PARAMETER :: rv           = 461.     !J/K/kg
  REAL    , PARAMETER :: mvoma        = 18./28.8 !unitless
  REAL    , PARAMETER :: evap0        = 611.     !Pa
  REAL    , PARAMETER :: pi           = 3.14159265358979
  REAL    , PARAMETER :: degrad       = pi/180.
!
!-----For landuse
!
  integer ncats
  integer, parameter :: lucats_usgs   = 28
  integer, parameter :: lucats_nlcd50 = 50
  integer, parameter :: lucats_nlcd40 = 40
  integer, parameter :: lucats_igbp   = 21
  integer mlu_usgs(lucats_usgs)
  integer mlu_nlcd50(lucats_nlcd50)
  integer mlu_nlcd40(lucats_nlcd40)
  integer mlu_igbp(lucats_igbp)
  real z0_CAMx(26)

  data z0_CAMx /0.001,0.01 ,0.001,0.9  ,2.0  ,0.9  ,1.0  ,2.5  ,0.6  ,0.2  , &
                0.2  ,0.2  ,0.04 ,0.1  ,0.1  ,0.1  ,0.1  ,0.1  ,0.2  ,0.05 , &
                1.0  ,0.03 ,0.1  ,0.04 ,0.9  ,0.9  /
!
!-----Data statements to map USGS 28 category landuse to CAMx 26 category
!     landuse (supports Z03 dry deposition)
!
  data mlu_usgs(1)  /21/ !Urban
  data mlu_usgs(2)  /15/ !Dry Crop/Pasture
  data mlu_usgs(3)  /20/ !Irrigated Crop/Pasture
  data mlu_usgs(4)  /15/ !Mixed Dry/Irrigated Crop/Pasture
  data mlu_usgs(5)  /15/ !Crop/Grassland Mosaic
  data mlu_usgs(6)  /25/ !Crop/Woodland Mosaic
  data mlu_usgs(7)  /14/ !Grassland
  data mlu_usgs(8)  /11/ !Shrubland
  data mlu_usgs(9)  /11/ !Mix Shrub/Grass
  data mlu_usgs(10) /14/ !Savanna
  data mlu_usgs(11) / 7/ !Deciduous Broadleaf Forest
  data mlu_usgs(12) / 6/ !Deciduous Needleleaf Forest
  data mlu_usgs(13) / 5/ !Evergreen Broadleaf Forest
  data mlu_usgs(14) / 4/ !Evergreen Needleleaf Forest
  data mlu_usgs(15) /25/ !Mixed Forest
  data mlu_usgs(16) / 1/ !Water
  data mlu_usgs(17) /23/ !Herbaceous Wetland
  data mlu_usgs(18) /25/ !Wooded Wetland
  data mlu_usgs(19) /24/ !Barren, Sparse Vegetation
  data mlu_usgs(20) /22/ !Herbaceous Tundra
  data mlu_usgs(21) /25/ !Wooded Tundra
  data mlu_usgs(22) /22/ !Mixed Tundra
  data mlu_usgs(23) /22/ !Bare Ground Tundra
  data mlu_usgs(24) / 2/ !Snow or Ice
  data mlu_usgs(25) /24/ !Playa
  data mlu_usgs(26) /24/ !Lava
  data mlu_usgs(27) /24/ !White Sand
  data mlu_usgs(28) / 3/ !Lakes
!
!-----Data statements to map NLCD 50 category landuse to CAMx 26 category
!     landuse (supports Z03 dry deposition)
!
  data mlu_nlcd50(1)  / 1/ !NLCD-MODIS: Open Water
  data mlu_nlcd50(2)  / 2/ !NLCD-MODIS: Perennial Ice-Snow
  data mlu_nlcd50(3)  /11/ !NLCD-MODIS: Developed Open Space
  data mlu_nlcd50(4)  /21/ !NLCD-MODIS: Developed Low Intensity
  data mlu_nlcd50(5)  /21/ !NLCD-MODIS: Developed Medium Intensity
  data mlu_nlcd50(6)  /21/ !NLCD-MODIS: Developed High Intensity
  data mlu_nlcd50(7)  /24/ !NLCD-MODIS: Barren Land (Rock-Sand-Clay)
  data mlu_nlcd50(8)  /24/ !NLCD-MODIS: Unconsolidated Shore
  data mlu_nlcd50(9)  / 7/ !NLCD-MODIS: Deciduous Forest
  data mlu_nlcd50(10) / 4/ !NLCD-MODIS: Evergreen Forest
  data mlu_nlcd50(11) /25/ !NLCD-MODIS: Mixed Forest
  data mlu_nlcd50(12) /13/ !NLCD-MODIS: Dwarf Scrub
  data mlu_nlcd50(13) /11/ !NLCD-MODIS: Shrub-Scrub
  data mlu_nlcd50(14) /14/ !NLCD-MODIS: Grassland-Herbaceous
  data mlu_nlcd50(15) /11/ !NLCD-MODIS: Sedge-Herbaceous
  data mlu_nlcd50(16) /22/ !NLCD-MODIS: Lichens
  data mlu_nlcd50(17) /22/ !NLCD-MODIS: Moss
  data mlu_nlcd50(18) /22/ !NLCD-MODIS: Tundra
  data mlu_nlcd50(19) /14/ !NLCD-MODIS: Pasture-Hay
  data mlu_nlcd50(20) /15/ !NLCD-MODIS: Cultivated Crops
  data mlu_nlcd50(21) /25/ !NLCD-MODIS: Woody Wetlands
  data mlu_nlcd50(22) /25/ !NLCD-MODIS: Palustrine Forested Wetland
  data mlu_nlcd50(23) /23/ !NLCD-MODIS: Palustrine Scrub-Shrub Wetland
  data mlu_nlcd50(24) /25/ !NLCD-MODIS: Estuarine Forested Wetland
  data mlu_nlcd50(25) /23/ !NLCD-MODIS: Estuarine Scrub-Shrub Wetland
  data mlu_nlcd50(26) /23/ !NLCD-MODIS: Emergent Herbaceous Wetlands
  data mlu_nlcd50(27) /23/ !NLCD-MODIS: Palustrine Emergent Wetland
  data mlu_nlcd50(28) /23/ !NLCD-MODIS: Estuarine Emergent Wetland
  data mlu_nlcd50(29) /23/ !NLCD-MODIS: Palustrine Aquatic Bed
  data mlu_nlcd50(30) /23/ !NLCD-MODIS: Estuarine Aquatic Bed
  data mlu_nlcd50(31) / 1/ !NLCD-MODIS: Water
  data mlu_nlcd50(32) / 4/ !NLCD-MODIS: Evergreen Needleleaf Forest
  data mlu_nlcd50(33) / 5/ !NLCD-MODIS: Evergreen Broadleaf Forest
  data mlu_nlcd50(34) / 6/ !NLCD-MODIS: Deciduous Needleleaf Forest
  data mlu_nlcd50(35) / 7/ !NLCD-MODIS: Deciduous Broadleaf Forest
  data mlu_nlcd50(36) /25/ !NLCD-MODIS: Mixed Forests
  data mlu_nlcd50(37) /11/ !NLCD-MODIS: Closed Shrublands
  data mlu_nlcd50(38) /11/ !NLCD-MODIS: Open Shrublands
  data mlu_nlcd50(39) /25/ !NLCD-MODIS: Woody Savannas
  data mlu_nlcd50(40) /14/ !NLCD-MODIS: Savannas
  data mlu_nlcd50(41) /14/ !NLCD-MODIS: Grasslands
  data mlu_nlcd50(42) /23/ !NLCD-MODIS: Permanent Wetlands
  data mlu_nlcd50(43) /15/ !NLCD-MODIS: Croplands
  data mlu_nlcd50(44) /21/ !NLCD-MODIS: Urban and Built Up
  data mlu_nlcd50(45) /15/ !NLCD-MODIS: Cropland-Natural Vegetation Mosaic
  data mlu_nlcd50(46) / 2/ !NLCD-MODIS: Permanent Snow and Ice
  data mlu_nlcd50(47) /24/ !NLCD-MODIS: Barren or Sparsely Vegetated
  data mlu_nlcd50(48) / 1/ !NLCD-MODIS: IGBP Water
  data mlu_nlcd50(49) /11/ !NLCD-MODIS: unclassified
  data mlu_nlcd50(50) /11/ !NLCD-MODIS: fill value
!
!-----Data statements to map NLCD 40 category landuse to CAMx 26 category
!     landuse (supports Z03 dry deposition)
!
  data mlu_nlcd40(1)  / 4/ !NLCD-MODIS: Evergreen Needleleaf Forest
  data mlu_nlcd40(2)  / 5/ !NLCD-MODIS: Evergreen Broadleaf Forest
  data mlu_nlcd40(3)  / 6/ !NLCD-MODIS: Deciduous Needleleaf Forest
  data mlu_nlcd40(4)  / 7/ !NLCD-MODIS: Deciduous Broadleaf Forest
  data mlu_nlcd40(5)  /25/ !NLCD-MODIS: Mixed Forest
  data mlu_nlcd40(6)  /11/ !NLCD-MODIS: Closed Shrublands
  data mlu_nlcd40(7)  /11/ !NLCD-MODIS: Open Shrublands
  data mlu_nlcd40(8)  /25/ !NLCD-MODIS: Woody Savanna
  data mlu_nlcd40(9)  /14/ !NLCD-MODIS: Savanna
  data mlu_nlcd40(10) /14/ !NLCD-MODIS: Grasslands
  data mlu_nlcd40(11) /23/ !NLCD-MODIS: Permanent Wetlands
  data mlu_nlcd40(12) /15/ !NLCD-MODIS: Croplands
  data mlu_nlcd40(13) /21/ !NLCD-MODIS: Urban Built Up
  data mlu_nlcd40(14) /15/ !NLCD-MODIS: Cropland-Natural Vegetation Mosaic
  data mlu_nlcd40(15) / 2/ !NLCD-MODIS: Snow and Ice
  data mlu_nlcd40(16) /24/ !NLCD-MODIS: Barren or Sparsely Vegetated
  data mlu_nlcd40(17) / 1/ !NLCD-MODIS: Water
  data mlu_nlcd40(18) /11/ !NLCD-MODIS: Reserved
  data mlu_nlcd40(19) /11/ !NLCD-MODIS: Reserved
  data mlu_nlcd40(20) /11/ !NLCD-MODIS: Reserved
  data mlu_nlcd40(21) / 1/ !NLCD-MODIS: Open Water
  data mlu_nlcd40(22) / 2/ !NLCD-MODIS: Perennial Ice and snow
  data mlu_nlcd40(23) /11/ !NLCD-MODIS: Developed Open Space
  data mlu_nlcd40(24) /21/ !NLCD-MODIS: Developed Low Intensity
  data mlu_nlcd40(25) /21/ !NLCD-MODIS: Developed Medium Intensity
  data mlu_nlcd40(26) /21/ !NLCD-MODIS: Developed High Intensity
  data mlu_nlcd40(27) /24/ !NLCD-MODIS: Barren Rock
  data mlu_nlcd40(28) / 7/ !NLCD-MODIS: Deciduous Forest
  data mlu_nlcd40(29) / 4/ !NLCD-MODIS: Evergreen Forest
  data mlu_nlcd40(30) /25/ !NLCD-MODIS: Mixed Forest
  data mlu_nlcd40(31) /13/ !NLCD-MODIS: Dwarf Scrub
  data mlu_nlcd40(32) /11/ !NLCD-MODIS: Shrub-Scrub
  data mlu_nlcd40(33) /14/ !NLCD-MODIS: Grassland-Herbaceous
  data mlu_nlcd40(34) /11/ !NLCD-MODIS: Sedge-Herbaceous
  data mlu_nlcd40(35) /22/ !NLCD-MODIS: Lichens
  data mlu_nlcd40(36) /22/ !NLCD-MODIS: Moss
  data mlu_nlcd40(37) /14/ !NLCD-MODIS: Pasture-Hay
  data mlu_nlcd40(38) /15/ !NLCD-MODIS: Cultivated Crops
  data mlu_nlcd40(39) /25/ !NLCD-MODIS: Woody Wetlands
  data mlu_nlcd40(40) /23/ !NLCD-MODIS: Emergent Herbaceous Wetlands
!
!-----Data statements to map IGBP MODIS 21 category landuse to CAMx 26 category
!     landuse (supports Z03 dry deposition)
!
  data mlu_igbp(1)   / 4/  !Evergreen Needleleaf Forest
  data mlu_igbp(2)   / 5/  !Evergreen Broadleaf Forest
  data mlu_igbp(3)   / 6/  !Deciduous Needleleaf Forest
  data mlu_igbp(4)   / 7/  !Deciduous Broadleaf Forest
  data mlu_igbp(5)   /25/  !Mixed Forests
  data mlu_igbp(6)   /11/  !Closed Shrublands
  data mlu_igbp(7)   /11/  !Open Shrublands
  data mlu_igbp(8)   /25/  !Woody Savannas
  data mlu_igbp(9)   /14/  !Savannas
  data mlu_igbp(10)  /14/  !Grasslands
  data mlu_igbp(11)  /23/  !Permanent Wetlands
  data mlu_igbp(12)  /15/  !Croplands
  data mlu_igbp(13)  /21/  !Urban and Built-Up
  data mlu_igbp(14)  /15/  !Cropland/Natural Vegetation Mosaic
  data mlu_igbp(15)  / 2/  !Snow and Ice
  data mlu_igbp(16)  /24/  !Barren or Sparsely Vegetated
  data mlu_igbp(17)  / 1/  !Water (assumed to be oceans)
  data mlu_igbp(18)  /25/  !Wooded Tundra
  data mlu_igbp(19)  /22/  !Mixed Tundra
  data mlu_igbp(20)  /22/  !Barren Tundra
  data mlu_igbp(21)  / 3/  !Lakes (fresh water)
!  
!-----First time through, allocate WRF and CAMx fields, initialize cloud
!     variables, calculate grid parameters, read terrain data
!
  if (lfirst) then
!
!----Check if WRF was run with the hybrid vertical coordinate
!    (hybrid terrain-following + isobaric aloft)
!
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'HYBRID_OPT',ihyb)
     if (rcode.ne.0) then
       write(*,*)'Cannot find attribute HYBRID_OPT in WRF file'
       write(*,*)'Assuming WRF was run with original', &
                 ' terrain-following eta coordinate'
       ihyb = 0
     elseif (ihyb.ne.0) then
       write(*,*)'WRF was run with the hybrid vertical coordinate'
     else
       write(*,*)'WRF was run with the original eta vertical coordinate'
     endif
!
!-----Get WRF projection and grid info at the first time of the 1st WRF file
!
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'MAP_PROJ',iproj)
     if (iproj.eq.1) then
       wrfproj = 'Lambert Conformal'
       if (project.ne.'LAMBERT' .and. project.ne.'UTM' .and. &
           project.ne.'LATLON') then
         write(*,*)'CAMx projection does not work with WRF projection'
         write(*,*)'WRF:  ',wrfproj
         write(*,*)'CAMx: ',project
         write(*,*)'Program stopping'
         stop
       endif
     elseif (iproj.eq.2) then
       wrfproj = 'Polar Stereographic'
       if (project.ne.'POLAR') then
         write(*,*)'CAMx projection does not work with WRF projection'
         write(*,*)'WRF:  ',wrfproj
         write(*,*)'CAMx: ',project
         write(*,*)'Program stopping'
         stop
       endif
     elseif (iproj.eq.3) then
       wrfproj = 'Mercator'
       if (project.ne.'MERCATOR' .and. project.ne.'UTM' .and. &
           project.ne.'LATLON') then
         write(*,*)'CAMx projection does not work with WRF projection'
         write(*,*)'WRF:  ',wrfproj
         write(*,*)'CAMx: ',project
         write(*,*)'Program stopping'
         stop
       endif
     else
       write(*,*)'Unrecognized WRF map projection'
       write(*,*)'WRFCAMx can process these WRF projections:'
       write(*,*)'1 = Lambert Conformal'
       write(*,*)'2 = Polar Stereographic'
       write(*,*)'3 = Mercator'
       stop
     endif
!
!-----Allocate arrays
!
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'WEST-EAST_GRID_DIMENSION',nxin)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'SOUTH-NORTH_GRID_DIMENSION',nyin)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'BOTTOM-TOP_GRID_DIMENSION',nzin)
     call alloc_wrf(nxin,nyin,nzin)
     nx = nxin - 1
     ny = nyin - 1
     nz = nzin - 1
     ql = 0.
     qi = 0.
     qc = 0.
     qcold = 0.
     kl = 0.
     ki = 0.
     kw = 0.
     kf = 0.
     kdf = 0.
     ksf = 0.
     kud = 0.
     kue = 0.
     kdd = 0.
     kde = 0.
     ke  = 0.
     kd  = 0.
     timec = 0.
     qpr = 0.
     qps = 0.
     qpg = 0.
     if (nzc.gt.nz) then
        write(*,*)'CAMx vertical dimension exceeds WRF dimension:'
        write(*,*)'WRF:  ',nz
        write(*,*)'CAMx: ',nzc
        write(*,*)'Program stopping'
        stop
     endif
     call alloc_camx(nxc,nyc,nzc,nz)
!
!----Get WRF map projection info
!
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'DX',deltax)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'MOAD_CEN_LAT',ylato)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'STAND_LON',xlono)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'TRUELAT1',tlat1)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'TRUELAT2',tlat2)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'CEN_LAT',clat)
     rcode = NF_GET_ATT_REAL(cdfid,nf_global,'CEN_LON',clon)
     call get_var_2d_real_cdf(cdfid,'XLAT',ylat,nx,ny,itimes,debug)
     call get_var_2d_real_cdf(cdfid,'XLONG',xlon,nx,ny,itimes,debug)
 
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'GRID_ID',nest)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'PARENT_ID',parent_id)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'BL_PBL_PHYSICS',ipbl)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'CU_PHYSICS',icu)
     rcode = NF_GET_ATT_INT(cdfid,nf_global,'SF_SURFACE_PHYSICS',isrf)
     nsrf = 0
     if (isrf.eq.2) then
       nsrf = 4
     elseif (isrf.eq.7) then
       nsrf = 2
     endif
     if (nsrf.eq.0 .and. ldiag) then
       write(*,*) 'Unrecognized soil model:'
       write(*,*) 'soil moisture and temperature not processed'
     endif
     if ((kvmeth.eq.'MYJ' .or. kvmeth.eq.'ALL') .and. &
         (ipbl.ne.2 .and. ipbl.ne.4)) then
       write(*,*) 'You chose the MYJ TKE method for Kv diagnosis.' 
       write(*,*) 'WRF was not run with a TKE scheme.'
       if (kvmeth.eq.'ALL') then
         write(*,*) 'MYJ Kv fields cannot be calculated.'
         kmyj = 0
       else
         write(*,*) 'Choose another Kv option.'
         stop
       endif
     endif
     if (scmeth.eq.'KF' .and. icu.ne.1 .and. icu.ne.11) then
       write(*,*) 'You chose the KF method for sub-grid clouds.'
       write(*,*) 'WRF was not run with the KF cumulus scheme.'
       write(*,*) 'Turn off the KF method.'
       stop
     endif

     dylato = abs(ylato - olatin)
     dxlono = abs(xlono - olonin)
     dtlat1a = abs(tlat1 - tlat1in)
     dtlat1b = abs(tlat1 - tlat2in)
     dtlat2a = abs(tlat2 - tlat2in)
     dtlat2b = abs(tlat2 - tlat1in)
     if (iproj.eq.1 .and. dxlono.le.0.001 .and.        &
         ((dtlat1a.le.0.001 .and. dtlat2a.le.0.001) .or.   &
          (dtlat1b.le.0.001 .and. dtlat2b.le.0.001))) ylato = olatin

     if (iproj.eq.1) then
      call lcpgeo(0,xlono,ylato,tlat1,tlat2,xmin,ymin,xlon(1,1),ylat(1,1))
      call lcpgeo(0,xlono,ylato,tlat1,tlat2,xmax,ymax,xlon(nx,ny),ylat(nx,ny))
     elseif (iproj.eq.2) then
       call pspgeo(0,xlono,ylato,tlat1,xmin,ymin,xlon(1,1),ylat(1,1))
       call pspgeo(0,xlono,ylato,tlat1,xmax,ymax,xlon(nx,ny),ylat(nx,ny))
     elseif (iproj.eq.3) then
       call mrcgeo(0,xlono,ylato,tlat1,xmin,ymin,xlon(1,1),ylat(1,1))
       call mrcgeo(0,xlono,ylato,tlat1,xmax,ymax,xlon(nx,ny),ylat(nx,ny))
     endif

     deltax = deltax/1000.
     dxdif = abs(deltax - dxcamx)
     x0wrf = xmin - deltax/2.
     y0wrf = ymin - deltax/2.
     xfwrf = xmax + deltax/2.
     yfwrf = ymax + deltax/2.
    
     write(*,*)
     write(*,*)'Grid parameters for the WRF domain'
     write(*,'(a,i2,1x,a20)')'        Projection:',iproj,wrfproj
     write(*,'(a,2f10.3)')   '    Origin Lat/lon:',ylato,xlono
     if (iproj.eq.1) then
       write(*,'(a,2f10.3)') '    True Latitudes:',tlat1,tlat2
     else
       write(*,'(a,2f10.3)') '     True Latitude:',tlat1
     endif
     write(*,*)
     write(*,'(a,2i10)')     '           GRID ID:',nest
     write(*,'(a,2i10)')     '         PARENT ID:',parent_id
     write(*,'(a,3i10)')     '          NX,NY,NZ:',nx,ny,nz
     write(*,'(a,f10.3)')    '           DX (km):',deltax
     write(*,'(a,2f10.3)')   ' Grid center (deg):',clat,clon
     write(*,'(a,2f10.3)')   'SW x/y corner (km):',x0wrf,y0wrf
     write(*,'(a,2f10.3)')   'NE x/y corner (km):',xfwrf,yfwrf
     write(*,*)
!
!-----Compare WRF projection/resolution to CAMx projection/resolution
!
!-----Matching projections and grids
!
     lwindow = .false.
     if (project.eq.'LAMBERT   ') then
       if (iproj.eq.1 .and. dxlono.le.0.001 .and.        &
         ((dtlat1a.le.0.001 .and. dtlat2a.le.0.001) .or.   &
          (dtlat1b.le.0.001 .and. dtlat2b.le.0.001)) .and. &
           dxdif.le.0.001) then
         lwindow = .true.
         write(*,'(a)')'CAMx Lambert/grid identical to WRF Lambert/grid'
       endif
     elseif (project.eq.'POLAR     ') then
       lwindow = .true.
       if (iproj.ne.2 .or. dxlono.gt.0.001 .or. dtlat1a.gt.0.001 .or. &
           dxdif.gt.0.001) then
         write(*,'(a)')'CAMx Polar/grid NOT identical to WRF Polar/grid'
         write(*,*)'Program stopping'
         stop
       endif
     elseif (project.eq.'MERCATOR  ') then
       if (iproj.eq.3 .and. dxlono.le.0.001 .and. dtlat1a.le.0.001 .and. &
           dxdif.le.0.001) then
         lwindow = .true.
         write(*,'(a)')'CAMx Mercator/grid identical to WRF Mercator/grid'
       endif
     endif

     if (lwindow) then
       xoffset = (x0camx - x0wrf)/deltax
       yoffset = (y0camx - y0wrf)/deltax
       ioffset = nint(xoffset)
       joffset = nint(yoffset)
       if (abs(xoffset-float(ioffset)).gt.0.01 .or. &
           abs(yoffset-float(joffset)).gt.0.01) then
         write(*,'(a)') 'CAMx SW grid corner does not align with WRF grid'
         write(*,*)'Program stopping'
         write(*,*)'x0,y0 CAMx: ',x0camx,y0camx
         write(*,*)'x0,y0 WRF : ',x0wrf,y0wrf
         write(*,*)'x/y offset: ',xoffset,yoffset
         stop
       endif
       write(*,'(a,2i10)')'      I,J offsets:',ioffset,joffset
       if (ioffset.lt.0 .or. joffset.lt.0) then
         write(*,*)'CAMx grid dimensions extend beyond WRF grid:'
         write(*,*)'I,J offsets cannot be < 0'
         write(*,*)'Program stopping'
         stop
       endif
       if (ioffset+nxc.gt.nx .or. joffset+nyc.gt.ny) then
         write(*,*)'CAMx grid dimensions extend beyond WRF grid:'
         write(*,*)'WRF dimensions:              ',nx,ny
         write(*,*)'CAMx dimensions with offset: ',ioffset+nxc,joffset+nyc
         write(*,*)'Program stopping'
         stop
       endif
       xsw = x0wrf + dxcamx*ioffset
       ysw = y0wrf + dxcamx*joffset
       xne = xsw + dxcamx*nxc
       yne = ysw + dycamx*nyc
       xse = xne
       yse = ysw
       xnw = xsw
       ynw = yne
     endif
!
!-----Different LCP or MERCATOR or interpolating to UTM or LATLON
!
     if (.not.lwindow) then
       if (project.eq.'LATLON    ') then
         write(*,'(a)')'Interpolating WRF to CAMx LATLON'
       elseif (project.eq.'UTM       ') then
         write(*,'(a)')'Interpolating WRF to CAMx UTM'
       elseif (project.eq.'LAMBERT   ') then
         write(*,'(a)')'Interpolating WRF to different CAMx LAMBERT'
       elseif (project.eq.'MERCATOR  ') then
         write(*,'(a)')'Interpolating WRF to different CAMx MERCATOR'
       endif
       do i = 1,nx
         xcrs(i) = x0wrf + deltax*(i - 0.5)
       enddo
       do j = 1,ny
         ycrs(j) = y0wrf + deltax*(j - 0.5)
       enddo
       do j = 1,ny+1
         ydot(j) = y0wrf + deltax*(j - 1.)
       enddo
       do ic = 1,nxc
         xc(ic) = x0camx + dxcamx*(ic - 0.5)
       enddo
       do jc = 1,nyc
         yc(jc) = y0camx + dycamx*(jc - 0.5)
       enddo
!
!-----Map WRF cross and Y-dot points to CAMx projection for later wind rotation
!
       do j = 1,ny+1
         do i = 1,nx
           if (iproj.eq.1) then
             call lcpgeo(1,xlono,ylato,tlat1,tlat2,xcrs(i),ydot(j),rlon,rlat)
           else
             call mrcgeo(1,xlono,ylato,tlat1,xcrs(i),ydot(j),rlon,rlat)
           endif
           if (project.eq.'LATLON    ') then
             xdprj(i,j) = rlon
             ydprj(i,j) = rlat
           elseif (project.eq.'LAMBERT   ') then
             call lcpgeo(0,olonin,olatin,tlat1in,tlat2in,xdprj(i,j), &
                         ydprj(i,j),rlon,rlat)
           elseif (project.eq.'UTM       ') then
             call utmgeo(0,izone,xdprj(i,j),ydprj(i,j),rlon,rlat)
           elseif (project.eq.'MERCATOR  ') then
             call mrcgeo(0,olonin,olatin,tlat1in,xdprj(i,j),ydprj(i,j), &
                         rlon,rlat)
           endif
         enddo
       enddo
!
!-----Map CAMx cell midpoints to WRF projection
!
       do jc = 1,nyc
         do ic = 1,nxc
           if (project.eq.'LATLON    ') then
             rlon = xc(ic)
             rlat = yc(jc)
           elseif (project.eq.'LAMBERT   ') then
             call lcpgeo(1,olonin,olatin,tlat1in,tlat2in,xc(ic),yc(jc), &
                         rlon,rlat)
           elseif (project.eq.'UTM       ') then
             call utmgeo(1,izone,xc(ic),yc(jc),rlon,rlat)
           elseif (project.eq.'MERCATOR  ') then
             call mrcgeo(1,olonin,olatin,tlat1in,xc(ic),yc(jc),rlon,rlat)
           endif
           if (iproj.eq.1) then
             call lcpgeo(0,xlono,ylato,tlat1,tlat2,xcwrf(ic,jc),ycwrf(ic,jc), &
                         rlon,rlat)
           else
             call mrcgeo(0,xlono,ylato,tlat1,xcwrf(ic,jc),ycwrf(ic,jc), &
                         rlon,rlat)
           endif
         enddo
       enddo
       xsw = xcwrf(1,1) - (xcwrf(2,1) - xcwrf(1,1))/2.
       ysw = ycwrf(1,1) - (ycwrf(1,2) - ycwrf(1,1))/2.
       xse = xcwrf(nxc,1) + (xcwrf(nxc,1) - xcwrf(nxc-1,1))/2.
       yse = ycwrf(nxc,1) - (ycwrf(nxc,2) - ycwrf(nxc,1))/2.
       xne = xcwrf(nxc,nyc) + (xcwrf(nxc,nyc) - xcwrf(nxc-1,nyc))/2.
       yne = ycwrf(nxc,nyc) + (ycwrf(nxc,nyc) - ycwrf(nxc,nyc-1))/2.
       xnw = xcwrf(1,nyc) - (xcwrf(2,nyc) - xcwrf(1,nyc))/2.
       ynw = ycwrf(1,nyc) + (ycwrf(1,nyc) - ycwrf(1,nyc-1))/2.
     endif
!
!-----Echo CAMx grid parameters
!
     write(*,*)
     write(*,*)'Grid parameters for the CAMx domain in WRF space'
     write(*,'(a,3i10)')  '         NX,NY,NZ:',nxc,nyc,nzc
     write(*,'(a,2f10.3)')'            DX,DY:',dxcamx,dycamx
     write(*,'(a,2f10.3)')'        SW corner:',xsw,ysw
     write(*,'(a,2f10.3)')'        SE corner:',xse,yse
     write(*,'(a,2f10.3)')'        NE corner:',xne,yne
     write(*,'(a,2f10.3)')'        NW corner:',xnw,ynw
     write(*,*)

     if (.not.lwindow .and. &
         (xsw.lt.xmin .or. xnw.lt.xmin .or. &
          ysw.lt.ymin .or. yse.lt.ymin .or. &
          xse.gt.xmax .or. xne.gt.xmax .or. &
          ynw.gt.ymax .or. yne.gt.ymax)) then
       write(*,*)'CAMx grid spans outside WRF data grid'
       write(*,*)'Data interpolation is not possible'
       write(*,*)'Program stopping'
       stop
     endif
!
!-----Get interpolation parameters for different LCP, UTM, MERCATOR or LATLON
!
     if (.not.lwindow) then
       do jc = 1,nyc
         do ic = 1,nxc
           do i = 1,nx
             if (xcrs(i).gt.xcwrf(ic,jc)) then
               icrs(ic,jc) = i
               goto 102
             endif
           enddo
           write(*,*)'Did not find icrs',ic,jc,xcwrf(ic,jc)
           stop
 102       do j = 1,ny
             if (ycrs(j).gt.ycwrf(ic,jc)) then
               jcrs(ic,jc) = j
               goto 104
             endif
           enddo
           write(*,*)'Did not find jcrs',ic,jc,ycwrf(ic,jc)
           stop
 104       continue
         enddo
       enddo
     endif
!
!-----Read terrestrial fields (topo, landuse)
!
     rcode = nf_inq_varid(cdfid,'HGT',id_var)
     if (rcode .eq. 0) then
       call get_var_2d_real_cdf(cdfid,'HGT',topo,nx,ny,itimes,debug)
     else
        write(*,*)'***************************************'
        write(*,*)'Cannot find HGT topography field'
        write(*,*)'Topo fields will not be processed!'
        write(*,*)'***************************************'
        topo = 0.
     endif

     lucat = '      '
     rcode = NF_GET_ATT_TEXT(cdfid,nf_global,'MMINLU',lucatin)
     if (lucatin(1:4).eq.'USGS')   lucat = 'USGS  '
     if (lucatin(1:6).eq.'MODIFI') lucat = 'IGBP  '
     if (lucatin(1:4).eq.'NLCD')   lucat = 'NLCD  '
     if (lucat.ne.'USGS  ' .and. lucat.ne.'NLCD  ' .and. lucat.ne.'IGBP  ') then
       write(*,*)'*********************************************'
       write(*,*)'USGS, NLCD, IGBP landuse fields were not used'
       write(*,*)'Landuse fields will not be processed!'
       write(*,*)'*********************************************'
       clu = 0.
     else
       write(*,*)'WRF landuse is: ',lucat
       rcode = NF_GET_ATT_INT(cdfid,nf_global,'NUM_LAND_CAT',ncats)
       write(*,*)'Number of landuse categories: ',ncats
       if (lucat.eq.'USGS  ' .and. ncats.ne.lucats_usgs   .and. &
                                   ncats.ne.lucats_usgs-1 .and. &
                                   ncats.ne.lucats_usgs-4) then
         write(*,*)'*********************************************'
         write(*,*)'Expected number of USGS categories:',lucats_usgs, &
                                                         lucats_usgs-1, &
                                                         lucats_usgs-4
         write(*,*)'STOPPING'
         write(*,*)'*********************************************'
         stop
       endif
       if (lucat.eq.'IGBP  ' .and. ncats.ne.lucats_igbp-1 .and. &
                                   ncats.ne.lucats_igbp) then
         write(*,*)'*********************************************'
         write(*,*)'Expected number of IGBP categories:',lucats_igbp-1, &
                                                         lucats_igbp
         write(*,*)'STOPPING'
         write(*,*)'*********************************************'
         stop
       endif
       if (lucat.eq.'NLCD  ' .and. ncats.ne.lucats_nlcd40 .and. &
                                   ncats.ne.lucats_nlcd50) then
         write(*,*)'*********************************************'
         write(*,*)'Expected number of NLCD categories:',lucats_nlcd40, &
                                                         lucats_nlcd50
         write(*,*)'STOPPING'
         write(*,*)'*********************************************'
         stop
       endif
       rcode = nf_inq_varid(cdfid,'LANDUSEF',id_var)
       if (rcode .eq. 0) then
         call get_var_3d_real_cdf(cdfid,'LANDUSEF',rluf,nx,ny,ncats, &
                                  itimes,debug)
         write(*,*)'Using fractional landuse fields'
         do j = 1,ny
           do i = 1,nx
             do n = 1,26
               clu(i,j,n) = 0.
             enddo
             do n = 1,ncats
               if (ncats.eq.lucats_usgs .or. ncats.eq.lucats_usgs-1 .or. &
                   ncats.eq.lucats_usgs-4) then
                 lu = mlu_usgs(n)
               elseif (ncats.eq.lucats_nlcd40) then
                 lu = mlu_nlcd40(n)
               elseif (ncats.eq.lucats_nlcd50) then
                 lu = mlu_nlcd50(n)
               elseif (ncats.eq.lucats_igbp-1 .or. ncats.eq.lucats_igbp) then
                 lu = mlu_igbp(n)
               endif
               clu(i,j,lu) = clu(i,j,lu) + rluf(i,j,n)
             enddo
           enddo
         enddo
       else
         rcode = nf_inq_varid(cdfid,'LU_INDEX',id_var)
         if (rcode .eq. 0) then
           call get_var_2d_real_cdf(cdfid,'LU_INDEX',rlu,nx,ny,itimes,debug)
           write(*,*)'Using dominant landuse fields'
           do j = 1,ny
             do i = 1,nx
               lu = int(rlu(i,j))
               do n = 1,26
                 clu(i,j,n) = 0.
               enddo
               if (ncats.eq.lucats_usgs .or. ncats.eq.lucats_usgs-1 .or. &
                   ncats.eq.lucats_usgs-4) then
                 clu(i,j,mlu_usgs(lu)) = 1.
               elseif (ncats.eq.lucats_nlcd40) then
                 clu(i,j,mlu_nlcd40(lu)) = 1.
               elseif (ncats.eq.lucats_nlcd50) then
                 clu(i,j,mlu_nlcd50(lu)) = 1.
               elseif (ncats.eq.lucats_igbp-1 .or. ncats.eq.lucats_igbp) then
                 clu(i,j,mlu_igbp(lu)) = 1.
               endif
             enddo
           enddo
         else
           write(*,*)'*********************************************'
           write(*,*)'Cannot find USGS, NLCD or IGBP landuse fields'
           write(*,*)'Landuse fields will not be processed!'
           write(*,*)'*********************************************'
           clu = 0.
         endif
       endif
     endif
!
!-----Read seaice fields and adjust CAMx LU accordingly
!
     if (lice) then
       rcode = nf_inq_varid(cdfid,'SEAICE',id_var)
       if (rcode .eq. 0) then
         call get_var_2d_real_cdf(cdfid,'SEAICE',si,nx,ny,itimes,debug)
         do j = 1,ny
           do i = 1,nx
             if (si(i,j).gt.0.) then
               clu(i,j,2) = min(clu(i,j,1),si(i,j)) + clu(i,j,2)
               clu(i,j,1) = max(0.,clu(i,j,1) - si(i,j))
             endif
           enddo
         enddo
       else
         write(*,*)'Cannot find seaice mask in the WRF file:'
         write(*,*)'Sea ice will not be applied to LUs #1 and #2'
       endif
     endif
!
!-----Set surface roughness
!
     if (lucat.eq.'USGS  ' .or. lucat.eq.'NLCD  ' .or. lucat.eq.'IGBP  ') then
       write(*,*) 'Setting surface rougness using mapped CAMx landuse'
       z0 = 0.
       do j = 1,ny
         do i = 1,nx
           do n = 1,26
             z0(i,j) = z0(i,j) + z0_CAMx(n)*clu(i,j,n)
           enddo
         enddo
       enddo
     else
       write(*,*)'Surface roughness cannot be set from CAMx landuse'
       if (kvmeth.eq.'CMAQ' .or. kvmeth.eq.'YSU' .or. kvmeth.eq.'ALL') then
         write(*,*)'CMAQ and YSU Kv fields cannot not be calculated.'
         if (kvmeth.eq.'ALL') then
           kcmaq = 0
           kysu  = 0
         else
           write(*,*) 'Choose another Kv option.'
           stop 
         endif
       endif
     endif
  endif
!
!-----Read time-variable fields
!
!-----Start with basic state fields
!
  rcode = nf_inq_varid(cdfid,'P',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Perturbation pressure is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'P',pa,nx,ny,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'PB',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Base pressure is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'PB',pab,nx,ny,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'PH',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Perturbation geopotential is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'PH',ph,nx,ny,nz+1,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'PHB',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Base geopotential is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'PHB',phb,nx,ny,nz+1,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'U',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'U-wind is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'U',ua,nx+1,ny,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'V',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'V-wind is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'V',va,nx,ny+1,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'T',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Temperature is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'T',ta,nx,ny,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'QVAPOR',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'Water vapor is not available in the WRF file'
    stop
  else
    call get_var_3d_real_cdf(cdfid,'QVAPOR',qa,nx,ny,nz,itimes,debug)
  endif
!
!-----Cloud hydrometeors
!
  rcode = nf_inq_varid(cdfid,'QCLOUD',id_var)
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'QCLOUD is not available in the WRF file'
    ql = 0.
  else
    call get_var_3d_real_cdf(cdfid,'QCLOUD',ql,nx,ny,nz,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'QICE',id_var)
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'QICE is not available in the WRF file'
    qi = 0.
  else
    call get_var_3d_real_cdf(cdfid,'QICE',qi,nx,ny,nz,itimes,debug)
  endif
!
!-----Sub-grid cloud parameters (K-F cumulus only!)
!
  if ((icu.eq.1 .or. icu.eq.11) .and. scmeth.eq.'KF') then
    rcode = nf_inq_varid(cdfid,'CLDFRA_DP',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'CLDFRA_DP is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'CLDFRA_DP',kdf,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'CLDFRA_SH',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'CLDFRA_SH is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'CLDFRA_SH',ksf,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'QC_CU',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'QC_CU is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'QC_CU',kl,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'QI_CU',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'QI_CU is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'QI_CU',ki,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'UDR_KF',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'UDR_KF is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'UDR_KF',kud,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'UER_KF',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'UER_KF is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'UER_KF',kue,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'DDR_KF',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'DDR_KF is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'DDR_KF',kdd,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'DER_KF',id_var)
    if (rcode .ne. 0) then
      write(*,*) 'DER_KF is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_3d_real_cdf(cdfid,'DER_KF',kde,nx,ny,nz,itimes,debug)
    endif
    rcode = nf_inq_varid(cdfid,'TIMEC_KF',id_var) !!!! Assume 2-D !!!!
    if (rcode .ne. 0) then
      write(*,*) 'TIMEC_KF is not available in the WRF file'
      write(*,*) 'You chose the KF method for sub-grid cumulus.'
      write(*,*) 'Turn off the KF option.'
      stop
    else
      call get_var_2d_real_cdf(cdfid,'TIMEC_KF',timec,nx,ny,itimes,debug)
    endif
  endif
!
!-----TKE for Kv diagnosis
!
  if (kmyj.ne.0) then
    rcode = nf_inq_varid(cdfid,'TKE_PBL',id_var)
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'TKE_PBL is not available in the WRF file.'
      if (lfirst) write(*,*) 'You chose the TKE method for Kv diagnosis.' 
      kmyj = 0
    endif
    rcode = nf_inq_varid(cdfid,'EL_PBL',id_var)
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'EL_PBL is not available in the WRF file.'
      if (lfirst) write(*,*) 'You chose the TKE method for Kv diagnosis.'
      kmyj = 0
    endif
    if (kmyj.ne.0) then
      call get_var_3d_real_cdf(cdfid,'TKE_PBL',tke,nx,ny,nz+1,itimes,debug)
      call get_var_3d_real_cdf(cdfid,'EL_PBL',el,nx,ny,nz+1,itimes,debug)
    else
      if (kvmeth.eq.'ALL') then
        if (lfirst) write(*,*) 'MYJ Kv fields cannot be calculated.'
      else
        write(*,*) 'Choose another Kv option.'
        stop
      endif
    endif
  endif
!
!-----Additional miscellaneous fields
!
  pcflag = .false.
  rcode = nf_inq_varid(cdfid,'PREC_ACC_C',id_var)
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'PREC_ACC_C is not available in the WRF file'
    rcode = nf_inq_varid(cdfid,'RAINC',id_var)
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'RAINC is not available in the WRF file'
      rainc = 0.
    else
      pcflag = .true.
      if (lfirst) write(*,*) 'Using RAINC for precip rate'
      call get_var_2d_real_cdf(cdfid,'RAINC',rainc,nx,ny,itimes,debug)
    endif
  else
    call get_var_2d_real_cdf(cdfid,'PREC_ACC_C',rainc,nx,ny,itimes,debug)
  endif

  prflag = .false.
  rcode = nf_inq_varid(cdfid,'PREC_ACC_NC',id_var)
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'PREC_ACC_NC is not available in the WRF file'
    rcode = nf_inq_varid(cdfid,'RAINNC',id_var)
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'RAINNC is not available in the WRF file'
      rainr = 0.
    else
      prflag = .true.
      if (lfirst) write(*,*) 'Using RAINNC for precip rate'
      call get_var_2d_real_cdf(cdfid,'RAINNC',rainr,nx,ny,itimes,debug)
    endif
  else
    call get_var_2d_real_cdf(cdfid,'PREC_ACC_NC',rainr,nx,ny,itimes,debug)
  endif

  rcode = nf_inq_varid (cdfid,'TSK',id_var )
  if (rcode .ne. 0) then
    write(*,*) 'TSK is not available in the WRF file'
    stop
  else
    call get_var_2d_real_cdf(cdfid,'TSK',tsrf,nx,ny,itimes,debug)
  endif

  rcode = nf_inq_varid (cdfid,'PBLH',id_var )
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'PBLH is not available in the WRF file'
    if (kvmeth.eq.'CMAQ' .or. &
        kvmeth.eq.'YSU' .or. kvmeth.eq.'ALL') then
      if (lfirst) &
        write(*,*)'CMAQ and YSU Kv fields cannot be calculated.'
      if (kvmeth.eq.'ALL') then
        kcmaq = 0
        kysu  = 0
      else
        write(*,*) 'Choose another Kv option.'
        stop
      endif
    endif
    if (scmeth .eq. 'DIAG') then
      write(*,*) 'Sub-grid clouds cannot be diagnosed.'
      write(*,*) 'Reset the cloud diagnosis flag to NONE'
      stop
    endif
  else
    call get_var_2d_real_cdf(cdfid,'PBLH',pbl,nx,ny,itimes,debug)
  endif

  rcode = nf_inq_varid (cdfid,'SNOW',id_var ) !Find snow water content (kg/m2)
  if (rcode .ne. 0) then
    if (lfirst) write(*,*) 'SNOW is not available in the WRF file'
    snocvr = 0.
  else
    call get_var_2d_real_cdf(cdfid,'SNOW',snocvr,nx,ny,itimes,debug)
    snocvr = snocvr/1000.                     !convert kg/m2 to m
  endif

!-----Vertical coordinate data

  rcode = nf_inq_varid (cdfid,'ZNW',id_var )
  if (rcode .ne. 0) then
    write(*,*) 'ZNW (Eta coordinate) is not available in the WRF file'
    stop
  else
    call get_var_1d_real_cdf(cdfid,'ZNW',eta,nz+1,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'MU',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'MU is not available in the WRF file'
    stop
  else
    call get_var_2d_real_cdf(cdfid,'MU',mu,nx,ny,itimes,debug)
  endif

  rcode = nf_inq_varid(cdfid,'MUB',id_var)
  if (rcode .ne. 0) then
    write(*,*) 'MUB is not available in the WRF file'
    stop
  else
    call get_var_2d_real_cdf(cdfid,'MUB',mub,nx,ny,itimes,debug)
  endif

  if (ihyb.eq.0) then
    do k = 1,nz
      mut(:,:,k) = mub(:,:) + mu(:,:)
    enddo
  else
    rcode = nf_inq_varid (cdfid,'C1H',id_var )
    if (rcode .ne. 0) then
      write(*,*) 'C1H is not available in the WRF file'
      stop
    else
      call get_var_1d_real_cdf(cdfid,'C1H',c1h,nz,itimes,debug)
    endif
    rcode = nf_inq_varid (cdfid,'C2H',id_var )
    if (rcode .ne. 0) then
      write(*,*) 'C2H is not available in the WRF file'
      stop
    else
      call get_var_1d_real_cdf(cdfid,'C2H',c2h,nz,itimes,debug)
    endif
    do k = 1,nz
      mut(:,:,k) = c1h(k)*(mub(:,:) + mu(:,:)) + c2h(k)
    enddo
  endif

!-----Diagnostic data

  if (ldiag) then
    rcode = nf_inq_varid (cdfid,'U10',id_var )
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'U10 is not available in the WRF file'
      u10in = 0.
    else
      call get_var_2d_real_cdf(cdfid,'U10',u10in,nx,ny,itimes,debug)
    endif
  
    rcode = nf_inq_varid (cdfid,'V10',id_var )
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'V10 is not available in the WRF file'
      v10in = 0.
    else
      call get_var_2d_real_cdf(cdfid,'V10',v10in,nx,ny,itimes,debug)
    endif

    rcode = nf_inq_varid (cdfid,'T2',id_var )
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'T2 is not available in the WRF file'
      t2in = 0.
    else
      call get_var_2d_real_cdf(cdfid,'T2',t2in,nx,ny,itimes,debug)
    endif

    rcode = nf_inq_varid (cdfid,'SWDOWN',id_var )
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'SWDOWN is not available in the WRF file'
      swd = 0.
    else
      call get_var_2d_real_cdf(cdfid,'SWDOWN',swd,nx,ny,itimes,debug)
    endif

    smin = 0.
    stin = 0.
    if (nsrf.gt.0) then
      rcode = nf_inq_varid (cdfid,'SMOIS',id_var )
      if (rcode .ne. 0) then
        if (lfirst) write(*,*) 'SMOIS is not available in the WRF file'
      else
        call get_var_3d_real_cdf(cdfid,'SMOIS',smois,nx,ny,nsrf,itimes,debug)
        smin(:,:) = smois(:,:,1)
      endif
      rcode = nf_inq_varid (cdfid,'TSLB',id_var )
      if (rcode .ne. 0) then
        if (lfirst) write(*,*) 'TSLB is not available in the WRF file'
      else
        call get_var_3d_real_cdf(cdfid,'TSLB',stemp,nx,ny,nsrf,itimes,debug)
        stin(:,:) = stemp(:,:,1)
      endif
    endif

    rcode = nf_inq_varid (cdfid,'ISLTYP',id_var )
    if (rcode .ne. 0) then
      if (lfirst) write(*,*) 'ISLTYP is not available in the WRF file'
      stypin = 0
    else
      call get_var_2d_int_cdf(cdfid,'ISLTYP',istin,nx,ny,itimes,debug)
      do j = 1,ny
        do i = 1,nx
          stypin(i,j) = int(istin(i,j))
        enddo
      enddo
    endif

  endif
!
!-----Calculate cartesian height from geopotential height, convert height MSL
!     to height AGL (m)
!
  ph = (ph + phb)/g
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        zh(i,j,k) = ph(i,j,k+1) - ph(i,j,1)
        if (zh(i,j,k).lt.0.) then
          write(*,*)'ZH<0: at (i,j,k):',i,j,k,zh(i,j,k)
          zh(i,j,k) = 0.
        endif
      enddo
    enddo
  enddo
!
!-----Get total pressure (Pa to mb) and convert temperature from potential
!     to actual in K
!
  pa  = pa + pab
  ta = (ta + 300.)*(pa/p1000mb)**rcp
  pa = pa*0.01     ! Pa to mb
!
!-----Determine hourly incremental total precipitation rates (mm/hr)
!
  if (prflag) then
    do j = 1,ny
      do i = 1,nx
        rtmpr = rainr(i,j)
        if (lfirst .or. lnew) rainro(i,j) = rtmpr
        rainr(i,j) = amax1(0.,(rtmpr - rainro(i,j)))/(float(dtout)/60.)
        rainro(i,j) = rtmpr
      enddo
    enddo
  endif
  if (pcflag) then
    do j = 1,ny
      do i = 1,nx
        rtmpc = rainc(i,j)
        if (lfirst .or. lnew) rainco(i,j) = rtmpc
        rainc(i,j) = amax1(0.,(rtmpc - rainco(i,j)))/(float(dtout)/60.)
        rainco(i,j) = rtmpc
      enddo
    enddo
  endif
!
!-----Convert water vapor to ppm,
!     Aggregate and convert resolved and sub-grid cloud water content to g/m3
!
  solmax = 0.
  do j = 1,ny
    do i = 1,nx
      solmax = max(swd(i,j),solmax)
      do k = 1,nz
        tv = ta(i,j,k)*(1. + 0.61*qa(i,j,k))
        rho = 1000.*100.*pa(i,j,k)/rd/tv
        qa(i,j,k) = qa(i,j,k)*1.e6*28.8/18.
        qc(i,j,k) = (ql(i,j,k) + qi(i,j,k))*rho
        if (qc(i,j,k).le.0.01) qc(i,j,k) = 0.

        if (scmeth.eq.'KF') then
          kw(i,j,k) = (kl(i,j,k) + ki(i,j,k))*rho
          kf(i,j,k) = kdf(i,j,k) + ksf(i,j,k)
          ke(i,j,k) = kue(i,j,k) - kde(i,j,k)
          kd(i,j,k) = kud(i,j,k) + kdd(i,j,k)
        endif

      enddo
    enddo
  enddo
!
!-----Cloud and precip processing
!
  do 100 j = 1,ny
    do 100 i = 1,nx
!
!-----Prepare time-averaged T, P, Q, Qc profiles for diagnosis
!     of sub-grid clouds and cloud/precip water profiles
!
      do k = 1,nz
        if (lfirst .or. lnew) then
          taold(i,j,k) = ta(i,j,k)
          paold(i,j,k) = pa(i,j,k)
          qvold(i,j,k) = qa(i,j,k)
          qcold(i,j,k) = qc(i,j,k)
        endif
        zz(k) = zh(i,j,k)
        dz(k) = zz(k)
        if (k.gt.1) dz(k) = zz(k) - zz(k-1)

        tt(k) = (ta(i,j,k) + taold(i,j,k))/2.
        taold(i,j,k) = ta(i,j,k)

        pr(k) = 100.*(pa(i,j,k) + paold(i,j,k))/2.
        paold(i,j,k) = pa(i,j,k)

        qv(k) = 1.e-6*18./28.8*(qa(i,j,k) + qvold(i,j,k))/2.
        qvold(i,j,k) = qa(i,j,k)
        esat  = evap0*exp((lv/rv)*(1./273. - 1./tt(k)))
        qsat  = mvoma/(pr(k)/esat - 1.)
        rh(k) = min(100.*qv(k)/qsat,100.)

        qq(k) = (qc(i,j,k) + qcold(i,j,k))/2.
        cfrac(k) = 1.
        if (qq(k).le.0.01) then
          qq(k) = 0.
          cfrac(k) = 0.
        endif
        qcold(i,j,k) = qc(i,j,k)
      enddo
!
!-----Diagnose convective and/or stratiform sub-grid clouds
!
      lconv = .false.
      ctopw(i,j) = 0.
      capew(i,j) = 0.
      if (rainc(i,j).gt.0.2) lconv = .true.
      if (scmeth.eq.'DIAG' .or. &
         (scmeth.eq.'KF' .and. lstrat)) then
        call clddiag(nz,lconv,scmeth,lstrat,pbl(i,j),rainc(i,j),zz,pr, &
                     tt,qv,qq,cfrac,ctopw(i,j))
      endif
!
!-----KF: determine profile-average cloud fraction and convert ent/det
!     fluxes to kg/m2/s (area is fraction of grid column)
!
      if (scmeth.eq.'KF') then
        cfr(i,j) = 0.
        tc(i,j) = 0.
        zsum = 0.
        do k = 1,nz
          if (kf(i,j,k).gt.0.) then
            zsum = zsum + dz(k)
            cfr(i,j) = cfr(i,j) + kf(i,j,k)*dz(k)
          endif
        enddo
        if (zsum.gt.0.) then
          cfr(i,j) = cfr(i,j)/zsum
          tc(i,j) = timec(i,j)
        endif
        if (cfr(i,j).le.0.1) then
          cfr(i,j) = 0.
          tc(i,j) = 0.
          do k = 1,nz
            kw(i,j,k) = 0.
            kp(i,j,k) = 0.
            kf(i,j,k) = 0.
            ke(i,j,k) = 0.
            kd(i,j,k) = 0.
          enddo
        else
          ctopw(i,j) = zsum/1000.
          do k = 1,nz
            ke(i,j,k) = ke(i,j,k)/(cfr(i,j)*(deltax*1000.)**2)
            kd(i,j,k) = kd(i,j,k)/(cfr(i,j)*(deltax*1000.)**2)
          enddo
        endif
      endif
!
!-----Find resolved convective cloud tops, and use largest among resolved and
!     DIAG/KF options
!
      if (rainr(i,j).gt.0.2) then
        kbot = 0
        ktop = 0
        do k = 1,nz
          if (kbot.eq.0 .and. qq(k).gt.0.01) kbot = k
          if (kbot.gt.0 .and. qq(k).lt.0.01) then
            ktop = k - 1
            if (zz(ktop)-zz(kbot).lt.3000.) then
              kbot = 0
              ktop = 0
              cycle
            endif
            goto 105
          endif
        enddo
        if (kbot.gt.0) ktop = nz
 105    if (kbot.gt.0 .and. ktop.gt.0 .and. &
            tt(kbot).gt.275. .and. tt(ktop).lt.255.) &
          ctopw(i,j) = max(zz(ktop)/1000.,ctopw(i,j))
      endif
!
!-----Calculate CAPE for cells with non-zero CONVTOP
!
      if (ctopw(i,j).gt.0.) then
        ostat = buoyancy(nz,tt,rh,pr,zz,1,capew(i,j),dum1,dum2,dum3,dum4,3)
        if (ostat.ne.0) then
          write(*,*)'Error in BUOYANCY calculating CAPE, cell: ',i,j
          stop
        endif
        capew(i,j) = max(capew(i,j),0.)
      endif
!
!-----Determine total (resolved+sub-grid) cloud optical depth
!
      tausum = 0.
      cfin(i,j) = 0.
      do k = nz,1,-1
        tauqq = 3.*qq(k)*dz(k)/(2.*15.)*cfrac(k)**(3./2.)
        taukw = 3.*kw(i,j,k)*dz(k)/(2.*15.)*kf(i,j,k)**(3./2.)
        tau(i,j,k) = tauqq + taukw
        tausum = tausum + tau(i,j,k)
        cfin(i,j) = max(cfin(i,j),cfrac(k),kf(i,j,k))
      enddo
!
!-----Adjust surface solar flux for diagnosed cloudiness
!
      if (ldiag .and. (scmeth.eq.'DIAG' .or. scmeth.eq.'KF')) then
        if (tausum.lt.5.) then
          ctrns = 1.
          fcld  = 0.
        else
          ctrns = (5. - exp(-tausum))/(4. + 0.42*tausum)
          fcld  = 1.
        endif
        if (swd(i,j).gt.solmax/2.) then
          ctrns = 1.6*ctrns*cos(degrad*45.)
        else
          ctrns = 1.
        endif
        swd(i,j) = swd(i,j)*(1. - fcld*(1. - ctrns))
      endif
!
!-----Calculate resolved precip profile for cells containing Pr+Pc > 0.2 mm/hr
!
      precip = rainr(i,j) + rainc(i,j)
      if (scmeth.EQ.'KF') precip = rainr(i,j)
      call pcpdiag(nz,lconv,precip,zz,tt,qq,qr,qg,qs)
!
!-----Load resolved precip water into resolved water arrays
!
      do k = 1,nz
        qc(i,j,k)  = qq(k)
        qpr(i,j,k) = qr(k)
        qpg(i,j,k) = qg(k)
        qps(i,j,k) = qs(k)
      enddo
!
!-----Calculate KF sub-grid precip profile for cells containing Pc > 0.2 mm/hr
!
      if (scmeth.EQ.'KF') then
        precip = rainc(i,j)/(cfr(i,j) + 1.e-10)
        do k = 1,nz
          qq(k) = kw(i,j,k)
        enddo
        call pcpdiag(nz,.true.,precip,zz,tt,qq,qr,qg,qs)
        do k = 1,nz
          kp(i,j,k)  = qr(k) + qg(k) + qs(k)
        enddo
      endif
!
100 continue
!
  lfirst = .false.

  end subroutine readwrf
