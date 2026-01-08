      subroutine ncf_set_varatt_3dmet(n3dmet,n2dmet,n3d,n2d,ldiag,kcmaq,
     &                                kmyj,kysu,x3dmet,x2dmet,
     &                                met_3dname,met_3dunits,
     &                                met_3dlong,met_3ddesc,met_3dcoords,
     &                                met_2dname,met_2dunits,
     &                                met_2dlong,met_2ddesc,met_2dcoords)
      implicit none
c
c-----This routine sets the output variable attributes
c
c      Argument description:
c       Inputs:
c       Outputs:
c            n3dmet        I max number of 3-D met variables
c            n2dmet        I max number of 2-D met variables
c            n3d           I number of output 3-D met variables
c            n2d           I number of output 2-D met variables
c            ldiag         L diagnostic variable output flag
c            kcmaq         I unit number for CMAQ Kv file
c            kmyj          I unit number for MYJ Kv file
c            kysu          I unit number for YSU Kv file
c            x3dmet        I vector of 3-D index x-ref
c            x2dmet        I vector of 2-D index x-ref
c            met_name      C array of variable names
c            met_units     C array of units
c            met_long      C array of "long names"
c            met_desc      C array of desciptions
c            met_coords    C array of coordinates
c
      include 'ncf_iodat.inc'
c
      integer n3dmet,n2dmet,n3d,n2d
      logical ldiag
      integer kcmaq,kmyj,kysu
      integer x3dmet(*),x2dmet(*)
      character*60 met_3dname(*)
      character*60 met_3dunits(*)
      character*60 met_3dlong(*)
      character*60 met_3ddesc(*)
      character*60 met_3dcoords(*)
      character*60 met_2dname(*)
      character*60 met_2dunits(*)
      character*60 met_2dlong(*)
      character*60 met_2ddesc(*)
      character*60 met_2dcoords(*)
c
c-----Entry point:
c
      n2d = 0
      n3d = 0

      n3d = n3d+1
      x3dmet(i3dz) = n3d
      met_3dname(n3d) = "z"
      met_3dunits(n3d) = "m AGL"
      met_3dlong(n3d) = "Layer height"
      met_3ddesc(n3d) = "Layer interface heights AGL"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dp) = n3d
      met_3dname(n3d) = "pressure"
      met_3dunits(n3d) = "mb"
      met_3dlong(n3d) = "pressure"
      met_3ddesc(n3d) = "pressure"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dt) = n3d
      met_3dname(n3d) = "temperature"
      met_3dunits(n3d) = "K"
      met_3dlong(n3d) = "temperature"
      met_3ddesc(n3d) = "temperature"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dq) = n3d
      met_3dname(n3d) = "humidity"
      met_3dunits(n3d) = "ppm"
      met_3dlong(n3d) = "humidity"
      met_3ddesc(n3d) = "humidity"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3du) = n3d
      met_3dname(n3d) = "uwind"
      met_3dunits(n3d) = "m s-1"
      met_3dlong(n3d) = "longitudinal wind speed"
      met_3ddesc(n3d) = "longitudinal wind speed"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dv) = n3d
      met_3dname(n3d) = "vwind"
      met_3dunits(n3d) = "m s-1"
      met_3dlong(n3d) = "latitudinal wind speed"
      met_3ddesc(n3d) = "latitudinal wind speed"
      met_3dcoords(n3d) = "latitude longitude"

      if (n3d.gt.n3dmet) then
        write(*,*)'Number of output 3-D Met variables > max in VARATT_3DMET'
        stop
      endif

      n2d = n2d+1
      x2dmet(i2dts) = n2d
      met_2dname(n2d) = "sfctemperature"
      met_2dunits(n2d) = "K"
      met_2dlong(n2d) = "surface temperature"
      met_2ddesc(n2d) = "surface temperature"
      met_2dcoords(n2d) = "latitude longitude"

      n2d = n2d+1
      x2dmet(i2dsd) = n2d
      met_2dname(n2d) = "snowewd"
      met_2dunits(n2d) = "m"
      met_2dlong(n2d) = "snow cover equivalent water depth"
      met_2ddesc(n2d) = "snow cover equivalent water depth"
      met_2dcoords(n2d) = "latitude longitude"

      n2d = n2d+1
      x2dmet(i2dsa) = n2d
      met_2dname(n2d) = "snowage"
      met_2dunits(n2d) = "hr"
      met_2dlong(n2d) = "snow cover age"
      met_2ddesc(n2d) = "snow cover age"
      met_2dcoords(n2d) = "latitude longitude"

      if (ldiag) then
        n2d = n2d + 1
        x2dmet(i2du10) = n2d
        met_2dname(n2d) = "u10"
        met_2dunits(n2d) = "m s-1"
        met_2dlong(n2d) = "longitudinal wind at 10 m"
        met_2ddesc(n2d) = "longitudinal wind at 10 m"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dv10) = n2d
        met_2dname(n2d) = "v10"
        met_2dunits(n2d) = "m s-1"
        met_2dlong(n2d) = "latitudinal wind at 10 m"
        met_2ddesc(n2d) = "latitudinal wind at 10 m"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dt2) = n2d
        met_2dname(n2d) = "t2"
        met_2dunits(n2d) = "K"
        met_2dlong(n2d) = "temperature at 2 m"
        met_2ddesc(n2d) = "temperature at 2 m"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dsw) = n2d
        met_2dname(n2d) = "swsfc"
        met_2dunits(n2d) = "W m-2"
        met_2dlong(n2d) = "SW surface flux"
        met_2ddesc(n2d) = "SW surface flux"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dsm) = n2d
        met_2dname(n2d) = "soilmoist"
        met_2dunits(n2d) = "m3 m-3"
        met_2dlong(n2d) = "volumetric soil moisture"
        met_2ddesc(n2d) = "volumetric soil moisture"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dst) = n2d
        met_2dname(n2d) = "soiltemp"
        met_2dunits(n2d) = "K"
        met_2dlong(n2d) = "soil temperature"
        met_2ddesc(n2d) = "soil temperature"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dty) = n2d
        met_2dname(n2d) = "soiltype"
        met_2dunits(n2d) = "unitless"
        met_2dlong(n2d) = "soil texture type"
        met_2ddesc(n2d) = "soil texture type"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dwrf) = n2d
        met_2dname(n2d) = "pblwrf"
        met_2dunits(n2d) = "m"
        met_2dlong(n2d) = "PBL depth from WRF"
        met_2ddesc(n2d) = "PBL depth from WRF"
        met_2dcoords(n2d) = "latitude longitude"

        if (kcmaq.ne.0) then
          n2d = n2d + 1
          x2dmet(i2dcmaq) = n2d
          met_2dname(n2d) = "pblcmaq"
          met_2dunits(n2d) = "m"
          met_2dlong(n2d) = "PBL depth from CMAQ Kv"
          met_2ddesc(n2d) = "PBL depth from CMAQ Kv"
          met_2dcoords(n2d) = "latitude longitude"
        endif

        if (kmyj.ne.0) then
          n2d = n2d + 1
          x2dmet(i2dmyj) = n2d
          met_2dname(n2d) = "pblmyj"
          met_2dunits(n2d) = "m"
          met_2dlong(n2d) = "PBL depth from MYJ Kv"
          met_2ddesc(n2d) = "PBL depth from MYJ Kv"
          met_2dcoords(n2d) = "latitude longitude"
        endif

        if (kysu.ne.0) then
          n2d = n2d + 1
          x2dmet(i2dysu) = n2d
          met_2dname(n2d) = "pblysu"
          met_2dunits(n2d) = "m"
          met_2dlong(n2d) = "PBL depth from YSU Kv"
          met_2ddesc(n2d) = "PBL depth from YSU Kv"
          met_2dcoords(n2d) = "latitude longitude"
        endif
      endif

      if (n2d.gt.n2dmet) then
        write(*,*)
        write(*,*)'Number of output 2-D Met variables > max in VARATT_3DMET'
        stop
      endif

      return
      end
