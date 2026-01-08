      subroutine ncf_set_varatt_3dcld(n3dmet,n2dmet,n3d,n2d,ldiag,scmeth,
     &                                x3dmet,x2dmet,
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
c            scmeth        C sub-grid cloud method
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
      character*10 scmeth
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
      x3dmet(i3dcwtr) = n3d
      met_3dname(n3d) = "cloudwater"
      met_3dunits(n3d) = "g m-3"
      met_3dlong(n3d) = "cloud water content"
      met_3ddesc(n3d) = "cloud water content" 
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3drwtr) = n3d
      met_3dname(n3d) = "rainwater"
      met_3dunits(n3d) = "g m-3"
      met_3dlong(n3d) = "rain water content"
      met_3ddesc(n3d) = "rain water content"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dswtr) = n3d
      met_3dname(n3d) = "snowater"
      met_3dunits(n3d) = "g m-3"
      met_3dlong(n3d) = "snow water content"
      met_3ddesc(n3d) = "snow water content"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dgwtr) = n3d
      met_3dname(n3d) = "grplwater"
      met_3dunits(n3d) = "g m-3"
      met_3dlong(n3d) = "graupel water content"
      met_3ddesc(n3d) = "graupel water content"
      met_3dcoords(n3d) = "latitude longitude"

      n3d = n3d+1
      x3dmet(i3dcod) = n3d
      met_3dname(n3d) = "cloudod"
      met_3dunits(n3d) = "unitless"
      met_3dlong(n3d) = "cloud optical depth"
      met_3ddesc(n3d) = "cloud optical depth"
      met_3dcoords(n3d) = "latitude longitude"

      if (scmeth.eq.'KF') then
        n3d = n3d + 1
        x3dmet(i3dkf_cwtr) = n3d
        met_3dname(n3d) = "kf_cldwater"
        met_3dunits(n3d) = "g m-3"
        met_3dlong(n3d) = "KF cloud water"
        met_3ddesc(n3d) = "KF cloud water"
        met_3dcoords(n3d) = "latitude longitude"

        n3d = n3d + 1
        x3dmet(i3dkf_pwtr) = n3d
        met_3dname(n3d) = "kf_pcpwater"
        met_3dunits(n3d) = "g m-3"
        met_3dlong(n3d) = "KF precipitation water"
        met_3ddesc(n3d) = "KF precipitation water"
        met_3dcoords(n3d) = "latitude longitude"

        n3d = n3d + 1
        x3dmet(i3dkf_ent) = n3d
        met_3dname(n3d) = "kf_entrain"
        met_3dunits(n3d) = "kg m-2 s-1"
        met_3dlong(n3d) = "KF entrainment flux"
        met_3ddesc(n3d) = "KF entrainment flux"
        met_3dcoords(n3d) = "latitude longitude"

        n3d = n3d + 1
        x3dmet(i3dkf_det) = n3d
        met_3dname(n3d) = "kf_detrain"
        met_3dunits(n3d) = "kg m-2 s-1"
        met_3dlong(n3d) = "KF detrainment flux"
        met_3ddesc(n3d) = "KF detrainment flux"
        met_3dcoords(n3d) = "latitude longitude"
      endif

      if (n3d.gt.n3dmet) then
        write(*,*)
        write(*,*)'Number of output 3-D cloud variables > max in VARATT_3DCLD'
        stop
      endif

      if (ldiag) then
        n2d = n2d + 1
        x2dmet(i2dprat) = n2d
        met_2dname(n2d) = "preciprate"
        met_2dunits(n2d) = "mm hr-1"
        met_2dlong(n2d) = "surface precipitation rate"
        met_2ddesc(n2d) = "surface precipitation rate"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dtcod) = n2d
        met_2dname(n2d) = "tcloudod"
        met_2dunits(n2d) = "unitless"
        met_2dlong(n2d) = "total cloud optical depth"
        met_2ddesc(n2d) = "total cloud optical depth"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dct) = n2d
        met_2dname(n2d) = "cloudtop"
        met_2dunits(n2d) = "km"
        met_2dlong(n2d) = "convective cloud top"
        met_2ddesc(n2d) = "convective cloud top"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dcp) = n2d
        met_2dname(n2d) = "cape"
        met_2dunits(n2d) = "J kg-1"
        met_2dlong(n2d) = "convective available potential energy"
        met_2ddesc(n2d) = "convective available potential energy"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dcf) = n2d
        met_2dname(n2d) = "cloudfrac"
        met_2dunits(n2d) = "unitless"
        met_2dlong(n2d) = "Cloud cover fraction"
        met_2ddesc(n2d) = "Cloud cover fraction"
        met_2dcoords(n2d) = "latitude longitude"
      endif

      if (scmeth.eq.'KF') then
        n2d = n2d + 1
        x2dmet(i2dkf_frc) = n2d
        met_2dname(n2d) = "kf_cldfrac"
        met_2dunits(n2d) = "unitless"
        met_2dlong(n2d) = "KF cloud fraction"
        met_2ddesc(n2d) = "KF cloud fraction"
        met_2dcoords(n2d) = "latitude longitude"

        n2d = n2d + 1
        x2dmet(i2dkf_tscl) = n2d
        met_2dname(n2d) = "kf_tscale"
        met_2dunits(n2d) = "s"
        met_2dlong(n2d) = "KF cloud time scale"
        met_2ddesc(n2d) = "KF cloud time scale"
        met_2dcoords(n2d) = "latitude longitude"
      endif

      if (n2d.gt.n2dmet) then
        write(*,*)
        write(*,*)'Number of output 2-D cloud variables > max in VARATT_3DCLD'
        stop
      endif

      return
      end
