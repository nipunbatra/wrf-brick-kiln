      subroutine ncf_set_varatt_kv(mxkv,nkv,met_name,met_units,met_long_name,
     &                             met_desc,met_coords)
      implicit none
c
c-----This routine sets the output variable attributes
c
c      Argument description:
c       Inputs:
c       Outputs:
c            mxkv          I max number of Kv variables
c            nkv           I number of output Kv variables
c            met_name      C array of variable names
c            met_units     C array of units
c            met_long_name C array of "long names"
c            met_desc      C array of desciption
c            met_coords    C array of coordinates
c
      integer mxkv,nkv
      character*60 met_name(*)
      character*60 met_units(*)
      character*60 met_long_name(*)
      character*60 met_desc(*)
      character*60 met_coords(*)
c
c-----Entry point:
c
      nkv = 1
      met_name(1) = "kv"
      met_units(1) = "m2 s-1"
      met_long_name(1) = "vertical diffusivity"
      met_desc(1) = "vertical diffusivity"
      met_coords(1) = "latitude longitude"

      if (nkv.gt.mxkv) then
        write(*,*)
        write(*,*)'Number of output 3-D Kv variables > max in VARATT_KV'
      endif
      return
      end
