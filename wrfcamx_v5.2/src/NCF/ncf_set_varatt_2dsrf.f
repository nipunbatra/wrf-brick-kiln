      subroutine ncf_set_varatt_2dsrf(n2dsrf,n2d,met_name,met_units,
     &                                met_long_name,met_desc,met_coords)
      implicit none
c
c-----This routine sets the output variable attributes
c
c      Argument description:
c       Inputs:
c       Outputs:
c            n2dsrf        I max number of 2-D met variables
c            n2d           I number of output 2-D met variables
c            met_name      C array of variable names
c            met_units     C array of units
c            met_long_name C array of "long names"
c            met_desc      C array of desciptions
c            met_coords    C array of coordinates
c
      integer n2dsrf,n2d
      character*60 met_name(*)
      character*60 met_units(*)
      character*60 met_long_name(*)
      character*60 met_desc(*)
      character*60 met_coords(*)
c
c-----Entry point:
c
      n2d = 1
      met_name(n2d) = "topo"
      met_units(n2d) = "m MSL"
      met_long_name(n2d) = "topographic elevation"
      met_desc(n2d) = "topographic elevation m above sea level"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "water"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "1 water (ocean)"
      met_desc(n2d) = "1 water (ocean)"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "ice"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "2 ice"
      met_desc(n2d) = "2 ice"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "lake"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "3 lake (fresh)"
      met_desc(n2d) = "3 lake (fresh)"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "eneedl"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "4 evergreen needle leaf forest"
      met_desc(n2d) = "4 evergreen needle leaf forest"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "ebroad"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "5 evregreen broad leaf forest"
      met_desc(n2d) = "5 evregreen broad leaf forest"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "dneedl"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "6 deciduous needle leaf forest"
      met_desc(n2d) = "6 deciduous needle leaf forest"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "dbroad"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "7 deciduous broad leaf forest"
      met_desc(n2d) = "7 deciduous broad leaf forest"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "tbroad"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "8 tropical broad leaf forest"
      met_desc(n2d) = "8 tropical broad leaf forest"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "ddecid"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "9 drought deciduous trees"
      met_desc(n2d) = "9 drought deciduous trees"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "eshrub"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "10 evergreen shrub"
      met_desc(n2d) = "10 evergreen shrub"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "dshrub"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "11 deciduous shrub"
      met_desc(n2d) = "11 deciduous shrub"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "tshrub"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "12 thorn shrub"
      met_desc(n2d) = "12 thorn shrub"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "sgrass"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "13 short grass"
      met_desc(n2d) = "13 short grass"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "lgrass"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "14 long grass"
      met_desc(n2d) = "14 long grass"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "crops"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "15 cropland"
      met_desc(n2d) = "15 cropland"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "rice"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "16 rice"
      met_desc(n2d) = "16 rice"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "sugar"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "17 sugar"
      met_desc(n2d) = "17 sugar"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "maize"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "18 maize"
      met_desc(n2d) = "18 maize"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "cotton"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "19 cotton"
      met_desc(n2d) = "19 cotton"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "icrops"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "20 irrigated cropland"
      met_desc(n2d) = "20 irrigated cropland"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "urban"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "21 urban"
      met_desc(n2d) = "21 urban"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "tundra"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "22 tundra"
      met_desc(n2d) = "22 tundra"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "swamp"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "23 swamp"
      met_desc(n2d) = "23 swamp"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "desert"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "24 desert (barren)"
      met_desc(n2d) = "24 desert (barren)"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "mwood"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "25 mixed woodland"
      met_desc(n2d) = "25 mixed woodland"
      met_coords(n2d) = "latitude longitude"

      n2d = n2d+1
      met_name(n2d) = "tforest"
      met_units(n2d) = "fraction"
      met_long_name(n2d) = "26 transitional forest"
      met_desc(n2d) = "26 transitional forest"
      met_coords(n2d) = "latitude longitude"

      if (n2d.gt.n2dsrf) then
        write(*,*)
        write(*,*)'Number of output 2-D surface variables > max in VARATT_2DSRF'
        stop
      endif

      return
      end
