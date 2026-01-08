      subroutine ncf_set_varatt_grid()
      implicit none
c
c-----This routine sets the variable definitions and descriptions for
c     the NetCDF file
c
c      Argument description:
c       Inputs:
c       Outputs:
c
      include 'ncf_iodat.inc'
c
c-----Entry point:
c
      if (ncf_cproj.eq.0) then
        ncf_x_units = "degrees"
        ncf_x_var_desc = "Longitude degrees east"
        ncf_y_units = "degrees"
        ncf_y_var_desc = "Latitude degrees north"
      else
        ncf_x_units = "km"
        ncf_x_var_desc = "X cartesian distance from projection origin"
        ncf_y_units = "km"
        ncf_y_var_desc = "Y cartesian distance from projection origin"
      endif
      ncf_x_long_name = "X coordinate"
      ncf_y_long_name = "Y coordinate"
      ncf_layer_units = "Layer index"
      ncf_layer_long_name = "Model layer"
      ncf_layer_var_desc = "Model layer"
      ncf_longitude_units = "Degrees east"
      ncf_longitude_long_name = "Longitude"
      ncf_longitude_var_desc = "Longitude degrees east"
      ncf_longitude_coordinates = "latitude longitude"
      ncf_latitude_units = "Degrees north"
      ncf_latitude_long_name = "Latitude"
      ncf_latitude_var_desc = "Latitude degrees north"
      ncf_latitude_coordinates = "latitude longitude"
c
      ncf_tflag_units = "YYYYDDD,HHMMSS"
      ncf_tflag_long_name = "Start time flag"
      ncf_tflag_var_desc = "Timestep start date and time"
      ncf_etflag_units = "YYYYDDD,HHMMSS"
      ncf_etflag_long_name = "End time flag"
      ncf_etflag_var_desc = "Timestep end date and time"
c
      return
      end
