      subroutine ncf_wrt_vars_grid(action,iounit)
      implicit none
c
c-----This routine writes the variable definitions and descriptions to 
c     the NetCDF file
c
c      Argument description:
c       Inputs:
c           action  C  name of file to open
c           iounit  I  NetCDF file ID of file
c       Outputs:
c
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
      character*(*) action
      integer       iounit
      integer       istrln
c
      integer ierr, grid_dimid(2), time_dimid(3), this_varid
c
c-----Entry point:
c
      grid_dimid(1) = ncf_col_dimid
      grid_dimid(2) = ncf_row_dimid
c
      time_dimid(1) = ncf_date_time_dimid
      time_dimid(2) = ncf_var_dimid
      time_dimid(3) = ncf_tstep_dimid
c
      ierr = nf_def_var(iounit, "X", NF_DOUBLE, 1, 
     &                                      ncf_col_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                                    istrln(ncf_x_units),ncf_x_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                              istrln(ncf_x_long_name),ncf_x_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                                istrln(ncf_x_var_desc),ncf_x_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "Y", NF_DOUBLE, 1, 
     &                                      ncf_row_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                                    istrln(ncf_y_units),ncf_y_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                            istrln(ncf_y_long_name),ncf_y_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                             istrln(ncf_y_var_desc),ncf_y_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "layer", NF_DOUBLE, 1, 
     &                                ncf_lay_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                           istrln(ncf_layer_units),ncf_layer_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                   istrln(ncf_layer_long_name),ncf_layer_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                     istrln(ncf_layer_var_desc),ncf_layer_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "TFLAG", NF_INT, 3, 
     &                                  time_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                          istrln(ncf_tflag_units),ncf_tflag_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                  istrln(ncf_tflag_long_name),ncf_tflag_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                     istrln(ncf_tflag_var_desc),ncf_tflag_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000

      ierr = nf_def_var(iounit, "ETFLAG", NF_INT, 3, 
     &                                 time_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                        istrln(ncf_etflag_units),ncf_etflag_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                istrln(ncf_etflag_long_name),ncf_etflag_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                  istrln(ncf_etflag_var_desc),ncf_etflag_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "longitude", NF_DOUBLE, 2, 
     &                              grid_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                  istrln(ncf_longitude_units),ncf_longitude_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &          istrln(ncf_longitude_long_name),ncf_longitude_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &             istrln(ncf_longitude_var_desc),ncf_longitude_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &        istrln(ncf_longitude_coordinates),ncf_longitude_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "latitude", NF_DOUBLE, 2, 
     &                                grid_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                   istrln(ncf_latitude_units),ncf_latitude_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &          istrln(ncf_latitude_long_name),ncf_latitude_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &             istrln(ncf_latitude_var_desc),ncf_latitude_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &       istrln(ncf_latitude_coordinates),ncf_latitude_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
 7000 continue
      write(*,'(//,a)') 'ERROR in NCF_WRT_VARS_GRID:'
      write(*,'(A)') trim(action)
      write(*,'(A)') 'Cannot create file attributes for grid variables.'
      stop
c
 9999 continue
      return
      end
