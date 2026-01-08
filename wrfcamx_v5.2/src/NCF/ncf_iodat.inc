c  
c     NC_IODAT.INC declares netCDF variables
c                            
c-----------------------------------------------------------------------
c   Output variable pointers:
c-----------------------------------------------------------------------
c
      integer, parameter :: i3dz = 1
      integer, parameter :: i3dp = 2
      integer, parameter :: i3dt = 3
      integer, parameter :: i3dq = 4
      integer, parameter :: i3du = 5
      integer, parameter :: i3dv = 6

      integer, parameter :: i2dts   = 1
      integer, parameter :: i2dsd   = 2
      integer, parameter :: i2dsa   = 3
      integer, parameter :: i2du10  = 4
      integer, parameter :: i2dv10  = 5
      integer, parameter :: i2dt2   = 6
      integer, parameter :: i2dsw   = 7
      integer, parameter :: i2dsm   = 8
      integer, parameter :: i2dst   = 9 
      integer, parameter :: i2dty   = 10
      integer, parameter :: i2dwrf  = 11
      integer, parameter :: i2dcmaq = 12
      integer, parameter :: i2dmyj  = 13
      integer, parameter :: i2dysu  = 14

      integer, parameter :: i3dcwtr    = 1
      integer, parameter :: i3drwtr    = 2
      integer, parameter :: i3dswtr    = 3
      integer, parameter :: i3dgwtr    = 4
      integer, parameter :: i3dcod     = 5
      integer, parameter :: i3dkf_cwtr = 6
      integer, parameter :: i3dkf_pwtr = 7
      integer, parameter :: i3dkf_ent  = 8
      integer, parameter :: i3dkf_det  = 9

      integer, parameter :: i2dprat    = 1
      integer, parameter :: i2dtcod    = 2
      integer, parameter :: i2dct      = 3
      integer, parameter :: i2dcp      = 4
      integer, parameter :: i2dcf      = 5
      integer, parameter :: i2dkf_frc  = 6
      integer, parameter :: i2dkf_tscl = 7
c
c-----------------------------------------------------------------------
c   Dimension definitions:
c-----------------------------------------------------------------------
c
c   ncf_date_time       -- number of variables in time-stamp (always 2 - 
c                          YYYYDDD and HHMMSS)
c   ncf_date_time_dimid -- NetCDF dimension ID for this dimension
c   ncf_lay             -- number of vertical layers in modeling grid
c   ncf_lay_dimid       -- NetCDF dimension ID for this dimension
c   ncf_col             -- number of columns in modeling grid
c   ncf_col_dimid       -- NetCDF dimension ID for this dimension
c   ncf_row             -- number of rows in modeling grid
c   ncf_row_dimid       -- NetCDF dimension ID for this dimension
c   ncf_var             -- number of variables in data section
c   ncf_var_dimid       -- NetCDF dimension ID for this dimension
c   ncf_tstep_dimid     -- NetCDF dimension ID for TSTEP dimension
c
c-----------------------------------------------------------------------
c
      integer ncf_date_time
      integer ncf_date_time_dimid
      integer ncf_lay
      integer ncf_lay_dimid
      integer ncf_col
      integer ncf_col_dimid
      integer ncf_row
      integer ncf_row_dimid
      integer ncf_var
      integer ncf_var_dimid
      integer ncf_tstep_dimid

      common /ncf_dom_dat/ ncf_date_time, ncf_date_time_dimid, ncf_lay,
     &                     ncf_lay_dimid, ncf_col, ncf_col_dimid, 
     &                     ncf_row, ncf_row_dimid, ncf_var, ncf_var_dimid, 
     &                     ncf_tstep_dimid
c
c-----------------------------------------------------------------------
c  Domain definition variables:
c-----------------------------------------------------------------------
c
c   ncf_x_units               -- units of coordinate values
c   ncf_x_long_name           -- long name of this variable
c   ncf_x_var_desc            -- description of this variable
c   ncf_y_units               -- units of coordinate values
c   ncf_y_long_name           -- long name of this variable
c   ncf_y_var_desc            -- description of this variable
c   ncf_layer_units           -- units of layer values
c   ncf_layer_long_name       -- long name of this variable
c   ncf_layer_var_desc        -- description of this variable
c   ncf_longitude_units       -- units of longitude values
c   ncf_longitude_long_name   -- long name of this variable
c   ncf_longitude_var_desc    -- description of this variable
c   ncf_longitude_coordnates  -- coordinates of this variable
c   ncf_latitude_units        -- units of longitude values
c   ncf_latitude_long_name    -- long name of this variable
c   ncf_latitude_var_desc     -- description of this variable
c   ncf_latitude_coordnates   -- coordinates of this variable
c
      character*60 ncf_x_units
      character*60 ncf_x_long_name
      character*60 ncf_x_var_desc
      character*60 ncf_y_units
      character*60 ncf_y_long_name
      character*60 ncf_y_var_desc
      character*60 ncf_layer_units
      character*60 ncf_layer_long_name
      character*60 ncf_layer_var_desc
      character*60 ncf_longitude_units
      character*60 ncf_longitude_long_name
      character*60 ncf_longitude_var_desc
      character*60 ncf_longitude_coordinates
      character*60 ncf_latitude_units
      character*60 ncf_latitude_long_name
      character*60 ncf_latitude_var_desc
      character*60 ncf_latitude_coordinates
c
      common /ncf_dom_chr/ ncf_x_units, ncf_x_long_name, ncf_x_var_desc, 
     &                     ncf_y_units, ncf_y_long_name, ncf_y_var_desc, 
     &                     ncf_layer_units, ncf_layer_long_name, 
     &                     ncf_layer_var_desc, ncf_longitude_units,  
     &                     ncf_longitude_long_name, ncf_longitude_var_desc, 
     &                     ncf_longitude_coordinates, ncf_latitude_units, 
     &                     ncf_latitude_long_name, ncf_latitude_var_desc, 
     &                     ncf_latitude_coordinates 
c
c-----------------------------------------------------------------------
c  Time stamp variables:
c-----------------------------------------------------------------------
c
c   ncf_tflag             -- array of date/times for each output interval
c   ncf_tflag_units       -- units of the time stamp variable
c   ncf_tflag_long_name   -- long name of the time stamp variable
c   ncf_tflag_var_desc    -- description of the time stamp variable
c   ncf_tflag_varid       -- variable ID for the time stamp variable
c   ncf_etflag            -- array of date/times for each output interval
c   ncf_etflag_units      -- units of the time stamp variable
c   ncf_etflag_long_name  -- long name of the time stamp variable
c   ncf_etflag_var_desc   -- description of the time stamp variable
c   ncf_etflag_varid      -- variable ID for the time stamp variable
c
      integer,     allocatable, dimension(:,:) :: ncf_tflag
      character*60 ncf_tflag_units 
      character*60 ncf_tflag_long_name
      character*60 ncf_tflag_var_desc
      integer      ncf_tflag_varid
      integer,     allocatable, dimension(:,:) :: ncf_etflag
      character*60 ncf_etflag_units 
      character*60 ncf_etflag_long_name
      character*60 ncf_etflag_var_desc
      integer      ncf_etflag_varid

      common /ncf_time_chr/ ncf_tflag_units, ncf_tflag_long_name,  
     &                      ncf_tflag_var_desc, 
     &                      ncf_etflag_units, ncf_etflag_long_name,
     &                      ncf_etflag_var_desc
c
      common /ncf_time_dat/ ncf_tflag_varid, ncf_etflag_varid
c
c-----------------------------------------------------------------------
c  Species name variables:
c-----------------------------------------------------------------------
c
c   ncf_var_units       -- units of variables
c   ncf_var_long_name   -- long_name of variables
c   ncf_var_desc        -- desciptions of variables
c   ncf_var_coordinates -- desciptions of variables
c
      character*60 ncf_var_units 
      character*60 ncf_var_long_name
      character*60 ncf_var_desc
      character*60 ncf_var_coordinates
c
      common /ncf_var_chr/ ncf_var_units, ncf_var_long_name, 
     &                     ncf_var_desc, ncf_var_coordinates
c
c-----------------------------------------------------------------------
c  Global attributes for IO/API (domain definition):
c-----------------------------------------------------------------------
c
c   ncf_xorig     -- X-coordinate of domain origin
c   ncf_yorig     -- Y-coordinate of domain origin
c   ncf_xorigbuf  -- X-coordinate of domain origin with nested buffer
c   ncf_yorigbuf  -- Y-coordinate of domain origin with nested buffer
c   ncf_xcell     -- cell width in X direction
c   ncf_ycell     -- cell width in Y direction
c   ncf_iutm      -- UTM zone of grid
c   ncf_istag     -- wind stagger flag
c   ncf_cproj     -- type of projection (CAMx codes)
c   ncf_gdtyp     -- type of projection (IOAPI codes)
c   ncf_xcent     -- longitude pole
c   ncf_ycent     -- latitude pole
c   ncf_p_alp     -- true latitude1
c   ncf_p_bet     -- true latitude2
c   ncf_p_gam     -- same as longitude pole
c   ncf_nlays     -- number of layers
c   ncf_nrows     -- number of rows
c   ncf_ncols     -- number of columns
c   ncf_nrowsbuf  -- number of rows with nested buffer
c   ncf_ncolsbuf  -- number of columns with nested buffer
c   ncf_nthik     -- IOAPI variable
c   ncf_nvars     -- number of variables
c   ncf_nsteps    -- number of time steps
c
      real    ncf_xorig
      real    ncf_yorig
      real    ncf_xorigbuf
      real    ncf_yorigbuf
      real    ncf_xcell
      real    ncf_ycell
      integer ncf_iutm
      integer ncf_istag
      integer ncf_cproj
      integer ncf_gdtyp
      real    ncf_xcent
      real    ncf_ycent
      real    ncf_p_alp
      real    ncf_p_bet
      real    ncf_p_gam
      integer ncf_nlays
      integer ncf_nrows
      integer ncf_ncols
      integer ncf_nrowsbuf
      integer ncf_ncolsbuf
      integer ncf_nthik
      integer ncf_nvars
      integer ncf_nsteps
c
      common /ncf_glob_dom_dat/ ncf_xorig, ncf_yorig, ncf_xorigbuf, ncf_yorigbuf,
     &                          ncf_xcell, ncf_ycell, ncf_iutm, ncf_istag,
     &                          ncf_cproj, ncf_gdtyp, ncf_xcent, ncf_ycent,
     &                          ncf_p_alp, ncf_p_bet, ncf_p_gam, ncf_nlays,
     &                          ncf_nrows, ncf_ncols, ncf_nrowsbuf, ncf_ncolsbuf,
     &                          ncf_nthik, ncf_nvars, ncf_nsteps
c
c-----------------------------------------------------------------------
c  Global attributes for IO/API (file description):
c-----------------------------------------------------------------------
c
c   ncf_cdate              -- date of file creation (YYYYJJJ)
c   ncf_ctime              -- time of file creation
c   ncf_wdate              -- date of file write (YYYYJJJ)
c   ncf_wtime              -- time of file write  
c   ncf_note               -- simulation note
c   ncf_itzon              -- time zone of simulation
c   ncf_ftype              -- IO-API file type (always 1 = CUSTOM3)
c   ncf_vgtyp              -- IO-API v-grid type (always 6 = H: m above ground) 
c   ncf_vgtop              -- IO-API grid top for sigma coordinates
c   ncf_vglvls             -- IO-API levels 
c   ncf_gdnam              -- 
c   ncf_upnam              -- 
c   ncf_sdate              -- date of simulation start (YYYYJJJ)
c   ncf_sdatec             -- date of simulation start (YYYYMMDD)
c   ncf_stime              -- time of simulation start
c   ncf_tstep              -- time step (hours times HHMMSS)
c   ncf_conventions        -- NetCDF convnetions in this file
c   ncf_history            -- history of file access
c
      integer       ncf_cdate
      integer       ncf_ctime
      integer       ncf_wdate
      integer       ncf_wtime
      character*200 ncf_note
      integer       ncf_itzon
      integer       ncf_ftype
      integer       ncf_vgtyp
      real          ncf_vgtop
      real          ncf_vglvls
      character*16  ncf_gdnam
      character*16  ncf_upnam
      integer       ncf_sdate
      integer       ncf_sdatec
      integer       ncf_stime
      integer       ncf_tstep
      character*60  ncf_conventions
      character*60  ncf_history

      common /ncf_glob_chr/ ncf_note, ncf_gdnam, ncf_upnam, 
     &                      ncf_conventions, ncf_history

      common /ncf_glob_dat/ ncf_cdate, ncf_ctime, ncf_wdate, ncf_wtime,  
     &                      ncf_itzon, ncf_ftype, ncf_vgtyp, ncf_vgtop, 
     &                      ncf_vglvls,ncf_sdate, ncf_stime, ncf_tstep,
     &                      ncf_sdatec
c
c-----------------------------------------------------------------------
c  Variables for chunking parameters:
c-----------------------------------------------------------------------
c
c   ncf_compress          -- logical flag for file compression
c   NCF_CACHESIZE         -- total size of the raw data chunk cache in MegaBytes
c   NCF_NELEMS            -- number slots in the per-variable chunk cache 
c                            (should be a prime number larger than the number
c                             of chunks in the cache)
c   NCF_PREEMPTION        -- indicates how much chunks that have been fully
c                            read are favored for preemption (0 means fully
c                            read chunks are treated no differently than other
c                            chunks, while a value of 100 means fully read
c                            chunks are always preempted before other chunks)
c   NCF_SHUFFLE           -- If non-zero, turn on the shuffle filter
c   NCF_DEFLATE           -- If non-zero, turn on the deflate filter at the
c                            level specified by the deflate_level parameter
c   NCF_DEFLATE_LEVEL     -- 0 = no deflate and 9 = slowest, but best deflate
c   NCF_CHUNK_SIZE_VAR_X  -- number of chunks in X direction
c   NCF_CHUNK_SIZE_VAR_Y  -- number of chunks in Y direction
c  
      logical ncf_compress
      integer NCF_CACHESIZE
      integer NCF_NELEMS
      integer NCF_PREEMPTION
      integer NCF_SHUFFLE
      integer NCF_DEFLATE
      integer NCF_DEFLATE_LEVEL
      integer NCF_CHUNK_SIZE_VAR_X
      integer NCF_CHUNK_SIZE_VAR_Y
c
      parameter( NCF_CACHESIZE        = 32 )
      parameter( NCF_NELEMS           = 37 )
      parameter( NCF_PREEMPTION       = 100 )
      parameter( NCF_SHUFFLE          = 0 )
      parameter( NCF_DEFLATE          = 1 )
      parameter( NCF_DEFLATE_LEVEL    = 2 )
      parameter( NCF_CHUNK_SIZE_VAR_X = 2 )
      parameter( NCF_CHUNK_SIZE_VAR_Y = 2 )

      common /ncf_comp_dat/ ncf_compress
