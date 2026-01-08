      subroutine ncf_wrt_vars(action,iounit,numcols,numrows,numlays, &
                              nvars,varnames,varunits,varlong,vardesc, &
                              varcoord,num_dims)
      use ncf_iomod
      implicit none
!
!-----This routine writes the variable definitions and descriptions to 
!     the NetCDF file
!
!      Argument description:
!       Inputs:
!           action   C  name of file to open
!           iounit   I  NetCDF file ID of file
!           numcols  I  number of columns in the grid
!           numrows  I  number of rows in the grid
!           numlays  I  number of layers in the grid
!           nvars    I  number of variables in the file
!           varnames C  names of each variable
!           varunits C  units for each variable
!           varlong  C  long name for each variable
!           vardesc  C  description of each variable
!           varcoord C  coordinates for each variable
!           num_dims I  number of dimensions in variable
!       Outputs:
!
      include 'netcdf.inc'
!
      character*(*) action
      integer       iounit
      integer       numcols
      integer       numrows
      integer       numlays
      integer       nvars
      character*60 varnames(nvars)
      character*60 varunits(nvars)
      character*60 varlong(nvars)
      character*60 vardesc(nvars)
      character*60 varcoord(nvars)
      integer       num_dims
      integer istrln
!
      character*60 this_var
      integer      ierr, ivar, var_dimid(4), var_chunk(4), var_id
!
!-----Entry point:
!
      var_chunk(1) = INT(REAL(numcols)/NCF_CHUNK_SIZE_VAR_X)
      if( var_chunk(1) .GT. numcols ) goto 7002
      var_chunk(2) = INT(REAL(numrows)/NCF_CHUNK_SIZE_VAR_Y)
      if( var_chunk(2) .GT. numrows ) goto 7003
      var_chunk(3) = MAX(1,INT(numlays/2))
      var_chunk(4) = 1
!
      var_dimid(1) = ncf_col_dimid
      var_dimid(2) = ncf_row_dimid
      var_dimid(3) = ncf_lay_dimid
      var_dimid(4) = ncf_tstep_dimid
!
      if( num_dims .EQ.  3) then
         var_dimid(3) = ncf_tstep_dimid
         var_chunk(3) = 1
      endif
!
      do ivar=1,nvars
!
!  --- define everything for this variable ---
!
        this_var = trim(varnames(ivar))
        ncf_var_units = varunits(ivar)
        ncf_var_long_name = varlong(ivar)
        ncf_var_desc = vardesc(ivar)
        ncf_var_coordinates = varcoord(ivar)
!
!  --- define the variable ---
!
        ierr = nf_def_var(iounit,trim(this_var),  &
                          NF_FLOAT, num_dims, var_dimid, var_id)
        if( ierr .NE. NF_NOERR ) goto 7000
!
!  --- add the attributes ---
!
        ierr = nf_put_att_text(iounit,var_id,'long_name', &
                istrln(ncf_var_long_name),ncf_var_long_name)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,var_id,'units', &
                                            60,ncf_var_units)
!                        istrln(ncf_var_units),ncf_var_units)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,var_id,'var_desc', &
                  istrln(ncf_var_desc),ncf_var_desc)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,var_id,'coordinates', &
                 istrln(ncf_var_coordinates),ncf_var_coordinates)
        if( ierr .NE. NF_NOERR ) goto 7000
!
!  --- chunkers go here ... Irie! ----
!
#ifdef CHUNK
        if( ncf_compress ) then
          ierr = nf_def_var_chunking(iounit, var_id, &
                                             NF_CHUNKED, var_chunk)
          if( ierr .NE. NF_NOERR ) goto 7001
!
          ierr = nf_def_var_deflate(iounit, var_id, NCF_SHUFFLE,  &
                                      NCF_DEFLATE, NCF_DEFLATE_LEVEL )
          if( ierr .NE. NF_NOERR ) goto 7001
        endif
#endif
!
      enddo
!
      goto 9999
!
 7000 continue
      write(*,'(//,A)') 'ERROR in NCF_WRT_VARS:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot create file variable: ',trim(this_var)
      stop
!
 7001 continue
      write(*,'(//,A)') 'ERROR in NCF_WRT_VARS:'
      write(*,'(A)') trim(action)
      write(*,'(2A)') 'Cannot set chunk parameters for variable: ', &
                                    trim(this_var)
      stop
!
 7002 continue
      write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
      write(*,'(A)') trim(action)
      write(*,'(A)') 'Cannot set chunk parameters for this file.'
      write(*,'(2A,/,A)') 'The NCF_CHUNK_SIZE_VAR_X parameter causes the ', &
           'chunk value to be larger','than the number of columns.'
      write(*,'(A,I3,A,I4)') 'Number of columns: ',numcols
      stop
!
 7003 continue
      write(*,'(//,A)') 'ERROR in NCF_IODAT.INC:'
      write(*,'(A)') trim(action)
      write(*,'(A)') 'Cannot set chunk parameters for this file.'
      write(*,'(2A,/,A)') 'The NCF_CHUNK_SIZE_VAR_Y parameter causes the ', &
           'chunk value to be larger','than the number of rows.'
      write(*,'(A,I3,A,I4)') 'Number of rows: ',numrows
      stop
!
 9999 continue
      return
      end
