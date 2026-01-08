c**** NCF_WRT_DATA_SPECIES_CSA
c
      subroutine ncf_wrt_data_species_csa(action,iounit,num_dims,
     &          num_cols,num_rows,num_lays_grid,num_lays_conc,
     &               num_lays_out,nspcs,spec_name,cncfld,height)
      use ncf_iomod
      use grid
      use chmstry
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the data for the timestep variables to the
c    NetCDF file
c
c     Copyright 1996 - 2025
c     Ramboll
c      Argument description:
c       Inputs:
c           action        C name of file to open
c           iounit        I NetCDF file ID of file
c           num_dims      I number of dimensions in this data 
c           num_cols      I number of columns in this grid
c           num_rows      I number of rows in this grid
c           num_lays_grid I number of layers in the array
c           num_lays_conc I number of layers in the conc array
c           num_lays_out  I number of layers in the array
c           nspcs         I number of species in the file
c           cncfld        R gridded array of concentrations
c           height        R gridded array of layer heights   
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'netcdf.inc'
      include 'filunit.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       num_dims
      integer       num_cols
      integer       num_rows
      integer       num_lays_grid
      integer       num_lays_conc
      integer       num_lays_out
      integer       nspcs
      character*(*) spec_name(nspcs)
      real          cncfld(num_cols,num_rows,num_lays_conc,nspcs)      
      real          height(num_cols,num_rows,num_lays_grid)

c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*14 this_var
      integer      data_start(4), data_count(4), ispc, i, j, k
      integer      ierr, this_varid
      real,        allocatable, dimension(:,:,:) :: array_3d
      real*8,      allocatable, dimension(:,:,:) :: darray_3d
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set the position in the NetCDF variable to write ---
c
      data_start(1) = 1
      data_count(1) = num_cols
      data_start(2) = 1
      data_count(2) = num_rows
      data_start(3) = 1
      data_count(3) = num_lays_out
      data_start(4) = ncf_cur_tstep
      data_count(4) = 1
c
c  --- allocate the array (double) that will be used to write the data ---
c
      allocate( darray_3d(num_cols,num_rows,num_lays_out) )
c
c  --- first do the z variable ---
c
      do k=1,num_lays_out
         do j=1,num_rows
            do i=1,num_cols
               darray_3d(i,j,k) = DBLE( height(i,j,k) )
            enddo
         enddo
      enddo
c
c  --- get the id for the z variable and write the data ---
c
      this_var = "z"
      ierr = nf_inq_varid(iounit,"z",this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_vara_double(iounit,this_varid,data_start,data_count,darray_3d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  ---- deallocate the array ---
c
      deallocate( darray_3d )
c
c  --- allocate the array that will be used to write the data ---
c
      allocate( array_3d(num_cols,num_rows,num_lays_out) )
c
c  --- loop over all variables in this file ---
c
      do ispc=1,nspcs
c
c  --- get name for this variable, and get it's variable ID
c      in this file ---
c
         this_var = spec_name(ispc)
         call jstlft(this_var)
         ierr = nf_inq_varid(iounit,this_var(:istrln(this_var)),this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- load the data into the local array to write ---
c
         if( num_dims .EQ. 4 ) then
             do k=1,num_lays_out
                do j=1,num_rows
                   do i=1,num_cols
                      array_3d(i,j,k) = cncfld(i,j,k,ispc)
                   enddo
                enddo
             enddo
c
           ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
           if( ierr .NE. NF_NOERR ) goto 7001
         else
             do j=1,num_rows
                do i=1,num_cols
                   array_3d(i,j,1) = cncfld(i,j,1,ispc)
                enddo
             enddo
c
           ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
           if( ierr .NE. NF_NOERR ) goto 7001
         endif
c
      enddo
c
c  ---- deallocate the array ---
c
      deallocate( array_3d )
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_SPECIES:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_SPECIES:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot write data for the variable: ',
     &                                      this_var(:istrln(this_var))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
