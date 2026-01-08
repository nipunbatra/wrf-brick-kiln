c**** NCF_SET_SPECIES_MAPPING
c
      subroutine ncf_set_species_mapping(iounit,action,
     &       species_names_model,nspecies_model,nspecies_map,lmap_array)
      use filunit
      implicit none
c
c----CAMx v7.32 250801
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine checks the inputs species list and finds the species
c   inside the file. If the species is found the mapping is addigned
c   accordingly.
c   the simulation period and that there is gridded data for each
c   tome period that will be accessed.
c
c      Copyright 1996 - 2025
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit              I NCF file ID
c         action              C string that describes file being read
c         species_names_model C arrays of species names to map into
c         nspecies_model      I number of species in list
c       Outputs:
c         lmap_array          I NetCDF variable if for each species
c         nspecies_map        I number of species mapped 
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
      include 'camx.prm'
      include 'soap.inc'
      include 'soap3.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) action
      character*10  species_names_model(*)
      integer       nspecies_model
      integer       nspecies_map
      integer       lmap_array(*)
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
      integer ispc, ierr, this_varid, iems, ilen, isoap
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over species in list ----
c
      do ispc=1,nspecies_model
         lmap_array(ispc) = 0
         this_varid = 0
c
c   --- seek this species ---
c
         ilen = istrln(species_names_model(ispc))
         ierr = nf_inq_varid(iounit,
     &               species_names_model(ispc)(:ilen),this_varid)
         if( ierr .EQ. NF_NOERR  ) then
            lmap_array(ispc) = this_varid
         else
            if( luse_soap3 ) then
               do isoap=1,NUM_SOAP_SPECIES
                  if( species_names_model(ispc) .NE. soap3_spec_names(isoap) ) cycle
                  ilen = istrln(soap2_spec_names(isoap))
                  ierr = nf_inq_varid(iounit,
     &                  soap2_spec_names(isoap)(:ilen),this_varid)
                  if( ierr .EQ. NF_NOERR  ) lmap_array(ispc) = this_varid
               enddo
            endif
          endif
      enddo
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
