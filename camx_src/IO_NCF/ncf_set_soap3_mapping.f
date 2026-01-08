c**** NCF_SET_SOAP3_MAPPING
c
      subroutine ncf_set_soap3_mapping(iounit,action,
     &        species_names_model,nspecies_model,nspecies_map,lmap_array)
      use filunit
      use chmstry
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
c   accordingly. This version is for the new species names introduced
c   with SOAP3.
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
      include 'netcdf.inc'
      include 'soap3.inc'
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
      integer ispc, isoap, ierr, this_varid, iems, ilen
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over species in list ----
c
      do ispc=1,nspecies_model
         this_varid = 0
         if( lmap_array(ispc) .GT. 0 ) cycle
c
c   --- loop over new names for match of SOAP3 species ---
c
         do isoap=1,NUM_SOAP_SPECIES
            if( species_names_model(ispc) .EQ. soap3_spec_names(isoap) ) then 
c
c   --- first look for new species names, these have already been set ---
c
               ilen = istrln(soap2_spec_names(isoap))
               ierr = nf_inq_varid(iounit,
     &                  soap2_spec_names(isoap)(:ilen),this_varid)
               if( ierr .EQ. NF_NOERR ) then
                   lmap_array(ispc) = this_varid
                   nemspc = nemspc + 1
                   emspcname(nemspc) = soap2_spec_names(isoap)
                   lemmap(ispc) = nemspc
               endif
            endif
         enddo
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
 
