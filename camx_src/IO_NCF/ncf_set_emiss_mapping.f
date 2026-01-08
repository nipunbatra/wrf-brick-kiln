c**** NCF_SET_EMISS_MAPPING
c
      subroutine ncf_set_emiss_mapping(iounit,action)
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
c   This routine checks which species are in the file and adds them to
c   the master emissions species list. It also sets the array for
c   the index of species into this file.
c
c      Copyright 1996 - 2025
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit        I NCF file ID
c         action        C string that describes file being read
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
      include 'camx.prm'
      include 'vbs.inc'
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
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
      integer NUMVBS
      parameter( NUMVBS = 11 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*20  this_var(NUM_SOAP_SPECIES), this_mod_species
      integer ispc, ierr, this_varid, idx_hoa, idx_svoa, i
      logical lfound
c
      integer isoa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      this_var = soap3_spec_names
      idx_hoa = -9
      idx_svoa = -9
c
c   --- loop ober species in list ----
c
      do ispc=1,nspec
         this_mod_species = spname(ispc)
         if( this_mod_species(1:4) .EQ. 'HOA ' ) then
             this_mod_species = 'POA       '
             idx_hoa = ispc
         endif
         if( this_mod_species(1:5) .EQ. 'SVOA ' ) idx_svoa = ispc
c
c   --- seek this species ---
c
         do isoa=1,NUM_SOAP_SPECIES
             if( this_mod_species .EQ. soap3_spec_names(isoa) ) 
     &                    this_mod_species = soap2_spec_names(isoa) 
         enddo
         ierr = nf_inq_varid(iounit,this_mod_species,this_varid)
         if( ierr .EQ. NF_NOERR .OR. this_mod_species(1:4) .EQ. 'POA '  ) then
c
c   --- now check if this species is already in master list ---
c
           lfound = .FALSE.
           do i=1,nemspc
             if(this_mod_species .EQ. emspcname(i) ) then
                lfound = .TRUE.
                lemmap(ispc) = i
             endif
           enddo
           if( .NOT. lfound ) then
             nemspc = nemspc + 1
             emspcname(nemspc) = this_mod_species
             lemmap(ispc) = nemspc
           endif
         endif
c
c  --- next species ---
c
      enddo
c
c  --- update mapping for HOA species
c
      if( luse_soap3 ) then
         do isoa=IDX_SOAP_POA_OP,IDX_SOAP_POA_BB
           lfound = .FALSE.
            do i=1,nemspc
              if(soap_emiss_names(isoa) .EQ. emspcname(i) ) 
     &                                            lfound = .TRUE.
            enddo
            if( .NOT. lfound ) then
              nemspc = nemspc + 1
              emspcname(nemspc) = soap_emiss_names(isoa)
            endif
         enddo
         lfound = .FALSE.
         do i=1,nemspc
            if(emspcname(i)(1:4) .EQ. 'POA ' ) then
              lfound = .TRUE.
              lemmap(idx_hoa) = i
            endif
         enddo
         if( .NOT. lfound ) then
            nemspc = nemspc + 1
            emspcname(nemspc) = 'POA       '
            lemmap(idx_hoa) = nemspc
         endif
         lfound = .FALSE.
         do i=1,nemspc
            if(emspcname(i)(1:5) .EQ. 'SVOA ' ) then
              lfound = .TRUE.
              lemmap(idx_svoa) = i
            endif
         enddo
         if( .NOT. lfound ) then
            nemspc = nemspc + 1
            emspcname(nemspc) = 'SVOA       '
            lemmap(idx_svoa) = nemspc
         endif
      endif
c
c  --- do VBS or SOAP3 as a separate case ---
c
      if( lvbs .AND. LVBSPREPROC ) then
         do ispc=1,NUMVBS
            ierr = nf_inq_varid(iounit,this_var(ispc),this_varid)
            if( ierr .NE. NF_NOERR  ) cycle
            lfound = .FALSE.
            do i=1,nemspc
              if(this_var(ispc) .EQ. emspcname(i) ) then
                 lfound = .TRUE.
                 lemmap(ispc) = i
                 if( luse_soap3 .AND. idx_hoa .GT. 0 ) lemmap(idx_hoa) = i
              endif
            enddo
            if( .NOT. lfound ) then
              nemspc = nemspc + 1
              emspcname(nemspc) = this_var(ispc)
              lemmap(ispc) = nemspc
              if( luse_soap3 .AND. idx_hoa .GT. 0 ) lemmap(idx_hoa) = nemspc
            endif
         enddo
      endif
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
 
