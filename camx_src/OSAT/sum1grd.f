      subroutine sum1grd(ncola,nrowa,numcols,numrows,numlays,
     &                  nspmod,nspmod_tot,nsptrac,igroup,igrd,idx,
     &                    this_species,emssum,emsgrd,emsbas,emsoth,
     &                                        emslft,emstot,lemit)
      use grid
      use chmstry
      use tracer
c
c
c      Copyright 1996 - 2025
c     Ramboll
c
c
c----CAMx v7.32 250801
c
c     SUM1GRD sums up the area emission of one species for a given group
c     in a given grid
c
c       07/19/02  --gwilson-- Added seperate source area map for each grids.
c       10/28/09  --gwilson-- Changed dimension of variables to accomodate
c                             the dynamic memory allocation
c       03/01/16  --gwilson-- Added partial source area map
c
c     Input argument:
c        numcols           number of columns in any grid
c        numrows           number of columns in any grid
c        numlays           number of columns in any grid
c        nspmod            number of model species
c        nspmod_tot        number of model species, including SOAP3 emissions species
c        nsptrac           number of tracer species
c        igroup            group ID
c        igrd              grid ID
c        idx               species ID
c        this_species      species name
c        emsgrd            the species emission in the grid
c
c     Output arguments:
c        emssum            emission summed over grid
c        emsbas            base emission
c        emsoth            "otherwise" emission
c        emslft            leftover emission
c        emstot            total emission
c        lemit             flag to determine if tracer class is emitted
c        
      include "camx.prm"
      include "soap.inc"
      include "soap3.inc"
c
      character*10 this_species
      integer   ncola
      integer   nrowa
      integer   numcols
      integer   numrows
      integer   numlays
      integer   nspmod
      integer   nspmod_tot
      integer   nsptrac
      integer   idx
      real      emssum(nspmod_tot,nsptrac)
      real      emsgrd(numcols,numrows,numlays)
      real      emsbas(nspmod_tot,nsptrac)
      real      emsoth(nspmod_tot,nsptrac)
      real*8    emslft(ncola,nrowa,nspmod_tot)
      real*8    emstot(ncola,nrowa,nspmod_tot)
      logical   lemit(*)
      integer   ipart
      real      frac
c
      logical   luse
c
c  --- Entry point ---
c
      luse = .FALSE.
      if( idx .LE. nspmod ) then
         do icls=1,ntrcls
           if( trspmap(idx,icls) .NE. 0.  .OR.
     &                         yhratmap(idx,icls) .NE. 0. .OR.
     &                         ylratmap(idx,icls) .NE. 0. ) then
             if( trspmap(idx,icls) .NE. 0. ) lemit(icls) = .TRUE. 
             luse = .TRUE.
           endif
         enddo
      endif
      if( idx .GT. nspmod  ) then
         do icls=ntrcls+1,soap3_ntrcls
             lemit(icls) = .TRUE. 
             luse = .TRUE.
         enddo
      endif
      if( .NOT. luse ) goto 9999
c
c   --- sum up the emissions excluding the boundary for this grid ---
c
      do 40 j=2,nrow(igrd)-1
        jcrs = (j - 2)/nmesh(igrd) + j1(igrd)
        if( igrd .eq. 1 ) jcrs = j
        do 50 i=2,ncol(igrd)-1
c
c   --- calculate the coarse grid offset ---
c
           icrs = (i - 2)/nmesh(igrd) + i1(igrd)
           if( igrd .eq. 1 ) icrs = i
c
c   --- skip cell if this grid has a child in this cell ---
c
           ijcl = i + (j-1)*ncol(igrd)
           if( idfin(iptr2d(igrd)-1+ijcl) .NE. 0 ) goto 50
c
c  --- get the region for this cell from mapping array,
c      the grid cell should be the coarse grid  ----
c
           do ipart=1,npartial(igroup,igrd)
              imap = igrmap(igroup,ipart,igrd,i,j)
              frac = frcmap(igroup,ipart,igrd,i,j)
              if( (imap .LE. 0 .AND. frac .GT. 0.) .OR. imap .GT. nregin ) goto 50
c
c  --- calculate the index into the tracer species for this gruoup/region ---
c
              if( ngroup .GT. 0 ) then
c
c   --- if group is base emissions, add to "leftover" group ----
c
                 if( igroup .EQ. 0 ) then
                   if( leftovr_area ) then
                      do icls=1,ntrcls
                        if( idx .LE. nspmod ) then
                           if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                              ipt = iemcls(icls)-1 + imap+ngroup*nregin
                              do k=1,numlays
                                emsbas(idx,ipt) = emsbas(idx,ipt) + emsgrd(i,j,k) * frac
                              enddo
                           endif
                        endif
                      enddo
                      if( idx .GT. nspmod ) then
                         isoa = nspmod
                         do icls=ntrcls+1,soap3_ntrcls
                           isoa = isoa + 1
                           iclass = icls - ntrcls + (IDX_SOAP_POA_OP-1)
                           ipt = soap3_iemcls(iclass)-1 + imap+ngroup*nregin
                           if( idx .EQ. isoa ) then
                              do k=1,numlays
                                  emsbas(isoa,ipt) = emsbas(isoa,ipt) +
     &                                                emsgrd(i,j,k) * frac
                              enddo
                           endif
                         enddo
                      endif
                   endif
                   do icls=1,ntrcls
                     if( idx .LE. nspmod ) then
                         if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                           do k=1,numlays
                             emstot(icrs,jcrs,idx) = 
     &                                  emstot(icrs,jcrs,idx) + emsgrd(i,j,k) * frac
                           enddo
                         endif
                     endif
                   enddo
                   if( idx .GT. nspmod ) then
                     isoa = nspmod
                     do icls=ntrcls+1,soap3_ntrcls
                         isoa = isoa + 1
                         if( idx .EQ. isoa ) then
                            do k=1,numlays
                              emstot(icrs,jcrs,isoa) = 
     &                                  emstot(icrs,jcrs,isoa) + emsgrd(i,j,k) * frac
                            enddo
                         endif
                      enddo
                   endif
c
c   --- otherwise, add to this group/region and subtract from "leftover" ---
c
                 else
                   do icls=1,ntrcls
                     if( idx .LE. nspmod ) then
                        if( trspmap(idx,icls) .NE. 0. .OR. 
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                          ipt = iemcls(icls)-1 + imap+(igroup-1)*nregin
                          do k=1,numlays
                              emssum(idx,ipt) = emssum(idx,ipt) + emsgrd(i,j,k) * frac
                              if( leftovr_area ) then
                                 ipt = iemcls(icls)-1 + imap+ngroup*nregin
                                 emsoth(idx,ipt) = emsoth(idx,ipt) + emsgrd(i,j,k) * frac
                              endif
                           enddo
                        endif
                     endif
                   enddo
                   if( idx .GT. nspmod ) then
                     isoa = nspmod
                     iclass = IDX_SOAP_POA_OP-1
                     do icls=ntrcls+1,soap3_ntrcls
                        iclass = iclass + 1
                        isoa = isoa + 1
                        if( this_species .NE. soap_emiss_names(iclass) ) cycle 
                        ipt = soap3_iemcls(iclass)-1 + imap+(igroup-1)*nregin
                        do k=1,numlays
                           emssum(isoa,ipt) = emssum(isoa,ipt) + emsgrd(i,j,k) * frac
                           if( leftovr_area ) then
                              ipt = soap3_iemcls(iclass)-1 + imap+ngroup*nregin
                              emsoth(isoa,ipt) = emsoth(isoa,ipt) + emsgrd(i,j,k) * frac
                           endif
                        enddo
                      enddo
                   endif
                   do icls=1,ntrcls
                     if( idx .LE. nspmod ) then
                        if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                          do k=1,numlays
                             emslft(icrs,jcrs,idx) = 
     &                               emslft(icrs,jcrs,idx) + emsgrd(i,j,k) * frac
                          enddo
                        endif
                     endif
                   enddo
                   if( idx .GT. nspmod ) then
                      iclass = IDX_SOAP_POA_OP-1
                      isoa = nspmod
                      do icls=ntrcls+1,soap3_ntrcls
                        iclass = iclass + 1
                        if( this_species .NE. soap_emiss_names(iclass) ) cycle 
                        do k=1,numlays
                             emslft(icrs,jcrs,idx) = emslft(icrs,jcrs,idx) +
     &                                          emsgrd(i,j,k) * frac
                        enddo
                      enddo
                   endif
                 endif
c
c   --- only using regular model emissions ---
c
              else
                 do icls=1,ntrcls
                   if( idx .LE. nspmod ) then
                      if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                         ipt = iemcls(icls) - 1 + imap
                         do k=1,numlays
                            emssum(idx,ipt) = emssum(idx,ipt) + emsgrd(i,j,k) * frac
                         enddo
                      endif
                   endif
                 enddo
                 if( idx .GT. nspmod ) then
                    isoa = nspmod
                    iclass = IDX_SOAP_POA_OP-1
                    do icls=ntrcls+1,soap3_ntrcls
                      iclass = iclass + 1
                      isoa = isoa + 1
                      if( this_species .NE. soap_emiss_names(iclass) ) cycle 
                      ipt = soap3_iemcls(iclass) - 1 + imap
                      do k=1,numlays
                         emssum(isoa,ipt) = emssum(isoa,ipt) + emsgrd(i,j,k) * frac
                      enddo
                    enddo
                 endif
              endif
           enddo
  50    continue
  40  continue
c
 9999 continue
      return
      end
