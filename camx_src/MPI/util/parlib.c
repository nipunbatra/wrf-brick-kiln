#include "util.h"

void par_assoc_buff_(buff,numbuff)
     int *numbuff;
     void *buff;
{
  extern char *ibuff;
  extern int ipos,nbuff;

  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;
}

void par_get_float8_(words,numwords)
     long *numwords;
     float *words;
{
  extern char *ibuff;
  extern int ipos,nbuff;
  int ierr=0;

  ierr= MPI_Unpack(ibuff,nbuff,&ipos,words,*numwords
		   ,MPI_FLOAT,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_get_float-%d \n",ierr);
}

void par_get_float_(words,numwords)
     int *numwords;
     float *words;
{
  extern char *ibuff;
  extern int ipos,nbuff;
  int ierr=0;

  ierr= MPI_Unpack(ibuff,nbuff,&ipos,words,*numwords
		   ,MPI_FLOAT,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_get_float-%d \n",ierr);
}

void par_get_int_(iwords,numwords)
     int *iwords,*numwords;
{
  extern char *ibuff;
  extern int ipos,nbuff;
  int ierr=0;

  ierr = MPI_Unpack(ibuff,nbuff,&ipos,iwords,*numwords
	     ,MPI_INT,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_get_int-%d \n",ierr);

}

extern int par_store_tag(MPI_Request msgtag);
extern int zero_tag_array();

void par_get_noblock_(void *buff,int *numbuff,int *mmtype
                               ,int *ihostnum,int *fmsgtag)
{
   MPI_Request msgtag;
   int ierr=0, nbuff;

   if (flag_msgtag == 0) ierr = zero_tag_array();

   nbuff=*numbuff*sizeof(float);
   ipos=0;

   ierr=MPI_Irecv(buff,nbuff,MPI_PACKED,*ihostnum,*mmtype
                 ,MPI_COMM_WORLD,&msgtag);

   *fmsgtag=par_store_tag(msgtag);

   if(ierr < 0) 
    printf("Error in par_get_noblock\n");

}

void par_init_flag_( )
{
   flag_msgtag = 0;
}

void par_init_put_(char *buff,int *numbuff) 
{
  extern char *ibuff;
  extern int ipos, nbuff;

  /*ibuff= buff;*/
  ibuff=(char*)buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0 ;

}

void par_put_float8_(words,numwords)
     long *numwords;
     float *words;
{
  extern char *ibuff;
  extern int ipos, nbuff;
  int ierr=0;

  ierr = MPI_Pack( words,*numwords,MPI_FLOAT,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_put_float8- %d \n",ierr);
}

void par_put_float_(words,numwords)
     int *numwords;
     float *words;
{
  extern char *ibuff;
  extern int ipos, nbuff;
  int ierr=0;

  ierr = MPI_Pack( words,*numwords,MPI_FLOAT,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_put_float- %d \n",ierr);
}

void par_put_int_(iwords,numwords)
     int *iwords, *numwords;
{
  extern char *ibuff;
  extern int ipos,nbuff;
  int ierr=0, isize=0;

  ierr = MPI_Pack( iwords,*numwords,MPI_INT,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);

  if(ierr < 0) printf("Error in par_put_int- %d %d %d \n",ierr,isize,ipos);
}

MPI_Request par_retrieve_tag(MPI_Request fmsgtag)
{
   extern int flag_msgtag;
   extern MPI_Request mpi_msgtags[MAX_MSGTAG];

   MPI_Request i;

   /*For Safety*/

/*** cgwilson Tue 01 Jul 2025 07:16:17 AM PDT

   if (fmsgtag >= MAX_MSGTAG) 
      printf("ERROR:par_retrieve_tag:fmsgtag >= MAX_MSGTAG\n");

***/

   if (fmsgtag >= MAX_MSGTAG) 
      printf("ERROR: retrieving message tag.\n");
    
   i = mpi_msgtags[fmsgtag];

   mpi_msgtags[fmsgtag] = 0;

   return (i);
}

extern int par_store_tag(MPI_Request msgtag);
extern int zero_tag_array();

void par_send_noblock_(int *mach, int *msgtype, int *fmsgtag)
{

   extern char *ibuff;
   extern int ipos, nbuff;

   MPI_Request msgtag;

   int ierr=0;

   if (flag_msgtag == 0) ierr = zero_tag_array();

   ierr= MPI_Isend( ibuff, ipos, MPI_PACKED, *mach, *msgtype, MPI_COMM_WORLD
    ,&msgtag);

   *fmsgtag=par_store_tag(msgtag);

/*** cgwilson Tue 01 Jul 2025 11:57:11 AM PDT

   if(ierr < 0) 
    printf("Error in par_send_noblock - %d %d %d \n",*mach, *msgtype,ierr);

****/

   if(ierr < 0) 
    printf("Error in par_send_noblock\n");

}

int par_store_tag(MPI_Request msgtag)
{
   extern int flag_msgtag;
   extern MPI_Request mpi_msgtags[MAX_MSGTAG];

   int i;

   for(i=0; i < MAX_MSGTAG; i++) {
      if (mpi_msgtags[i] == 0) {
         mpi_msgtags[i]= msgtag;
         return(i);
      }
   }

/*** cgwilson Tue 01 Jul 2025 11:56:00 AM PDT

  printf("par_store_tag error:msgtag,max_msgtag= %d %d\n",msgtag,MAX_MSGTAG);

***/

  printf("par_store_tag error.\n");

  return (i);
}

extern MPI_Request par_retrieve_tag(MPI_Request msgtag);
extern int zero_tag_array();

void par_wait_(int *fmsgtag,int *ibytes,int *msgtype,int *ihostnum)
{
   MPI_Request msgtag;
   MPI_Status status;
   int ierr=0;

   if (flag_msgtag == 0) ierr = zero_tag_array();
   
   msgtag=par_retrieve_tag(*fmsgtag);
   
   ierr=MPI_Wait(&msgtag,&status);
  
   if(ierr < 0) 
      printf("Error in par_wait\n");
    
   MPI_Get_count(&status,MPI_PACKED,ibytes);
   *msgtype=status.MPI_TAG;
   *ihostnum=status.MPI_SOURCE;
   
}

int zero_tag_array()
{
  extern int flag_msgtag;
  extern MPI_Request mpi_msgtags[MAX_MSGTAG];

  int i;
  
  flag_msgtag = 1;

  for(i=0;i<MAX_MSGTAG;i++)
    {
      mpi_msgtags[i] = 0;
    }
  return (0);

}
