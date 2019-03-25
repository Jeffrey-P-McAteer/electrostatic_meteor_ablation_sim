/* This is a simple routine to take advantage of
   the bin_output_ghost function in the ArrayNd 
   class. It only writes data for the current 
   time step, but it could be easily modified to
   append data to an existing file by putting the
   definition of openbyte in an IF statement. 

   Created for debugging efield_quasineut.cc and
   related files. 06Jan2017 (may)
*/

#include <cstdio>
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"

void write_local_bin(FArrayND data,const char *name)
{
  char path[128];
  long asize = sizeof(OTYPE)
    *max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1);
  char *openbtype;
  FILE *fp;
  FILE* fopensafe(char* filename, char* mode, unsigned long skip);
  openbtype="wb\0";
  sprintf(path,"%sdomain%03d/%s%02d.bin",
	  outdir,subdomain.id_number,name,mpi_rank);
  fp=fopensafe(path,openbtype,asize*((it-1)+1));
  data.bin_output(fp);
  fclose(fp);
}

void write_local_bin(FArrayND_ranged data,const char *name)
{

  int nx_ghost[2]={0,0};

  char path[128];
  long asize = sizeof(OTYPE)
    *max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1);
  char *openbtype;
  FILE *fp;
  FILE* fopensafe(char* filename, char* mode, unsigned long skip);
  openbtype="wb\0";
  sprintf(path,"%sdomain%03d/%s%02d.bin",
	  outdir,subdomain.id_number,name,mpi_rank);
  fp=fopensafe(path,openbtype,asize*((it-1)+1));
  data.bin_output_ghost(fp,nx_ghost);
  fclose(fp);
}

void write_local_bin(FArrayND_ranged data,const char *name,int *nx_ghost)
{
  char path[128];
  long asize = sizeof(OTYPE)
    *max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1);
  char *openbtype;
  FILE *fp;
  FILE* fopensafe(char* filename, char* mode, unsigned long skip);
  openbtype="wb\0";
  sprintf(path,"%sdomain%03d/%s%02d.bin",
	  outdir,subdomain.id_number,name,mpi_rank);
  fp=fopensafe(path,openbtype,asize*((it-1)+1));
  data.bin_output_ghost(fp,nx_ghost);
  fclose(fp);
}

#endif
