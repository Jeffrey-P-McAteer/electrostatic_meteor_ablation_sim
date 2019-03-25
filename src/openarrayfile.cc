// A routine to attach the directory name (from outdir), open, 
// advance to the correct position, if necessary,  and return a file pointer.

#include <stdio.h> 
#include "eppic.h"

FILE* openarrayfile(char* fileroot, int it, int subcycle)
{
    FILE* fopensafe(char* filename, char* mode);
    FILE* fopensafe(char* filename, char* mode, unsigned long skip);
    unsigned long asize=sizeof(OTYPE)*max(nx/nout_avg,1)*max(ny/nout_avg,1)*
      max(nz/nout_avg,1);

    char name[256], *openbtype="wb\0 ";
    if (it > 0) openbtype="ab\0";

    sprintf(name,"%s%s.bin",outdir,fileroot);
    FILE* fname=fopensafe(name,openbtype,
			  asize*((static_cast<unsigned long>(it)-1)/
				 (nout*subcycle)+1));

    return fname;
   }


