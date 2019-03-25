

#include "timerFile.h"
#include "eppic-io.h"
#include "eppic-mpi.h"
#include "eppic.h"

void init_timerFile(timerFile &tf, const char *filename, 
		    int nout, int restart) {

  // setup file
  char name[256];
  char *opentype;
  if (restart == 0) opentype="w\0 ";
  else opentype="at\0";
  sprintf(name,"%sdomain000/%s",outdir,filename);
  if (mpi_rank==0) tf.ftime=fopensafe(name, opentype);
  
  // and everything else
  tf.nout = nout;
  if (restart) tf.nflushed = 0;
  else tf.nflushed=-1;

}
