
#include "timerFile.h"
#include <iostream>
#include "eppic-mpi.h"

void flush_timerFile(timerFile &tf, int it) {
  if (mpi_rank == 0) {
    if (tf.nflushed == -1) {
      // then setup headings
      
      tf.nflushed++;
      fprintf(tf.ftime,"%11s ","time");
      for (timerList::iterator iter=tf.times.begin();
	   iter != tf.times.end();
	   ++iter) {
	fprintf(tf.ftime,"%11.11s ",iter->first.c_str());
      } 
      fprintf(tf.ftime,"\n");
    }
  }
  tf.nflushed++;

  if (it == -1) it = tf.nflushed;
  if (!(it%tf.nout)) {
    if (mpi_rank == 0) {
      // time to flush
      fprintf(tf.ftime,"%11d ",it);
      for (timerList::iterator iter=tf.times.begin();
	   iter != tf.times.end();
	   ++iter) {
	fprintf(tf.ftime,"%11ld ",iter->second);
      } 
      fprintf(tf.ftime,"\n");
    }
    reset_timerFile(tf);
  }
}
