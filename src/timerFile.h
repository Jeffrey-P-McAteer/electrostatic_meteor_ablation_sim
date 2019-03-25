
#ifndef EPPIC_TIMERFILE_H
#define EPPIC_TIMERFILE_H

#include <string>
#include <sys/times.h>
#include <ctime>
#include <map>
#include <cstdio>

typedef std::map<std::string,clock_t> timerList;

typedef struct {
  FILE* ftime; // file to write times to
  timerList times; // a list of string,clock_t pairs
  tms refTimer; // reference time
  int nout; // how frequently to flush to file
  int nflushed; // how many times flushed 
} timerFile;

// a routine to open the file
void init_timerFile(timerFile &, const char *, int nout=1, int restart=0);

// a routine to write times to file at time step it
void flush_timerFile(timerFile &, int it);

// a routine to zero all the times of a timerFile
void reset_timerFile(timerFile &);
#endif
