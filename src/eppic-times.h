
/* 
   a header file containing all global variables related to clocking each
   routine
*/


#ifndef EPPIC_TIMES_H
#define EPPIC_TIMES_H
#include "eppic-types.h"

//Clock related functions: Let's gather timing information
#include <sys/times.h> 
#include <ctime>
extern tms times_buf;
extern clock_t 
  vadvance_time,     // Time spent scattering and advancing velocities
  xadvance_time,     // Time spent advancing positions
  charge_time,       // Time spent gathering densities
  collect_time,      // Time spent collecting densities across processors
  efield_time,       // Time spent calculating electric field
  output_time,       // Time spent in output routine
  fluid_time;        // Time spent in fluid routine

extern time_t sim_start_time;
extern FTYPE time_limit;
// NEED TO MAKE AN INITIALIZATION ROUTINE!!!


#endif

