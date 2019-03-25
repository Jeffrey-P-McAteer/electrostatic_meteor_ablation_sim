/* This will calculate beta (the ionization probability) based on various
   models for a neutral-neutral ionizing collision.  The inputs must be
   the relative collision velocity in km/sec, and an integer related to
   the model the user wants to use.
*/

#include <csignal> /* for notify_and_trap */
#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"
      
FTYPE get_beta(FTYPE vel, int model)
{
  enum {UNITY, JONES, VONDRAK};

  if (model == UNITY) {
    return 1.0;
  } else if (model == JONES) {
    FTYPE vterm = 9.4*pow(10,-6)*pow(vel-10,2)*pow(vel,0.8);
    return vterm/(1.0+vterm);
  } else if (model == VONDRAK) {
    FTYPE vterm = 0.933*pow(vel-8.86,2)/pow(vel,1.94);
    return vterm;
  }
  
  printf("WARNING: No appropriate beta model elected, assuming all collisions are ionizing!\n");

  return 1.0;
}
