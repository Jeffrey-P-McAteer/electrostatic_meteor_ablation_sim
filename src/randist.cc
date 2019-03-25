/* Generate a random distribution of particle under the distribution passed as
   (*func) using the rejection method 
   (*func) must apply on the domain [0,1] and have a maximum value of 1.*/
#include "eppic.h"
#include "eppic-types.h"

FTYPE randist(FTYPE (*func)(FTYPE, int), int id, int &iran)
  /* id is passed on to the distribution function */
{
  
  FTYPE x, y, z;
  extern FTYPE ran3(int *);

  do {
  /* Choose a value within the range */
  x=ran3(&iran);
  y=ran3(&iran);
  z=(*func)(x,id);
  } while ( y > z);
    
  return x;
  
}

/* A 2D version of ranvel */
void randist(FTYPE (*func)(FTYPE, FTYPE, int),  FTYPE &x, FTYPE &y, 
	     int id,int &iran)
  /* id is passed on to the distribution function */
{
  
  FTYPE p, z;
  extern FTYPE ran3(int *);

  do {
  /* Choose a value within the range */
  x=ran3(&iran);
  y=ran3(&iran);
  p=ran3(&iran);
  z=(*func)(x, y, id);
  } while ( p > z);
  
}
