/* This routine calculates the distribution function of a 
   2-D bump-on-tail system 
   Note: when modifying this section, one must also modify gatherdf_pwc.cc*/

#include <math.h>
#include "eppic.h"

FTYPE f0dist(FTYPE vx, FTYPE vy, int id)
{

  FTYPE arg, a, b, fd, fb, f0;

  // Bulk distribution:
  arg=0;
  if (nx >1) arg += -Sqr( (vx-vx0d[id]) / (vxthd[id]*SQRT2));
  if (ny >1) arg += -Sqr( (vy-vy0d[id]) / (vythd[id]*SQRT2));
  fd = exp(arg);

  //Beam distribution:
  arg=0;
  if (nx >1 && vxthb[id]>0.) arg += -Sqr( (vx-vx0b[id]) / (vxthb[id]*SQRT2));
  if (ny >1 && vythb[id]>0.) arg += -Sqr( (vy-vy0b[id]) / (vythb[id]*SQRT2));
  fb = exp(arg);
  
  //Normalization:
  if (ny == 1) 
    a=n0d[id]/(SQRT2PI*(n0d[id]*vxthd[id]+n0b[id]*vxthb[id]));
  else if (nx == 1) 
    a=n0d[id]/(SQRT2PI*(n0d[id]*vythd[id]+n0b[id]*vythb[id]));
  else
    a=n0d[id]/(2*PI*(n0d[id]*vxthd[id]*vythd[id]+n0b[id]*vxthb[id]*vythb[id]));

  b=a*n0b[id]/n0d[id];

  f0 = a*fd + b*fb;


  return f0;
}
