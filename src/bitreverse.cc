/* This routine generates a repeatable bit reversed sequence of base, b, for 
   a more uniform distribution than true random numbers would generate.  
   n should be prime */

#include "eppic-types.h"

FTYPE bitreverse(int i, int b)
{
  static PTYPE step, start, t;
  PTYPE x;
  static int istore=-1;
  
  if (i != istore) {
    step=1./b;
    t=step;
    start=step;
  }
      
  x=t;
  t=t+step;
  if (t>=1) {
    start /= b;
    t=start;
  }

  istore=i+1;
  return x;
}

  
    
