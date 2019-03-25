/*  Shuffle a set of numbers */
#include "eppic.h"

void shuffle (PTYPEAVec &x, int &seed)
{
  FTYPE ran3(int *);
  PTYPE y;
  int j;
    
  for (int i=0; i<x.length(); i++) {
    do{
      j= (int)( ran3(&seed)*x.length() );
    } while (j<0 || j>=x.length());
    y=x[j];
    x[j]=x[i];
    x[i]=y;
  }
}

      
    
