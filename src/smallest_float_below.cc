#include "eppic.h"

float smallest_float_below(float x)
{
  
  // Calculates the smallest distinct float less than x

  // 
  float epsilon=x;

  do {
    epsilon /= 2.0;
      }
  while ( !(float(x-epsilon/2.0) ==  x) ); 
  
  return x-epsilon;
}

double smallest_float_below(double x)
{
  
  // Calculates the smallest distinct double less than x

  // 
  double epsilon=x;

  do {
    epsilon /= 2.0;
      }
  while ( !(double(x-epsilon/2.0) ==  x) ); 
  
  return x-epsilon;
}
