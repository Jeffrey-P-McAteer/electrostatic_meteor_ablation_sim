#include "eppic-types.h"
#include "eppic-math.h"

float smallest_float_above(float x)
{
  // Calculates the smallest distinct float less than x

  // 
  float epsilon=fabs(x);

  do {
    epsilon /= 2.0;
      }
  while ( !(float(x+epsilon/2.0) ==  x) ); 
  return x+epsilon;
}

double smallest_float_above(double x)
{
  // Calculates the smallest distinct double less than x

  // 
  double epsilon=fabs(x);

  do {
    epsilon /= 2.0;
      }
  while ( !(double(x+epsilon/2.0) ==  x) ); 
  return x+epsilon;
}
