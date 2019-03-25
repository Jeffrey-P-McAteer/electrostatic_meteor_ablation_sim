
#ifndef EPPIC_MATH_H
#define EPPIC_MATH_H

#include "eppic-types.h"
#include <math.h>

inline FTYPE max(FTYPE x,FTYPE y) {return x > y ? x : y;}
inline FTYPE min(FTYPE x,FTYPE y) {return x < y ? x : y;}
inline int max(int x,int y) {return x > y ? x : y;}
inline int min(int x,int y) {return x < y ? x : y;}

inline double Sqr(double x) {return x*x;}
inline float Sqr(float x) {return x*x;}
inline FTYPE sign(FTYPE x) {return (x<0?-1:1);}

int prime(int);
FTYPE gasdev(int &);
FTYPE reverse(long long int,int);

#include "Random.h"
FTYPE ran3(int *idum);

#ifndef PI
#define PI  3.1415926535897932385
#endif


float smallest_float_above(float x);
double smallest_float_above(double x);


void transpose2d(ArrayNd<FTYPE,2> &in_array);

#endif
