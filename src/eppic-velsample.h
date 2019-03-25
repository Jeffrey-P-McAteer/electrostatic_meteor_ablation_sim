


#ifndef EPPIC_VELSAMPLE_H
#define EPPIC_VELSAMPLE_H


#include "eppic-types.h"

typedef struct {
  int ran3Seed;
  ArrayNd<FTYPE,1> inv_cdf;
  PTYPE dv;
} velsample;

void init_velsample(velsample &veldev, PTYPE vth, PTYPE v0, 
		    int precision, int direction, int seed);
void init_velsample_ablating(velsample &veldev, PTYPE vth, PTYPE v0, 
                             int precision, int direction, int seed);

float rand_normal(unsigned long long int iran);

PTYPE sampleVel(velsample &veldev);

#endif
