
#ifndef EPPIC_FLUID_H
#define EPPIC_FLUID_H
#include "eppic-types.h"

/* // A structure defining fluid parameters (no arrays) */
/* // -->Intended to save space when using quasineutral electrons. */
/* //    If we keep this, we could have the fluid struct inherit  */
/* //    the fparams struct. Probably a minor detail. */
/* typedef struct { */
/*   FTYPE m, q, n0;        // mass, charge and number density  */
/*   FTYPE gamma, T0;       // gamma (isothermal=1, adiabatic=5/3), Temperature  */
/*   FTYPE nu, diffc;       // collision rate, hyperdiffusion coeff (1.0 is effective) */
/* } fparams; */

#ifdef USE_DOMAINS
// A structure defining a fluid 
typedef struct { 
  FTYPE m, q, n0;        // mass, charge and number density 
  FTYPE gamma, T0;       // gamma (isothermal=1, adiabatic=5/3), Temperature 
  FTYPE nu, diffc;       // collision rate, hyperdiffusion coeff (1.0 is effective)
  FArrayND_ranged den;          // array of densities normalized to n0
  FArrayND_ranged vx, vy, vz;   // array of velocities   
  //  FArrayND den_1, den_2; // Intermediate density (or other) arrays, if needed.
} fluid;

#else



// A structure defining a fluid 
typedef struct { 
  FTYPE m, q, n0;        // mass, charge and number density 
  FTYPE gamma, T0;       // gamma (isothermal=1, adiabatic=5/3), Temperature 
  FTYPE nu, diffc;       // collision rate, hyperdiffusion coeff (1.0 is effective)
  FArrayND den;          // array of densities normalized to n0
  FArrayND vx, vy, vz;   // array of velocities   
  //  FArrayND den_1, den_2; // Intermediate density (or other) arrays, if needed.
} fluid;


#endif // USE_DOMAINS

#endif
