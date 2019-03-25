
/*
  a header file that for now defines info about the grid and boundary 
  conditions.
*/

#ifndef EPPIC_SYSTEM_H
#define EPPIC_SYSTEM_H

#include "eppic-types.h"

typedef struct {
  int ndim;
  int nx, ny, nz;    // number of grid cells 
  FTYPE dx, dy, dz;  // grid step size in meters 
  FTYPE eps;    // dielectric constant (default=8.8542e-12) 
  int boundary_type[3]; /* Type of boundary condition: */
} eppic_system;

extern int boundary_type[3]; /* Type of particle boundary condition:
			        PERIODIC = 0
			        INJECT = 1
			        OPEN = 2
			  */
enum {PERIODIC, INJECT, OPEN};

extern int field_boundary_type[3][2];  /* Type of field boundary condition:
                                       PERIODIC_FIELD = 0
                                       DIRICHLET = 1
                                       NEUMANN = 2
                                    */
enum {PERIODIC_FIELD, DIRICHLET, NEUMANN};

extern FTYPE field_boundary_values[3][2];


#endif
