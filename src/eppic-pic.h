

#ifndef EPPIC_PIC_H
#define EPPIC_PIC_H
#include "eppic-types.h"

/// A structure defining particle distributions 
typedef struct { 
  long int np;            ///< Number of particles in species 
  long long int np_all;   ///< Number of particles across all domains
  PTYPE m,q,n0avg;        ///< mass, charge and number density 
  PTYPEAVec x, y, z;      ///< array of x positions (normalized to dx)
  PTYPEAVec vx, vy, vz;   ///< array of velocity positions
                          ///< (normalized to dt/dx, dt/dy and dt/dy)  
  intAVec absent;         ///< List of indicies of absent particles
  int nabsent;            ///< highest valid index of absent
} particle_dist;


/** A structure containing miscellaneous information relating to a 
    distribution of particles.
    Used for velocity damping information
    Used for charge variability information */
typedef struct { 
  PTYPEAVec vx0;/**< v values (normalized to dx/dt) to restore to 
		    when damping when Bx0 != 0 then vy0 stores v_perp*/
  PTYPEAVec vy0;/**< v values (normalized to dx/dt) to restore to 
		   when damping when Bx0 != 0 then vy0 stores v_perp*/
  PTYPEAVec vz0; /**< v values (normalized to dx/dt) to restore to 
		    when damping when Bx0 != 0 then vy0 stores v_perp*/
  PTYPEAVec charge;  /**< If the distribution requires each particle to
			have a distinct charge store it here (the 
			final charge is still multiplied by pic[id].q
			and the algorithms assume charge=0 or 1 (for now)*/
} particle_misc;


/** defines the particles and other miscellaneous parameters */
void init_particles(particle_dist *&pic, particle_misc *&misc);

/** This routine lays out n_non_uniform particles along a line, in a 
    psedo-random fashion, to reproduce a 1D probablility distribution 
    funciton.

    This is done for (np particles)/n_non_uniform identical lines, uniformly
    spaced in the other directions. direction defines which direction
    should be non-uniform

*/

#include "eppic-system.h"

void non_uniform_sampling1D(particle_dist &pic, /**< distribution to init */
			    FTYPEAVec &prob_dist, /**< 1D probability dist */
			    int n_non_uniform, /**< num particles along 1D */
			    eppic_system &sys_params,
			    FTYPE xstart,
			    int common_seed /**< seed (id) for pseudo-random sampling */			     );


 /** This routine lays out particles non-uniformly according to the product
     of two 1D probability distribution functions.

     \todo Use a switch to choose the direction

 */

void non_uniform_sampling2x1D(particle_dist &pic, /**< distribution to init */
			      FTYPEAVec &prob_distx,/**< 1D probability dist */
			      FTYPEAVec &prob_disty,/**< 1D probability dist */
			      int npic_x,/**< num particles along x 1D */
			      int npic_y,/**< num particles along y 1D */
			      FTYPE right_area,
			      FTYPE left_area,
			      long long int total_particles,
			      eppic_system &sys_params,
			      FTYPE xstart,
			      int common_seed/**< seed (id) for pseudo-random sampling */			     );


/** A routine to initialize particle distributions. This routine only prepares the
    arrays for the particles, it does not load them. See load_particle_dist
*/

void init_particle_dist(particle_dist &pic,
			FTYPE md,
			FTYPE qd,
			FTYPE n0d,
			long int npd,
			FTYPE part_pad,
			int vel_dim,
			eppic_system &sys_params
			);

/** This routine initializes a particle distribution so that its
    velocity distribution matches a Maxwellian distribution in 
    3D. 
    
    Currently, it assumes that the thermal velocities are equal
    in all directions. It uses a pseudo-random number generator for
    each direction, because this hopefully limits correlations between
    the particle positions and their velocities. 
*/

void maxwel_vel_dist_loading(particle_dist &pic, 
			     int vel_dim,
			     FTYPE vthx, FTYPE vthy, FTYPE vthz,
			     int common_seed);



#endif
