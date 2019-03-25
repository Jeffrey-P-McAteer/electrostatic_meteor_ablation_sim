#include "eppic.h"
#include "math.h"
#include "Random.h"
//inline double Sqr(double x) {return x*x;}

// Return the velocity of a particle with an average speed of vradius, a vthermal width of vradius
void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz)
{
  static RanFast rnd(233048+mpi_rank); // Initialize the random number generator
  static GaussDev Grnd(233048+1+mpi_rank); // Initialize the random number generator
  PTYPE vr= Grnd.dev()*vth + vradius;
  // We need 2 more angles to rotate the distribution around in 3D
  do { // Begin rejection loop until we find a place for this particle
    PTYPE theta2 = (M_PI) * rnd.dbl();
    PTYPE phi = (2.*M_PI) * rnd.dbl();
    vx = vr * sin(theta2)*cos(phi);
    vy = vr * sin(theta2)*sin(phi);
    vz = vr * cos(theta2);
    // This creates a gaussian ring around (vx,vy) but then 
    // rotates it around the vz axis, causing a greatly enhanzed 
    // density on this axis.  We can repair this by rejecting 
    // particles in proportion to (vr^2-vz^2)/vr^2
  } while (rnd.dbl() > sqrt((Sqr(vr)-Sqr(vz))/Sqr(vr)) );

}

