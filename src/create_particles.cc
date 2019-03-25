//----------------------------------------------------

// create particles particles throughout the system
// by placing new particles into the absent arrays.

// Meers Oppenheim June 2009
// 

//----------------------------------------------------

#include "eppic.h"
#include "eppic-mpi.h"
#include "math.h"
#include "Random.h"

void create_particles(int id, particle_dist &adv, FTYPE creation_rate)
{

  //initialize the random number generator:
  const unsigned long long int iran = 1042830+mpi_rank;
  static RanFast rnd(iran);

  static PTYPE fnx = (PTYPE) nx;
  PTYPE smallest_float_below(PTYPE x);
  static PTYPE fnx_minus_epsilon=smallest_float_below(fnx);

  // Number of particles to create this time step is:
  FTYPE n_part_create = creation_rate*(nx*dx*ny*dy*nz*dz)/dt;
  int int_part_create=static_cast<int>(n_part_create);
  // The fractional particle will contribute:
  if ( (n_part_create-int_part_create) > (rnd.dbl()) ) int_part_create += 1;
  // Only a shell of high energy particles are currently implemented...
  
  int ndim_place=ndim+vel_dim[id];
  PTYPEAVec place(int_part_create*ndim_place);
  int nplace=0;
  for (int i=0; i<int_part_create; i++) {
    // Position the particles uniformally through the system.
    place(nplace++)=rnd.dbl()*nx;
    if(ndim >= 2) place(nplace++)=rnd.dbl()*ny;
    if(ndim >= 3) place(nplace++)=rnd.dbl()*nz;
    
    // Calculating the velocity distribution is done in a seperate routine:
    void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz);
    PTYPE vx,vy,vz;
    create_shell_particle(vx0d[id],  vxthd[id], vx, vy, vz);

    place(nplace++)=vx*dt/dx;
    if (vel_dim[id] >= 2) place(nplace++)=vy*dt/dy;
    if (vel_dim[id] >= 3) place(nplace++)=vz*dt/dz;
  }

  void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                           bool add_absent=false);
  fill_particle_holes(place, nplace, adv, vel_dim(id));  


}

// Create a pair of particles, typically one co-located electron and ion:
void create_particles(int id1, int id2, particle_dist *adv, FTYPE creation_rate)
{


  if (id2 == -1 || id1 == id2 || id2 >= ndist) 
    terminate(-1,"create_pair_dist does not give appropriate id");

  //initialize the random number generator:
  const unsigned long long int iran = 1042830+mpi_rank;
  static RanFast rnd(iran);

  static PTYPE fnx = (PTYPE) nx;
  PTYPE smallest_float_below(PTYPE x);
  static PTYPE fnx_minus_epsilon=smallest_float_below(fnx);

  // Number of particles to create this time step is:
  FTYPE n_part_create = creation_rate*(nx*dx*ny*dy*nz*dz)/dt;
  int int_part_create=static_cast<int>(n_part_create);
  // The fractional particle will contribute:
  if ( (n_part_create-int_part_create) > (rnd.dbl()) ) int_part_create += 1;
  // Only a shell of high energy particles are currently implemented...
  
  int ndim_place=ndim+vel_dim[id1];
  PTYPEAVec place(int_part_create*ndim_place);
  int nplace=0;
  for (int i=0; i<int_part_create; i++) {
    // Position the particles uniformally through the system.
    place(nplace++)=rnd.dbl()*nx;
    if(ndim >= 2) place(nplace++)=rnd.dbl()*ny;
    if(ndim >= 3) place(nplace++)=rnd.dbl()*nz;
    
    // Calculating the velocity distribution is done in a seperate routine:
    void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz);
    PTYPE vx,vy,vz;
    create_shell_particle(v0_shell[id1],  vth_shell[id1], vx, vy, vz);

    place(nplace++)=vx*dt/dx;
    if (vel_dim[id1] >= 2) place(nplace++)=vy*dt/dy;
    if (vel_dim[id1] >= 3) place(nplace++)=vz*dt/dz;

  }

  void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                           bool add_absent=false);
  fill_particle_holes(place, int_part_create, adv[id1], vel_dim(id1));  

  // The second distribution will have particles placed in the same places 
  // but with equal and opposite momenta plus the momentum of the neutral created.

//  FTYPE gasdev(int &);
//  static int iran2 = -3389 - mpi_rank;
  static GaussDev vel(3954+1+mpi_rank);
  FTYPE inv_mass = 1.0/md[id2];
  FTYPE scalex=dt/dx;
  FTYPE scaley=dt/dy;
  FTYPE scalez=dt/dz;
  nplace=0;
  for (int i=0; i<int_part_create; i++) {
    
    nplace += NDIM; // Use the same positions as the original particle
    PTYPE vn=(vxthd_neutral[id2]*vel.dev() + vx0d_neutral[id2]) * scalex;
    place(nplace) = (m_neutral*vn-md[id1]*place(nplace))*inv_mass;
    nplace++;
    if (vel_dim[id1] >= 2) {
      vn=(vythd_neutral[id2]*vel.dev() + vy0d_neutral[id2]) * scaley;
      place(nplace) = (m_neutral*vn-md[id1]*place(nplace))*inv_mass;
      nplace++;
    }
    if (vel_dim[id1] >= 3) {
      vn=(vzthd_neutral[id2]*vel.dev() + vz0d_neutral[id2]) * scalez;
      place(nplace) = (m_neutral*vn-md[id1]*place(nplace))*inv_mass;
      nplace++;
    }

  }

  fill_particle_holes(place, int_part_create, adv[id2], vel_dim(id2));
}
