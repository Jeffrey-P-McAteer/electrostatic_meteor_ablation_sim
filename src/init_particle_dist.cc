// Initialize the paticle quantities 
#include "eppic-pic.h"
#include "eppic-mpi.h"
#include "eppic-system.h"
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
			)
{
    // padded number of particles
  int npad = (int)(npd*part_pad); 
  if (npd == 0){
    npad = (int)(part_pad); 
  }
  
    // Load particle characteristics 
    pic.m=md;
    pic.q=qd;
    pic.n0avg=n0d;
    pic.np=npd;
    pic.np_all = npd*nsubdomains;
    
    // Injection: If this is the leftmost domain then we must 
    // include enough particles for the extra cell on the left:
    if (sys_params.boundary_type[0]==INJECT && subdomain.id_number==0) {
      FTYPE scale = (sys_params.nx*nsubdomains+1);
      scale /= sys_params.nx*nsubdomains;
      pic.np=static_cast<int>(npd*scale);
    }

    long long int np_all = npd*nsubdomains;
    if ((sys_params.boundary_type[0] == INJECT)) {
      if (subdomain.id_number == 0) 
	np_all += pic.np - npd;
      else
	np_all += 
	  static_cast<int>(npd*(sys_params.nx*nsubdomains+nsubdomains)
			   /(sys_params.nx*nsubdomains))-npd;
    }
    pic.np_all = np_all;
      
    // Initiailize particle position arrays:
    pic.x=PTYPEAVec(npad) = 0.;
    if (sys_params.ndim >= 2) pic.y=PTYPEAVec(npad) =0.; 
    if (sys_params.ndim == 3) pic.z=PTYPEAVec(npad) = 0.;

    // Allocate the absent array to keep track of particles
    // that have migrated out of the domain
    if (nsubdomains > 1 || sys_params.boundary_type[0] == INJECT) {
      pic.absent=intAVec((int) (pic.np*(part_pad-1.0)))=-1;
      pic.nabsent=0;
    }

    // Initiailize particle velocity arrays:
    pic.vx=PTYPEAVec(npad) = 0.;
    if (sys_params.ndim >= 2) pic.vy=PTYPEAVec(npad) = 0.;
    if (sys_params.ndim == 3 || vel_dim == 3) pic.vz=PTYPEAVec(npad) = 0.;
    


}
 
