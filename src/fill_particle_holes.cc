//---------------------------------------------------------------------

// Assign indices to particles that have migrated into local subdomain.
// First check "absent" array for any holes, and fill these.  If more
// locations are required, increase np and add the additional particles
// to the end of the list.

// "iabsent" is index of last particle in "absent" array

// Meers Oppenheim & Doug Sondak
// 8/21/06

//---------------------------------------------------------------------

#include "eppic.h"
#include "eppic-mpi.h"

// Add optional absent variable to add particles but label them absent (Glenn 3/14/2018)
void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim, 
                         bool add_absent=false)
{

  int n, i, ind;
  int tot_dim = ndim + vel_dim;
  
  // particles imported from recv

  for(n=0; n<nrecv; n+=tot_dim){
    // Check for holes in current "absent" array,
    // starting from last value.  If there's a hole,
    // fill it.  If not, add to end of list.
    
    ind = n;
    if (adv.nabsent>0) {
      adv.nabsent -= 1;
      i=adv.absent(adv.nabsent);
      adv.absent(adv.nabsent) = -1; // Not necessary but cleaner
      // adv.np++;
    } else {
      i = adv.np++;
    }
  
    if (adv.np > adv.x.size()){
      printf("\n adv.np = %d , adv.x.size() = %d \n", adv.np, adv.x.size());
      terminate (-1,"Error: Array sizes need to be larger in fill_particle_holes.  Increase part_pad.");
    }
    
    adv.x(i)  = recv(ind++);
    if(ndim > 1) adv.y(i)  = recv(ind++);
    if(ndim > 2) adv.z(i)  = recv(ind++);
    
    adv.vx(i) = recv(ind++);
    if(vel_dim > 1) adv.vy(i) = recv(ind++);
    if(vel_dim > 2) adv.vz(i) = recv(ind++);
    if (add_absent){
      adv.absent(adv.nabsent++) = i;
    }
  }
}
