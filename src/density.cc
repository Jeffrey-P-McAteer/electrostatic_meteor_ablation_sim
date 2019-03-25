// Calculate the density of distribution id scaled by "scaler"

#include <math.h>
#include <stdio.h>
#include "eppic.h"
#include "eppic-mpi.h"


void density(FArrayND &den, int id, particle_dist *pic, fluid *fspecie, 
	     FTYPE scaler)
{

  extern void gather(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &), 
		     const int, FTYPE);
  extern void gather_domain(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &),  		     const int, const int, FTYPE);
  extern void gather_domain_inject(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &),  		     const int, const int, FTYPE);
  
  // The standard PIC gather: 
  if (method[id] == 0) {
    if (chargeon[id] == 0) {

#ifndef USE_DOMAINS
      gather(den, INDICIES(pic[id].x, pic[id].y, pic[id].z), pic[id].np, 
	     scaler*pic[id].n0avg);
#else
      int nxmin = 0;

      if (unscale_density){
        // Option added by Glenn to unscale the density.  By default unscale_density=false
        // If true, we don't want the scaled density, pass the charge (scaler) 
        // divided by (dx*dy*dz) times n0d for the scaler in gather
        // n0d[id] is the number of real particles each simulation particle represents
        gather_domain_inject(den, INDICIES(pic[id].x, pic[id].y, pic[id].z), 
                             pic[id].np, nxmin, n0d[id]*scaler/(dx*dy*dz));
      }
      else{
        // Default behavior.  Scale the density by n0avg
        gather_domain_inject(den, INDICIES(pic[id].x, pic[id].y, pic[id].z), 
                             pic[id].np, nxmin,
                             scaler*pic[id].n0avg*(nx*nsubdomains)*ny*nz/
                             (npd[id]*nsubdomains) );
      }

      // Adjust the densities if there are inject boundaries
      if (boundary_type[0] != PERIODIC){
        // loop through x=0
        if (subdomain.id_number==0){
          for (int iy=0; iy<ny+yguard_size; iy++){
            for (int iz=0; iz<nz+zguard_size; iz++){
              den(0,iy,iz) *= 2;
            }
          }
        }
        if (subdomain.id_number==nsubdomains-1){
          for (int iy=0; iy<ny+yguard_size; iy++){
            for (int iz=0; iz<nz+zguard_size; iz++){
              den(nx,iy,iz) *= 2;
            }
          }          
        }
      }      

      if (boundary_type[1] != PERIODIC){
        // loop through y=0 and y=ny+yguard_size planes
        for (int ix=0; ix<nx+1; ix++){
          for (int iz=0; iz<nz+zguard_size; iz++){
            den(ix,0,iz) *= 2;
            den(ix,ny+yguard_size-1,iz) *= 2;
          }
        }
      }

      if (boundary_type[2] != PERIODIC){
        for (int ix=0; ix<nx+1; ix++){
          for (int iy=0; iy<ny+yguard_size; iy++){
            den(ix,iy,0) *= 2;
            den(ix,iy,nz+zguard_size-1) *= 2;
          }
        }
      }

#endif // USE_DOMAINS
    } // end if chargeon[id]==0
  } else if (method[id] == -2) {
    // The fluid density exists - copy it over 
    den =  fspecie[id].den;
    den *= scaler;

  } else if (method[id] == -3) {
    // The fluid density exists - copy it over 
#ifdef USE_DOMAINS
    void sub_array(FArrayND_ranged &in,FArrayND &out,
		   INDICIES(int startx, int starty, int startz),
		   INDICIES(int endbeforex, int endbeforey, int endbeforez));
    den = 0.;
    //    if (subdomain.rank == 0) {
      sub_array(fspecie[id].den,den,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz));
      //    }
#else
    den =  fspecie[id].den;
#endif
    den *= scaler;
  } else if (method[id] == -4) {
    // Do nothing (constant fluid parameters only)
  } else {
    char message[130];
    sprintf(message,"Method=%d not implemented",method[id]);
    terminate(-1,message);
  }

}
