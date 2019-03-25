#include "eppic-pic.h"
#include "eppic-misc.h"
#include "eppic-math.h"
#include "eppic-system.h"
#include "eppic-mpi.h"

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
			      int common_seed/**< seed (id) for pseudo-random sampling */			     ) {
  
  FTYPEAVec cdfx,cdfy;
  calc_cdf(prob_distx,cdfx);
  calc_cdf(prob_disty,cdfy);
  FTYPEAVec inv_cdfx,inv_cdfy;
  invert_function(cdfx,inv_cdfx,npic_x+1);
  invert_function(cdfy,inv_cdfy,npic_y+1);

  int xprime_seed=prime(common_seed+2);
  int yprime_seed=prime(common_seed+3);
  int zprime_seed=prime(common_seed+4);
  //   int random_seed = -33;
  FTYPE x_revnum, y_revnum, posx, posy;
#if NDIM==3
  FTYPE z_revnum;
#endif
   int ilow=0, ihigh=1;
   int curr_part=0;
   FTYPE local_xrange=right_area-left_area;
   long long int ipart=0;


   for (long long int i=1;i<total_particles;i++) {
     ipart=i+total_particles*subdomain.rank+1;
     // calculate 3 reverse prime numbers
     
     

     DONDIM(x_revnum = reverse(ipart+1,xprime_seed);,		\
	    y_revnum = reverse(ipart+1,yprime_seed);,			\
	    z_revnum = reverse(ipart+1,zprime_seed););
     /*
     DONDIM(x_revnum = ran3(&random_seed);,		\
	    y_revnum = ran3(&random_seed);,		\
	    z_revnum = ran3(&random_seed););
     */
     // if non uniform direction, sample inverse function

     if ((x_revnum>=left_area)&&(x_revnum<right_area)) {

       ilow=int((x_revnum-left_area)/local_xrange*(npic_x));
       ihigh = ilow+1;
       posx = 
	 // shift if first domain and injection
	 xstart + 
	 // baseline x position
	 inv_cdfx(ilow) +
	 // liniearly interpolated x position
	 (inv_cdfx(ihigh)-inv_cdfx(ilow))*(((x_revnum-left_area)/local_xrange*(npic_x))-ilow);
       
       
       ilow=int(y_revnum*(npic_y));
       ihigh = ilow+1;
       posy =
	 // baseline x position
	 inv_cdfy(ilow) +
	 // liniearly interpolated y position
	 (inv_cdfy(ihigh)-inv_cdfy(ilow))*((y_revnum*(npic_y))-ilow);
       
       DONDIM(pic.x(curr_part) = posx,					\
	      pic.y(curr_part) = posy,					\
	      pic.z(curr_part) = sys_params.nz*z_revnum);
       
       curr_part++;
      
	 
     }
    
   }
   pic.np = curr_part-1;
}


/* old code that needs to be considered

   // loop over each particle
   FTYPE x_revnum, y_revnum, z_revnum;
   // quicker to do switch outside particle loop
   switch(direction){
   case 1: // x direction
     for (long long int ipart=0;ipart<pic.np;ipart++) {
       
       // calculate 3 reverse prime numbers

       DONDIM(x_revnum = reverse(ipart,xprime_seed);,\
	      y_revnum = reverse(ipart,yprime_seed);,		\
	      z_revnum = reverse(ipart,zprime_seed););
       // if non uniform direction, sample inverse function

       DONDIM(pic.x(ipart) = inv_cdf(int(floor(x_revnum*n_non_uniform))), \
	      pic.y(ipart) = ny*y_revnum,				\
	      pic.z(ipart) = nz*z_revnum)

     }
     break;
   case 2: // y direction
     for (int ipart=0;ipart<pic.np;ipart++) {

       // calculate 3 reverse prime numbers

       DONDIM(x_revnum = reverse(ipart,xprime_seed),		\
	      y_revnum = reverse(ipart,yprime_seed),		\
	      z_revnum = reverse(ipart,zprime_seed))

       // if non uniform direction, sample inverse function

       DONDIM(pic.y(ipart) = inv_cdf(int(floor(y_revnum*n_non_uniform)));, \
		pic.x(ipart) = nx*x_revnum;,				\
		pic.z(ipart) = nz*z_revnum;);


     }
     break;
   case 3: // z direction
     for (int ipart=0;ipart<pic.np;ipart++) {
       // calculate 3 reverse prime numbers

       DONDIM(x_revnum = reverse(ipart,xprime_seed),	 \
		y_revnum = reverse(ipart,yprime_seed),  \
		z_revnum = reverse(ipart,zprime_seed));

       // if non uniform direction, sample inverse function
       DONDIM(pic.z(ipart) = inv_cdf(int(floor(z_revnum*n_non_uniform)));, \
		pic.y(ipart) = ny*y_revnum;,				\
		pic.x(ipart) = nx*x_revnum;);

		}
     break;
     }
   
*/
