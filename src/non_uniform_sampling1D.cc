#include "eppic-pic.h"
#include "eppic-misc.h"
#include "eppic-math.h"
#include "eppic-system.h"

 /** This routine lays out n_non_uniform particles along a line, in a 
     psedo-random fashion, to reproduce a 1D probablility distribution 
     funciton.

     This is done for (np particles)/n_non_uniform identical lines, uniformly
     spaced in the other directions. Currently, the only non-uniform
     direction is x.

     \todo Use a switch to choose the direction

 */

void non_uniform_sampling1D(particle_dist &pic, /**< distribution to init */
			    FTYPEAVec &prob_dist, /**< 1D probability dist */
			    int n_non_uniform, /**< num particles along 1D */
			    eppic_system &sys_params,
			    FTYPE xstart,
			    int common_seed /**< seed (id) for pseudo-random sampling */			     ) {
  
   FTYPEAVec cdf;
   calc_cdf(prob_dist,cdf);
   FTYPEAVec inv_cdf;
   invert_function(cdf,inv_cdf,n_non_uniform+1);
   
   int xprime_seed=prime(common_seed+2);
   int yprime_seed=prime(common_seed+3);
   int zprime_seed=prime(common_seed+4);
   //   int random_seed = -33;
   FTYPE x_revnum, y_revnum, pos;
#if NDIM==3
   FTYPE z_revnum;
#endif
   int ilow=0, ihigh=1;
   for (long long int ipart=0;ipart<pic.np;ipart++) {
       
     // calculate 3 reverse prime numbers
     
     

     DONDIM(x_revnum = reverse(ipart+1,xprime_seed);,		\
	    y_revnum = reverse(ipart+1,yprime_seed);,		\
	    z_revnum = reverse(ipart+1,zprime_seed););
     /*
     DONDIM(x_revnum = ran3(&random_seed);,		\
	    y_revnum = ran3(&random_seed);,		\
	    z_revnum = ran3(&random_seed););
     */
     // if non uniform direction, sample inverse function

     ilow=int(x_revnum*(n_non_uniform));
     ihigh = ilow+1;
     pos = 
       // shift if first domain and injection
       xstart + 
       // baseline x position
       inv_cdf(ilow) +
       // liniearly interpolated x position
       (inv_cdf(ihigh)-inv_cdf(ilow))*((x_revnum*(n_non_uniform))-ilow);

       DONDIM(pic.x(ipart) = pos,					\
	    pic.y(ipart) = sys_params.ny*y_revnum,			\
	    pic.z(ipart) = sys_params.nz*z_revnum)

     }
    

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
