/* Advance the 1D unmagnatized particle velocities to the next time step
with either local or global damping. The equation modeled is:
dv
-- = q E - nu v ( v - v0 )/vth
dt
The input constant c is usually defined as c = damp_nu / vth * dx
 */

#if NDIM == 2
#include <math.h>
#include "eppic.h"
//#include <signal.h>
#include <stdio.h>

void vpush1D_bc1(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     particle_misc &misc, FTYPE c, FTYPE qmrat, FTYPE alpha)
{
  
  FTYPE Exscale, Exiscale;
  int i;
  int np = old.x.length();
  int ixl, ixh;
  int iyl, iyh;
  FTYPE wxh;
  FTYPE wyh;
  FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;
  FTYPE vx,vxmag;

  FTYPE charge_frac;

  /* Rescale the electric field so that it may simply be added to
	particle velocities in the particle loop. */
  Exscale =  alpha * qmrat * (dt*dt/dx) ;
  Exiscale = 1. / Exscale;
  Efield.Ex *= Exscale;

     /* For each particle, update the velocity ... */
  for (i = 0; i < np; ++i) {

    charge_frac=misc.charge[i];
    if ( charge_frac != 0) {

      /* Define the nearest grid points: */
      ixl = (int) cur.x(i);
      ixh = ixl + 1;
      if (ixh == nx) ixh = 0;
      
      iyl = (int) cur.y(i);
      iyh = iyl + 1;
      if (iyh == ny) iyh = 0;
      
      /* And the corresponding linear weighting factors: */
      wxh = cur.x(i) - ixl;
      wyh = cur.y(i) - iyl;

      /* Calculate Exi, the linearly-weighted Ex at the particle: */
      Exll = Efield.Ex(ixl,iyl);
      Exlh = Efield.Ex(ixl,iyh);
      Exhl = Efield.Ex(ixh,iyl);
      Exhh = Efield.Ex(ixh,iyh);
      Exl = Exll + wyh*(Exlh-Exll);
      Exh = Exhl + wyh*(Exhh-Exhl);
      Exi = Exl + wxh*(Exh-Exl);
      
      Exi *= charge_frac;

      /* Advance the particle velocities: */
      if (cur.x(i) < damp_nx) {
	adv.vx(i) = old.vx(i) + Exi;
      } else {
	/* Add damping - we're in the damped region: */
	vx=cur.vx(i);
	vxmag=fabs(vx);
	adv.vx(i) = old.vx(i)+Exi-c*vxmag*(vx-misc.vx0(i));
#ifdef CHECK
      // Particle damping can violate the Courant conditions.
      // Let's look for outrageous velocities.  
      if (adv.vx(i) > 128) {
	printf("Error: particle %d velocity = %g too fast\n",i,adv.vx(i));
	printf("Check damping rate (possibly too high)\n");
	raise(SIGTRAP);
      }
#endif      
      }

    } else {  /* if (charge_frac == 0) */
      adv.vx(i)=old.vx(i);
    }
    
    
  } /* End for */

  /*  Restore original values to the electric field array: */

  Efield.Ex *= Exiscale;
  
  return ;
  
} /* vpush */
#endif
