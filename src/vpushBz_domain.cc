/* 3 routines below for 1-D 2-D and 3-D to
   Advance the magnetized particle velocities to the next time step. 
   Assumes a constant Bz.*/
#include <cmath>
#include "eppic.h"
#include "interpolate_efield.h"

#if NDIM == 1
void vpushBz_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
                    FTYPE qmrat, FTYPE alpha, FTYPE Bz)
{
  
  // Just use the unmagnetized routine:
  void vpush_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, 
		    field &Efield, FTYPE qmrat, FTYPE alpha);

  vpush_domain(adv, old, cur, Efield, qmrat, alpha);
  
}
#endif


#if NDIM == 2
void vpushBz_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
                    FTYPE qmrat, FTYPE alpha, FTYPE Bz)
{
  
  const PTYPE nxmin=0; // This is the minimum acceptable value for particle positions

  FTYPE Exscale;
  FTYPE Eyscale;
  int i;
  int np = old.np;
  int ixl, ixh;
  int iyl, iyh;
  FTYPE wxh;
  FTYPE wyh;
  FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;
  FTYPE Eyll, Eylh, Eyhl, Eyhh, Eyl, Eyh, Eyi;

  // A vpush routine for non-zero Bz 
  FTYPE wcidt, bt, bs;
  FTYPE vya, vxa, vxb, vyc, vxc;
    
  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop.  
  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Exscale =  alpha * qmrat * (dt*dt/dx) / (-2*dx);

  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Eyscale =  0.5 * alpha * qmrat * (dt*dt/dy) / (-2*dy);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2*dx);
  FTYPE Ey0dy=Ey0_external*(-2*dy);

  // Rotation information - Boris constants 
  wcidt = alpha * dt/2 * (qmrat*Bz);
  bt = wcidt;
  bs = 2*bt/(1+bt*bt);
  // Scale constants to put velocity variables into appropriately normalized system
  bt *= dy/dx;
  bs *= dx/dy;
      
  // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {

    if (cur.x(i) >= nxmin) { // This excludes absent particles which should always be <0

      // Define the nearest grid points: 
      ixl = (int) cur.x(i);
      ixh = ixl + 1;

      iyl = (int) cur.y(i);
      iyh = iyl + 1;
      if (iyh == ny) iyh = 0;

      // And the corresponding linear weighting factors: 
      wxh = cur.x(i) - ixl;
      wyh = cur.y(i) - iyl;

      // Calculate Exi, the linearly-weighted Ex at the particle: 
      int ixl_m1 = ixl-1;
      //    if (ixl_m1 < 0) ixl_m1 = nx-1;
      int ixl_p1=ixh;
      int ixh_m1 = ixl;
      int ixh_p1=ixh+1;
    
      // Calculate Exi, the linearly-weighted Ex at the particle: 
      // Calculate Exi, the linearly-weighted Ex at the particle: 
      Exll = (Efield.phi_rho(ixl_p1,iyl) - Efield.phi_rho(ixl_m1,iyl)+Ex0dx)*Exscale;
      Exlh = (Efield.phi_rho(ixl_p1,iyh) - Efield.phi_rho(ixl_m1,iyh)+Ex0dx)*Exscale;
      Exhl = (Efield.phi_rho(ixh_p1,iyl) - Efield.phi_rho(ixh_m1,iyl)+Ex0dx)*Exscale;
      Exhh = (Efield.phi_rho(ixh_p1,iyh) - Efield.phi_rho(ixh_m1,iyh)+Ex0dx)*Exscale;
      Exl = Exll + wyh*(Exlh-Exll);
      Exh = Exhl + wyh*(Exhh-Exhl);
      Exi = Exl + wxh*(Exh-Exl);

      // Calculate Eyi, the linearly-weighted Ey at the particle: 
      int iyl_m1 = iyl-1;
      if (iyl_m1 < 0) iyl_m1 = ny-1;
      int iyl_p1=iyh;
      int iyh_m1 = iyl;
      int iyh_p1=iyh+1;
      if (iyh_p1 == ny) iyh_p1 = 0;

      Eyll = (Efield.phi_rho(ixl,iyl_p1) - Efield.phi_rho(ixl,iyl_m1)+Ey0dy)*Eyscale;
      Eylh = (Efield.phi_rho(ixl,iyh_p1) - Efield.phi_rho(ixl,iyh_m1)+Ey0dy)*Eyscale;
      Eyhl = (Efield.phi_rho(ixh,iyl_p1) - Efield.phi_rho(ixh,iyl_m1)+Ey0dy)*Eyscale;
      Eyhh = (Efield.phi_rho(ixh,iyh_p1) - Efield.phi_rho(ixh,iyh_m1)+Ey0dy)*Eyscale;

      Eyl = Eyll + wyh*(Eylh-Eyll);
      Eyh = Eyhl + wyh*(Eyhh-Eyhl);
      Eyi = Eyl + wxh*(Eyh-Eyl);

      // Apply the Boris mover (Ezi is assumed 0 in 2D). 
      vya = old.vy(i) + Eyi;
      vxa = old.vx(i) + Exi;
      vxb = vxa + vya * bt;
      vyc = vya - vxb * bs;
      vxc = vxb + vyc * bt;
      
      // Advance the particle velocities: 
      adv.vx(i) = vxc + Exi;
      adv.vy(i) = vyc + Eyi;

      adv.vz(i) = old.vz(i);

    } //end if (cur.x(i) >= nxmin) 

  } // End for 

  
  return ;

    
} // vpush 
#endif

#if NDIM == 3
/*
void interpolate_efield(FTYPE &Ex, FTYPE &Ey, FTYPE &Ez,
                        FTYPE x, FTYPE y, FTYPE z, 
                        FArrayND_ranged &phi_rho);
*/

void vpushBz_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
                    FTYPE qmrat, FTYPE alpha, FTYPE Bz)
{
  
  const PTYPE nxmin=0; // This is the minimum acceptable value for particle positions

  FTYPE Exscale;
  FTYPE Eyscale;
  FTYPE Ezscale;
  int i;
  int np = old.np;
  FTYPE Exi, Eyi, Ezi;
  
  // A vpush routine for non-zero Bz 
  FTYPE wcidt, bt, bs;
  FTYPE vxa, vya, vxb, vxc, vyc;
    
  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop.  
  Exscale =  0.5 * alpha * qmrat * (dt*dt/dx) / (-2*dx);

  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Eyscale =  0.5 * alpha * qmrat * (dt*dt/dy) / (-2*dy);

  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Ezscale =  alpha * qmrat * (dt*dt/dz) / (-2*dz);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2*dx);
  FTYPE Ey0dy=Ey0_external*(-2*dy);

  // Rotation information - Boris constants 
  wcidt = alpha * dt/2 * (qmrat*Bz);
  bt = wcidt;
  bs = 2*bt/(1+bt*bt);
  // Scale constants to put velocity variables into appropriately normalized system
  bt *= dy/dx;
  bs *= dx/dy;

  // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {
    if (cur.x(i) >= nxmin){
      // Get the E field at the particle
      interpolate_efield(Exi, Eyi, Ezi, 
                         cur.x(i), cur.y(i), cur.z(i), 
                         Efield.phi_rho);
      
      // Do appropriate scaling and add external fields
      Exi += Ex0dx;
      Exi *= Exscale;
      Eyi += Ey0dy;
      Eyi *= Eyscale;
      Ezi *= Ezscale;

      // Apply the Boris mover 
      vxa = adv.vx(i) + Exi;
      vya = adv.vy(i) + Eyi;
      vxb = vxa + vya * bt;
      vyc = vya - vxb * bs;
      vxc = vxb + vyc * bt;
      
      // Advance the particle velocities: 
      adv.vx(i) = vxc + Exi;
      adv.vy(i) = vyc + Eyi;
      adv.vz(i) += Ezi;
    } //end if (cur.x(i) >= nxmin) 

  } // End for 

  //  Restore original values to the electric field array: 
  
  
  return ;

    
} // vpush3dBz


#endif

