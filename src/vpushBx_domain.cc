/* 3 routines below for 1-D 2-D and 3-D to
   
Advance the magnatized particle velocities to the next time step. 

Assumes periodic in y and 
*/
#include <math.h>
#include "eppic.h"
#include "interpolate_efield.h"

#if NDIM == 1
void vpushBx_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     FTYPE qmrat, FTYPE alpha, FTYPE Bx)
{
  

  // Just use the unmagnetized routine:
  // code not written yet:
  void vpush_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, 
		    field &Efield,  FTYPE qmrat, FTYPE alpha);

  vpush_domain(adv, old, cur, Efield, qmrat, alpha);
  
}
#endif


#if NDIM == 2
void vpushBx_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     FTYPE qmrat, FTYPE alpha, FTYPE Bx)
{
  
  const PTYPE nxmin=-1.0; // This is the minimum acceptable value for particle positions

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

  // A vpush routine for non-zero Bx 
  FTYPE wcidt, bt, bs;
  FTYPE vya, vza, vzb, vyc, vzc;
    
  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop.  
  Exscale =  alpha * qmrat * (dt*dt/dx) / (-2*dx) ;

  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Eyscale =  0.5 * alpha * qmrat * (dt*dt/dy) / (-2*dy);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2*dx);
  FTYPE Ey0dy=Ey0_external*(-2*dy);

  // Rotation information - Boris constants 
  wcidt = alpha * dt/2 * (qmrat*Bx);
  bt = wcidt;
  bs = 2*bt/(1+bt*bt);
      
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
    Exll = (Efield.phi_rho(ixl_p1,iyl) - Efield.phi_rho(ixl_m1,iyl)+Ex0dx)*Exscale;
    Exlh = (Efield.phi_rho(ixl_p1,iyh) - Efield.phi_rho(ixl_m1,iyh)+Ex0dx)*Exscale;
    Exhl = (Efield.phi_rho(ixh_p1,iyl) - Efield.phi_rho(ixh_m1,iyl)+Ex0dx)*Exscale;
    Exhh = (Efield.phi_rho(ixh_p1,iyh) - Efield.phi_rho(ixh_m1,iyh)+Ex0dx)*Exscale;

    Exl = Exll + wyh*(Exlh-Exll);
    Exh = Exhl + wyh*(Exhh-Exhl);
    Exi = Exl + wxh*(Exh-Exl);

    // Advance the parallel velocity here
    adv.vx(i) = old.vx(i) + Exi;

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
    vza = old.vz(i);
    vzb = vza - vya * bt;
    vyc = vya + vzb * bs;
    vzc = vzb - vyc * bt;
      
    // Advance the particle velocities: 
    adv.vy(i) = vyc + Eyi;
    adv.vz(i) = vzc;

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

void vpushBx_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     FTYPE qmrat, FTYPE alpha, FTYPE Bx)
{
  FTYPE Exscale;
  FTYPE Eyscale;
  FTYPE Ezscale;
  int i;
  int np = old.np;
  int ixl, ixh;
  int iyl, iyh;
  int izl, izh;
  int izl2, izh2; // For one-sided differences in z
  FTYPE wxh, wyh, wzh;
  FTYPE Elll2, Ellh2, Elhl2, Elhh2, Ehhl2, Ehhh2; // For one-sided differences in z
  FTYPE Elll, Ellh, Elhl, Elhh, Ehll, Ehlh, Ehhl, Ehhh;
  FTYPE Ell, Elh, Ehl, Ehh, El, Eh;
  FTYPE Exi, Eyi, Ezi;

  // A vpush routine for non-zero Bx
  FTYPE wcidt, bt, bs;
  FTYPE vya, vza, vzb, vyc, vzc;
    
  const PTYPE nxmin=-1.0; // This is the minimum acceptable value for particle positions


  Exscale =  alpha * qmrat * (dt*dt/dx) / (-2*dx);
  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Eyscale =  0.5 * alpha * qmrat * (dt*dt/dy) / (-2*dy);
  // The extra 0.5 comes from the "half"-accel in the Boris mover.) 
  Ezscale =  0.5 * alpha * qmrat * (dt*dt/dz) / (-2*dz);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2*dx);
  FTYPE Ey0dy=Ey0_external*(-2*dy);

  // Rotation information - Boris constants 
  wcidt = alpha * dt/2 * (qmrat*Bx);
  bt = wcidt;
  bs = 2*bt/(1+bt*bt);
  // Scale constants to put velocity variables into appropriately normalized system
  bt *= dy/dz;
  bs *= dz/dy;

  // finite_diff method = -2 (backwards for zh and zhp1) -1 (backwards zhp1 only), 
  //                      0 (central), 1 (forwards for zlm1 only)
  // Variables added by Glenn on 4/17/2018 to deal with non-periodic z boundary
  int finite_diff_method;
  // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {

    if (cur.x(i) >= nxmin) { // This excludes absent particles which should always be <0
      // Get the E field at the particle
      interpolate_efield(Exi, Eyi, Ezi, 
                         cur.x(i), cur.y(i), cur.z(i), 
                         Efield.phi_rho);

      Exi += Ex0dx;
      Eyi += Ey0dy;
      Exi *= Exscale;
      Eyi *= Eyscale;
      Ezi *= Ezscale;

      
      // Apply the Boris mover (Ezi is assumed 0 in 2D). 
      vya = adv.vy(i) + Eyi;
      vza = adv.vz(i) + Ezi;
      vzb = vza - vya * bt;
      vyc = vya + vzb * bs;
      vzc = vzb - vyc * bt;
      
      // Advance the particle velocities: 
      adv.vx(i) += Exi;
      adv.vy(i) = vyc + Eyi;
      adv.vz(i) = vzc + Ezi;
    } // end if (cur.x(i) >= nxmin
  } // end for loop over particles
  
  return ;

    
} // vpush3dBx


#endif

