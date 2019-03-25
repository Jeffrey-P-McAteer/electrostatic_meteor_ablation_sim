/* 3 Routines below for 

   Advancing the unmagnatized particle velocities to the next time step
   in 1-D, 2-D and 3-D.

*/
#include "eppic.h"
#include "interpolate_efield.h"

#if NDIM == 1
void vpush_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	   FTYPE qmrat, FTYPE alpha)
{
  
   
  const PTYPE nxmin=0; // This is the minimum acceptable value for particle positions
  FTYPE Exscale, Exiscale;
  int i;
  int np = old.x.size();
  int ixl, ixh;
  FTYPE wxh;
  FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;

  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop. 
  Exscale =  alpha * qmrat * (dt*dt/dx)/(-2*dx);

  // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {
    
    if (cur.x(i) >= nxmin) { // This excludes absent particles which should always be <0

    // Define the nearest grid points: 
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    //    if (ixh == nx) ixh = 0; // not periodic

    // And the corresponding linear weighting factors: 
    wxh = cur.x(i) - ixl;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    int ixl_m1 = ixl-1;
    int ixl_p1=ixh;
    int ixh_m1 = ixl;
    int ixh_p1=ixh+1;

    Exl = (Efield.phi_rho(ixl_p1) - Efield.phi_rho(ixl_m1))*Exscale;
    Exh = (Efield.phi_rho(ixh_p1) - Efield.phi_rho(ixh_m1))*Exscale;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    Exi = Exl + wxh*(Exh-Exl);

    // Advance the particle velocities: 
    adv.vx(i) = old.vx(i) + Exi;
    } // end absent particle test
  } // End for 

  //  Restore original values to the electric field array: 

 
  return ;

    
} // End 1-D vpush 

#endif
#if NDIM == 2
void vpush_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	   FTYPE qmrat, FTYPE alpha)
{
  
  const PTYPE nxmin=0; // This is the minimum acceptable value for particle positions
  FTYPE Exscale;
  FTYPE Eyscale;
  int i;
  int np = old.x.length();
  int ixl, ixh;
  int iyl, iyh;
  FTYPE wxh;
  FTYPE wyh;
  FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;
  FTYPE Eyll, Eylh, Eyhl, Eyhh, Eyl, Eyh, Eyi;

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2*dx);
  FTYPE Ey0dy=Ey0_external*(-2*dy);

  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop. 
  Exscale =  alpha * qmrat * (dt*dt/dx)/(-2*dx) ;
  Eyscale =  alpha * qmrat * (dt*dt/dy)/(-2*dy) ;

     // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {

    if (cur.x(i) >= nxmin) { // This excludes absent particles which should always be <0
    // Define the nearest grid points: 
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    //    if (ixh == nx) ixh = 0;

    iyl = (int) cur.y(i);
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;

    // And the corresponding linear weighting factors: 
    wxh = cur.x(i) - ixl;
    wyh = cur.y(i) - iyl;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    int ixl_m1 = ixl-1;
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

    // Advance the particle velocities: 
    adv.vx(i) = old.vx(i) + Exi;
    adv.vy(i) = old.vy(i) + Eyi;

    } // end asbsent particle test
  } // End for 

  //  Restore original values to the electric field array: 


 
  return ;

    
} // End 2-D vpush 

#endif

#if NDIM == 3
/*void interpolate_efield(FTYPE &Ex, FTYPE &Ey, FTYPE &Ez,
                        FTYPE x, FTYPE y, FTYPE z, 
                        FArrayND_ranged &phi_rho);
*/

void vpush_domain(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     FTYPE qmrat, FTYPE alpha)
{

  
  const PTYPE nxmin=0; // This is the minimum acceptable value for particle positions
  FTYPE Exscale, Eyscale, Ezscale;

  const int np = old.x.length();

  FTYPE Exi, Eyi, Ezi;


  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop. 
  // dv_x = q*E_x*dt/m * (dt/dx) * (1/(-2*dx))
  // (dt/dx) to normalize velocity to grid/timestep units
  // (1/(-2*dx)) to complete interpolation problem in interpolate_efield
  Exscale =  alpha * qmrat * (dt*dt/dx) / (-2*dx);
  Eyscale =  alpha * qmrat * (dt*dt/dy) / (-2*dy);
  Ezscale =  alpha * qmrat * (dt*dt/dz) / (-2*dz);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2.*dx);
  FTYPE Ey0dy=Ey0_external*(-2.*dy);

  // For each particle, update the velocity ... 
  for (int i = 0; i < np; ++i) {
    
    if (cur.x(i) >= nxmin) { // This excludes absent particles which should always be <0
      interpolate_efield(Exi, Eyi, Ezi,
                         cur.x(i), cur.y(i), cur.z(i),
                         Efield.phi_rho);
      
      adv.vx(i) = old.vx(i) + (Exi+Ex0dx)*Exscale;
      adv.vy(i) = old.vy(i) + (Eyi+Ey0dy)*Eyscale;
      //adv.vz(i) = old.vz(i) + Ezi;
      adv.vz(i) = old.vz(i) + Ezi*Ezscale;
    } // end absent particle test
  } // End for 

 
  return ;

    
} // vpush

#endif
