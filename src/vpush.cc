/* 3 Routines below for 

   Advancing the unmagnatized particle velocities to the next time step
   in 1-D, 2-D and 3-D.

*/
#include "eppic.h"

#if NDIM == 1
void vpush(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	   FTYPE qmrat, FTYPE alpha)
{
  
  FTYPE Exscale, Exiscale;
  int i;
  int np = old.x.size();
  int ixl, ixh;
  FTYPE wxh;
  FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;

  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop. 
  Exscale =  alpha * qmrat * (dt*dt/dx)/(-2*dx) ;

  // For each particle, update the velocity ... 
  for (i = 0; i < np; ++i) {

    // Define the nearest grid points: 
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;

    // And the corresponding linear weighting factors: 
    wxh = cur.x(i) - ixl;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    int ixl_m1 = ixl-1;
    if (ixl_m1 < 0) ixl_m1 = nx-1;
    int ixl_p1=ixh;
    int ixh_m1 = ixl;
    int ixh_p1=ixh+1;
    if (ixh_p1 == nx) ixh_p1 = 0;

    Exl = (Efield.phi_rho(ixl_p1) - Efield.phi_rho(ixl_m1))*Exscale;
    Exh = (Efield.phi_rho(ixh_p1) - Efield.phi_rho(ixh_m1))*Exscale;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    Exi = Exl + wxh*(Exh-Exl);

    // Advance the particle velocities: 
    adv.vx(i) = old.vx(i) + Exi;

  } // End for 

  //  Restore original values to the electric field array: 

 
  return ;

    
} // End 1-D vpush 

#endif
#if NDIM == 2
void vpush(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	   FTYPE qmrat, FTYPE alpha)
{
  
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

    // Define the nearest grid points: 
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;

    iyl = (int) cur.y(i);
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;

    // And the corresponding linear weighting factors: 
    wxh = cur.x(i) - ixl;
    wyh = cur.y(i) - iyl;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    int ixl_m1 = ixl-1;
    if (ixl_m1 < 0) ixl_m1 = nx-1;
    int ixl_p1=ixh;
    int ixh_m1 = ixl;
    int ixh_p1=ixh+1;
    if (ixh_p1 == nx) ixh_p1 = 0;

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

  } // End for 

  //  Restore original values to the electric field array: 

 
  return ;

    
} // End 2-D vpush 

#endif

#if NDIM == 3
void vpush(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	     FTYPE qmrat, FTYPE alpha)
{
  
  FTYPE Exscale, Eyscale, Ezscale;

  const int np = old.x.length();
  int ixl, ixh;
  int iyl, iyh;
  int izl, izh;
  FTYPE wxh, wyh, wzh;
  FTYPE Elll, Ellh, Elhl, Elhh, Ehll, Ehlh, Ehhl, Ehhh;
  FTYPE Ell, Elh, Ehl, Ehh, El, Eh;
  FTYPE Exi, Eyi, Ezi;


  // Rescale the electric field so that it may simply be added to
  // particle velocities in the particle loop. 
  Exscale =  alpha * qmrat * (dt*dt/dx) / (-2*dx);
  Eyscale =  alpha * qmrat * (dt*dt/dy) / (-2*dy);
  Ezscale =  alpha * qmrat * (dt*dt/dz) / (-2*dz);

  // An external, Ey0_external, driver needs scaling for rapid use
  FTYPE Ex0dx=Ex0_external*(-2.*dx);
  FTYPE Ey0dy=Ey0_external*(-2.*dy);

  // For each particle, update the velocity ... 
  for (int i = 0; i < np; ++i) {

    // Define the nearest grid points: 
    // And the corresponding linear weighting factors: 
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;
    wxh = cur.x(i) - ixl;

    iyl = (int) cur.y(i);
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;
    wyh = cur.y(i) - iyl;

    izl = (int) cur.z(i);
    izh = izl + 1;
    if (izh == nz) izh = 0;
    wzh = cur.z(i) - izl;

    // Calculate Exi, the linearly-weighted Ex at the particle: 
    int ixl_m1 = ixl-1;
    if (ixl_m1 < 0) ixl_m1 = nx-1;
    int ixl_p1=ixh;
    int ixh_m1 = ixl;
    int ixh_p1=ixh+1;
    if (ixh_p1 == nx) ixh_p1 = 0;
    
    Elll = (Efield.phi_rho(ixl_p1,iyl,izl) - Efield.phi_rho(ixl_m1,iyl,izl)+Ex0dx)*Exscale;
    Ellh = (Efield.phi_rho(ixl_p1,iyl,izh) - Efield.phi_rho(ixl_m1,iyl,izh)+Ex0dx)*Exscale;
    Elhl = (Efield.phi_rho(ixl_p1,iyh,izl) - Efield.phi_rho(ixl_m1,iyh,izl)+Ex0dx)*Exscale;
    Elhh = (Efield.phi_rho(ixl_p1,iyh,izh) - Efield.phi_rho(ixl_m1,iyh,izh)+Ex0dx)*Exscale;
    Ehll = (Efield.phi_rho(ixh_p1,iyl,izl) - Efield.phi_rho(ixh_m1,iyl,izl)+Ex0dx)*Exscale;
    Ehlh = (Efield.phi_rho(ixh_p1,iyl,izh) - Efield.phi_rho(ixh_m1,iyl,izh)+Ex0dx)*Exscale;
    Ehhl = (Efield.phi_rho(ixh_p1,iyh,izl) - Efield.phi_rho(ixh_m1,iyh,izl)+Ex0dx)*Exscale;
    Ehhh = (Efield.phi_rho(ixh_p1,iyh,izh) - Efield.phi_rho(ixh_m1,iyh,izh)+Ex0dx)*Exscale;

    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
      
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Exi = El + wzh*(Eh-El);
    // Advance the particle velocities: 
    adv.vx(i) = old.vx(i) + Exi;

    // Calculate Eyi, the linearly-weighted Ey at the particle: 
    int iyl_m1 = iyl-1;
    if (iyl_m1 < 0) iyl_m1 = ny-1;
    int iyl_p1=iyh;
    int iyh_m1 = iyl;
    int iyh_p1=iyh+1;
    if (iyh_p1 == ny) iyh_p1 = 0;

    Elll = (Efield.phi_rho(ixl,iyl_p1,izl) - Efield.phi_rho(ixl,iyl_m1,izl)+Ey0dy)*Eyscale;
    Ellh = (Efield.phi_rho(ixl,iyl_p1,izh) - Efield.phi_rho(ixl,iyl_m1,izh)+Ey0dy)*Eyscale;
    Elhl = (Efield.phi_rho(ixl,iyh_p1,izl) - Efield.phi_rho(ixl,iyh_m1,izl)+Ey0dy)*Eyscale;
    Elhh = (Efield.phi_rho(ixl,iyh_p1,izh) - Efield.phi_rho(ixl,iyh_m1,izh)+Ey0dy)*Eyscale;
    Ehll = (Efield.phi_rho(ixh,iyl_p1,izl) - Efield.phi_rho(ixh,iyl_m1,izl)+Ey0dy)*Eyscale;
    Ehlh = (Efield.phi_rho(ixh,iyl_p1,izh) - Efield.phi_rho(ixh,iyl_m1,izh)+Ey0dy)*Eyscale;
    Ehhl = (Efield.phi_rho(ixh,iyh_p1,izl) - Efield.phi_rho(ixh,iyh_m1,izl)+Ey0dy)*Eyscale;
    Ehhh = (Efield.phi_rho(ixh,iyh_p1,izh) - Efield.phi_rho(ixh,iyh_m1,izh)+Ey0dy)*Eyscale;

    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
      
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Eyi = El + wzh*(Eh-El);
    // Advance the particle velocities: 
    adv.vy(i) = old.vy(i) + Eyi;

    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    int izl_m1 = izl-1;
    if (izl_m1 < 0) izl_m1 = nz-1;
    int izl_p1=izh;
    int izh_m1 = izl;
    int izh_p1=izh+1;
    if (izh_p1 == nz) izh_p1 = 0;

    Elll = (Efield.phi_rho(ixl,iyl,izl_p1) - Efield.phi_rho(ixl,iyl,izl_m1))*Ezscale;
    Ellh = (Efield.phi_rho(ixl,iyl,izh_p1) - Efield.phi_rho(ixl,iyl,izh_m1))*Ezscale;
    Elhl = (Efield.phi_rho(ixl,iyh,izl_p1) - Efield.phi_rho(ixl,iyh,izl_m1))*Ezscale;
    Elhh = (Efield.phi_rho(ixl,iyh,izh_p1) - Efield.phi_rho(ixl,iyh,izh_m1))*Ezscale;
    Ehll = (Efield.phi_rho(ixh,iyl,izl_p1) - Efield.phi_rho(ixh,iyl,izl_m1))*Ezscale;
    Ehlh = (Efield.phi_rho(ixh,iyl,izh_p1) - Efield.phi_rho(ixh,iyl,izh_m1))*Ezscale;
    Ehhl = (Efield.phi_rho(ixh,iyh,izl_p1) - Efield.phi_rho(ixh,iyh,izl_m1))*Ezscale;
    Ehhh = (Efield.phi_rho(ixh,iyh,izh_p1) - Efield.phi_rho(ixh,iyh,izh_m1))*Ezscale;

    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
      
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ezi = El + wzh*(Eh-El);
    // Advance the particle velocities: 
    adv.vz(i) = old.vz(i) + Ezi;

  } // End for 

 
  return ;

    
} // vpush

#endif
