/* Advance the 1D unmagnatized particle velocities to the next time step. */

#include "eppic.h"
void vpush1D(particle_dist &adv, particle_dist &old, particle_dist &cur, field &Efield, 
	   FTYPE qmrat, FTYPE alpha)
{
  
  FTYPE Exscale;
  int i;
  int np = old.x.length();
  int ixl, ixh;
  int iyl, iyh;
  FTYPE wxh;
  FTYPE wyh;
  
  void gradient(FArrayND &, INDICIES(FArrayND &, FArrayND &, FArrayND &), 
		INDICIES(FTYPE, FTYPE, FTYPE), FTYPE scaler);

  /* Rescale the electric field so that it may simply be added to
	particle_dist velocities in the particle loop. */
  Exscale =  alpha * qmrat * (dt*dt/dx);  
  FArrayND Ex_scaled(INDICIES(nx, ny, nz));
  if (Efield.Ex.length()==Ex_scaled.length()) Ex_scaled=Efield.Ex;
  else
    gradient(Efield.phi_rho, INDICIES(Ex_scaled, Efield.Ey, Efield.Ez), 
	     INDICIES(dx, dy, dz), (FTYPE) -1.); 
  Ex_scaled *= Exscale;


     /* For each particle, update the velocity ... */
  for (i = 0; i < np; ++i) {

    /* Define the nearest grid points: */
    /* And the corresponding linear weighting factors: */
    ixl = (int) cur.x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;
    wxh = cur.x(i) - ixl;

#if NDIM > 1    
    // 2-D fields
    iyl = (int) cur.y(i);
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;
    wyh = cur.y(i) - iyl;
#endif

#if NDIM == 2
    /* Calculate Exi, the linearly-weighted Ex at the particle: */
    FTYPE Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;
    Exll = Ex_scaled(ixl,iyl);
    Exlh = Ex_scaled(ixl,iyh);
    Exhl = Ex_scaled(ixh,iyl);
    Exhh = Ex_scaled(ixh,iyh);
    Exl = Exll + wyh*(Exlh-Exll);
    Exh = Exhl + wyh*(Exhh-Exhl);
    Exi = Exl + wxh*(Exh-Exl);

#elif NDIM == 3
    // 3-D fields
    int izl,izh;
    FTYPE wzh;
    FTYPE Elll,Ellh,Elhl,Elhh,Ehll,Ehlh,Ehhl,Ehhh;
    FTYPE Ell, Ehl, Elh, Ehh, Eil, Eih, Exi;

    izl = (int) cur.z(i);
    izh = izl + 1;
    if (izh == nz) izh = 0;
    wzh = cur.z(i) - izl;

    /* Calculate Exi, the linearly-weighted Ex at the particle: */
    Elll = Ex_scaled(ixl,iyl,izl);
    Ellh = Ex_scaled(ixl,iyl,izh);
    Elhl = Ex_scaled(ixl,iyh,izl);
    Elhh = Ex_scaled(ixl,iyh,izh);
    Ehll = Ex_scaled(ixh,iyl,izl);
    Ehlh = Ex_scaled(ixh,iyl,izh);
    Ehhl = Ex_scaled(ixh,iyh,izl);
    Ehhh = Ex_scaled(ixh,iyh,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
      
    Eil = Ell + wxh*(Ehl-Ell);
    Eih = Elh + wxh*(Ehh-Elh);
    Exi = Eil + wzh*(Eih-Eil);
    /* Advance the particle velocities: */

#else
    // 1-D fields:
    /* Calculate Exi, the linearly-weighted Ex at the particle: */
    FTYPE Exl, Exh, Exi;
    Exl = Ex_scaled(ixl);
    Exh = Ex_scaled(ixh);
    Exi = Exl + wxh*(Exh-Exl);

#endif

    /* Advance the particle velocities: */
    adv.vx(i) = old.vx(i) + Exi;

  } /* End for */
  
  return ;

    
} /* vpush */
