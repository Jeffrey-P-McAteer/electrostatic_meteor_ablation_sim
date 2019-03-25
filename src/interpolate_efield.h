// This will calculate the E field and interpolate the results to a paricle's
// location.
// Glenn Sugar 20180419
#ifndef INTERPOLATE_EFIELD_H
#define INTERPOLATE_EFIELD_H

#include <math.h>
#include "eppic.h"

#if NDIM == 1
FTYPE interpolate_efield(FTYPE x, FArrayND_ranged &phi_rho);
#endif

#if NDIM == 2
FTYPE* interpolate_efield(FTYPE x, FTYPE y, FArrayND_ranged &phi_rho);
#endif

#if NDIM == 3
/*void interpolate_efield(FTYPE &Ex, FTYPE &Ey, FTYPE &Ez,
                        FTYPE x, FTYPE y, FTYPE z, 
                        FArrayND_ranged &phi_rho);
*/
inline void interpolate_efield(FTYPE &Ex, FTYPE &Ey, FTYPE &Ez,
                               FTYPE x, FTYPE y, FTYPE z, 
                               FArrayND_ranged &phi_rho)
{
  int ixl, ixh, iyl, iyh, izl, izh;
  FTYPE wxh, wyh, wzh;
  FTYPE Elll, Ellh, Elhl, Elhh, Ehll, Ehlh, Ehhl, Ehhh;
  FTYPE Ell, Elh, Ehl, Ehh, El, Eh;

  // Define the nearest grid points: 
  // And the corresponding linear weighting factors: 
  ixl = (int) x;
  ixh = ixl + 1;
  wxh = x - ixl;
  
  iyl = (int) y;
  iyh = iyl + 1;
  if (iyh == ny+yguard_size) iyh = 0;
  wyh = y - iyl;
  
  izl = (int) z;
  izh = izl + 1;
  if (izh == nz+zguard_size) izh = 0;
  wzh = z - izl;

  // Calculate Exi, the linearly-weighted Ex at the particle:
  // Will be central difference if periodic field BC, middle subdomain,
  // or middle grid cells
  if (field_boundary_type[0][0] == periodic || 
      (subdomain.id_number > 0 && subdomain.id_number < nsubdomains-1) ||
      ((ixh < nx) && (ixl > 0))){
    int ixl_m1 = ixl-1;
    int ixl_p1 = ixh;
    int ixh_m1 = ixl;
    int ixh_p1 = ixh+1;
    Elll = phi_rho(ixl_p1,iyl,izl) - phi_rho(ixl_m1,iyl,izl);
    Ellh = phi_rho(ixl_p1,iyl,izh) - phi_rho(ixl_m1,iyl,izh);
    Elhl = phi_rho(ixl_p1,iyh,izl) - phi_rho(ixl_m1,iyh,izl);
    Elhh = phi_rho(ixl_p1,iyh,izh) - phi_rho(ixl_m1,iyh,izh);
    Ehll = phi_rho(ixh_p1,iyl,izl) - phi_rho(ixh_m1,iyl,izl);
    Ehlh = phi_rho(ixh_p1,iyl,izh) - phi_rho(ixh_m1,iyl,izh);
    Ehhl = phi_rho(ixh_p1,iyh,izl) - phi_rho(ixh_m1,iyh,izl);
    Ehhh = phi_rho(ixh_p1,iyh,izh) - phi_rho(ixh_m1,iyh,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ex = El + wzh*(Eh-El);
  }
  else if (ixh >= nx){
    // Nonperiodic boundary at rhs boundary.  Use backward difference
    // Get the needed indices for interpolation
    int ixl_m1 = ixl-1;
    int ixl_p1 = ixh;
    int ixh_m1 = ixl;
    int ixh_m2 = ixl_m1;

    // Backwards difference the x=h variables since we don't know rho(nx)
    // Calculate Exi, the linearly-weighted Ex at the particle: 
    Elll = phi_rho(ixl_p1,iyl,izl) - phi_rho(ixl_m1,iyl,izl);
    Ellh = phi_rho(ixl_p1,iyl,izh) - phi_rho(ixl_m1,iyl,izh);
    Elhl = phi_rho(ixl_p1,iyh,izl) - phi_rho(ixl_m1,iyh,izl);
    Elhh = phi_rho(ixl_p1,iyh,izh) - phi_rho(ixl_m1,iyh,izh);
    Ehll = 3*phi_rho(ixh,iyl,izl) - 4*phi_rho(ixh_m1,iyl,izl)
      + phi_rho(ixh_m2,iyl,izl);
    Ehlh = 3*phi_rho(ixh,iyl,izh) - 4*phi_rho(ixh_m1,iyl,izh)
      + phi_rho(ixh_m2,iyl,izh);
    Ehhl = 3*phi_rho(ixh,iyh,izl) - 4*phi_rho(ixh_m1,iyh,izl)
      + phi_rho(ixh_m2,iyh,izl);
    Ehhh = 3*phi_rho(ixh,iyh,izh) - 4*phi_rho(ixh_m1,iyh,izh)
      + phi_rho(ixh_m2,iyh,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ey = El + wzh*(Eh-El);       
  }
  else {
    // Forward Difference in x since we don't know ixl_m1
    int ixh_m1 = ixl;
    int ixh_p1 = ixh+1;
    int ixl_p1 = ixh;
    int ixl_p2 = ixh_p1;

    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = -3*phi_rho(ixl,iyl,izl) + 4*phi_rho(ixl_p1,iyl,izl)
      -phi_rho(ixl_p2,iyl,izl);
    Ellh = -3*phi_rho(ixl,iyl,izh) + 4*phi_rho(ixl_p1,iyl,izh)
      -phi_rho(ixl_p2,iyl,izh);
    Elhl = -3*phi_rho(ixl,iyh,izl) + 4*phi_rho(ixl_p1,iyh,izl)
      -phi_rho(ixl_p2,iyh,izl);
    Elhh = -3*phi_rho(ixl,iyh,izh) + 4*phi_rho(ixl_p1,iyh,izh)
      -phi_rho(ixl_p2,iyh,izh);
    Ehll = phi_rho(ixh_p1,iyl,izl) - phi_rho(ixh_m1,iyl,izl);
    Ehlh = phi_rho(ixh_p1,iyl,izh) - phi_rho(ixh_m1,iyl,izh);
    Ehhl = phi_rho(ixh_p1,iyh,izl) - phi_rho(ixh_m1,iyh,izl);
    Ehhh = phi_rho(ixh_p1,iyh,izh) - phi_rho(ixh_m1,iyh,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ez = El + wzh*(Eh-El);
  }

  // Calculate Eyi, the linearly-weighted Ey at the particle:
  // Ey could be central, forward, or backward difference
  if (field_boundary_type[1][0] == periodic || ((iyh < ny+yguard_size) &&
                                                (iyl > 0))){
    // Periodic boundary or internal node.  Use central difference
    int iyl_m1 = iyl-1;
    if (iyl_m1 < 0) iyl_m1 = ny+yguard_size-1;
    int iyl_p1 = iyh;
    int iyh_m1 = iyl;
    int iyh_p1 = iyh+1;
    if (iyh_p1 == ny+yguard_size) iyh_p1 = 0;

    // Calculate Eyi, the linearly-weighted Ey at the particle: 
    Elll = phi_rho(ixl,iyl_p1,izl) - phi_rho(ixl,iyl_m1,izl);
    Ellh = phi_rho(ixl,iyl_p1,izh) - phi_rho(ixl,iyl_m1,izh);
    Elhl = phi_rho(ixl,iyh_p1,izl) - phi_rho(ixl,iyh_m1,izl);
    Elhh = phi_rho(ixl,iyh_p1,izh) - phi_rho(ixl,iyh_m1,izh);
    Ehll = phi_rho(ixh,iyl_p1,izl) - phi_rho(ixh,iyl_m1,izl);
    Ehlh = phi_rho(ixh,iyl_p1,izh) - phi_rho(ixh,iyl_m1,izh);
    Ehhl = phi_rho(ixh,iyh_p1,izl) - phi_rho(ixh,iyh_m1,izl);
    Ehhh = phi_rho(ixh,iyh_p1,izh) - phi_rho(ixh,iyh_m1,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ey = El + wzh*(Eh-El);
  }
  else if (iyh >= ny+yguard_size){
    // Nonperiodic boundary at rhs boundary.  Use backward difference
    // Get the needed indices for interpolation
    int iyl_m1 = iyl-1;
    int iyl_p1 = iyh;
    int iyh_m1 = iyl;
    int iyh_m2 = iyl_m1;

    // Backwards difference the z=h variables since we don't know rho(nz)
    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = phi_rho(ixl,iyl_p1,izl) - phi_rho(ixl,iyl_m1,izl);
    Ellh = phi_rho(ixl,iyl_p1,izh) - phi_rho(ixl,iyl_m1,izh);
    Elhl = 3*phi_rho(ixl,iyh,izl) - 4*phi_rho(ixl,iyh_m1,izl)
      + phi_rho(ixl,iyh_m2,izl);
    Elhh = 3*phi_rho(ixl,iyh,izh) - 4*phi_rho(ixl,iyh_m1,izh)
      + phi_rho(ixl,iyh_m2,izh);
    Ehll = phi_rho(ixh,iyl_p1,izl) - phi_rho(ixh,iyl_m1,izl);
    Ehlh = phi_rho(ixh,iyl_p1,izh) - phi_rho(ixh,iyl_m1,izh);
    Ehhl = 3*phi_rho(ixh,iyh,izl) - 4*phi_rho(ixh,iyh_m1,izl)
      + phi_rho(ixh,iyh_m2,izl);
    Ehhh = 3*phi_rho(ixh,iyh,izh) - 4*phi_rho(ixh,iyh_m1,izh)
      + phi_rho(ixh,iyh_m2,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ey = El + wzh*(Eh-El);       
  }
  else{
    // Nonperiodic boundary, at lhs boundary.  Use forward difference
    int iyh_m1 = iyl;
    int iyh_p1 = iyh+1;
    int iyl_p1 = iyh;
    int iyl_p2 = iyh_p1;

    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = -3*phi_rho(ixl,iyl,izl) + 4*phi_rho(ixl,iyl_p1,izl)
      -phi_rho(ixl,iyl_p2,izl);
    Ellh = -3*phi_rho(ixl,iyl,izh) + 4*phi_rho(ixl,iyl_p1,izh)
      -phi_rho(ixl,iyl_p2,izh);
    Elhl = phi_rho(ixl,iyh_p1,izl) - phi_rho(ixl,iyh_m1,izl);
    Elhh = phi_rho(ixl,iyh_p1,izh) - phi_rho(ixl,iyh_m1,izh);
    Ehll = -3*phi_rho(ixh,iyl,izl) + 4*phi_rho(ixh,iyl_p1,izl)
      -phi_rho(ixh,iyl_p2,izl);
    Ehlh = -3*phi_rho(ixh,iyl,izh) + 4*phi_rho(ixh,iyl_p1,izh)
      -phi_rho(ixh,iyl_p2,izh);
    Ehhl = phi_rho(ixh,iyh_p1,izl) - phi_rho(ixh,iyh_m1,izl);
    Ehhh = phi_rho(ixh,iyh_p1,izh) - phi_rho(ixh,iyh_m1,izh);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ey = El + wzh*(Eh-El);
  }

  if (field_boundary_type[2][0] == periodic || ((izh < nz+zguard_size) &&
                                                (izl > 0))){
    // Periodic boundary or internal node.  Use central difference
    int izl_m1 = izl-1;
    if (izl_m1 < 0) izl_m1 = nz-1;
    int izl_p1=izh;
    int izh_m1 = izl;
    int izh_p1=izh+1;
    if (izh_p1 == nz+zguard_size) izh_p1 = 0;
    
    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = phi_rho(ixl,iyl,izl_p1) - phi_rho(ixl,iyl,izl_m1);
    Ellh = phi_rho(ixl,iyl,izh_p1) - phi_rho(ixl,iyl,izh_m1);
    Elhl = phi_rho(ixl,iyh,izl_p1) - phi_rho(ixl,iyh,izl_m1);
    Elhh = phi_rho(ixl,iyh,izh_p1) - phi_rho(ixl,iyh,izh_m1);
    Ehll = phi_rho(ixh,iyl,izl_p1) - phi_rho(ixh,iyl,izl_m1);
    Ehlh = phi_rho(ixh,iyl,izh_p1) - phi_rho(ixh,iyl,izh_m1);
    Ehhl = phi_rho(ixh,iyh,izl_p1) - phi_rho(ixh,iyh,izl_m1);
    Ehhh = phi_rho(ixh,iyh,izh_p1) - phi_rho(ixh,iyh,izh_m1);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ez = El + wzh*(Eh-El);
  }
  else if (izh >= nz+zguard_size){
    // Backwards difference to get Ez at nz-1 since we don't know rho(nz)
    // Get the needed indices for interpolation
    int izl_m1 = izl-1;
    int izl_p1 = izh;
    int izh_m1 = izl;
    int izh_m2 = izl_m1;

    // Backwards difference the z=h variables since we don't know rho(nz)
    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = phi_rho(ixl,iyl,izl_p1) - phi_rho(ixl,iyl,izl_m1);
    Ellh = 3*phi_rho(ixl,iyl,izh) - 4*phi_rho(ixl,iyl,izh_m1)
      + phi_rho(ixl,iyl,izh_m2);
    Elhl = phi_rho(ixl,iyh,izl_p1) - phi_rho(ixl,iyh,izl_m1);
    Elhh = 3*phi_rho(ixl,iyh,izh) - 4*phi_rho(ixl,iyh,izh_m1)
      + phi_rho(ixl,iyh,izh_m2);
    Ehll = phi_rho(ixh,iyl,izl_p1) - phi_rho(ixh,iyl,izl_m1);
    Ehlh = 3*phi_rho(ixh,iyl,izh) - 4*phi_rho(ixh,iyl,izh_m1)
      + phi_rho(ixh,iyl,izh_m2);
    Ehhl = phi_rho(ixh,iyh,izl_p1) - phi_rho(ixh,iyh,izl_m1);
    Ehhh = 3*phi_rho(ixh,iyh,izh) - 4*phi_rho(ixh,iyh,izh_m1)
      + phi_rho(ixh,iyh,izh_m2);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ez = El + wzh*(Eh-El);       
  }
  else{
    // Forward Difference in z since we don't know izl_m1
    int izh_m1 = izl;
    int izh_p1 = izh+1;
    int izl_p1 = izh;
    int izl_p2 = izh_p1;

    // Calculate Ezi, the linearly-weighted Ez at the particle: 
    Elll = -3*phi_rho(ixl,iyl,izl) + 4*phi_rho(ixl,iyl,izl_p1)
      -phi_rho(ixl,iyl,izl_p2);
    Ellh = phi_rho(ixl,iyl,izh_p1) - phi_rho(ixl,iyl,izh_m1);
    Elhl = -3*phi_rho(ixl,iyh,izl) + 4*phi_rho(ixl,iyh,izl_p1)
      -phi_rho(ixl,iyh,izl_p2);
    Elhh = phi_rho(ixl,iyh,izh_p1) - phi_rho(ixl,iyh,izh_m1);
    Ehll = -3*phi_rho(ixh,iyl,izl) + 4*phi_rho(ixh,iyl,izl_p1)
      -phi_rho(ixh,iyl,izl_p2);
    Ehlh = phi_rho(ixh,iyl,izh_p1) - phi_rho(ixh,iyl,izh_m1);
    Ehhl = -3*phi_rho(ixh,iyh,izl) + 4*phi_rho(ixh,iyh,izl_p1)
      -phi_rho(ixh,iyh,izl_p2);
    Ehhh = phi_rho(ixh,iyh,izh_p1) - phi_rho(ixh,iyh,izh_m1);
    Ell = Elll + wyh*(Elhl-Elll);
    Ehl = Ehll + wyh*(Ehhl-Ehll);
    Elh = Ellh + wyh*(Elhh-Ellh);
    Ehh = Ehlh + wyh*(Ehhh-Ehlh);
    El = Ell + wxh*(Ehl-Ell);
    Eh = Elh + wxh*(Ehh-Elh);
    Ez = El + wzh*(Eh-El);
  }

}


#endif
#endif
