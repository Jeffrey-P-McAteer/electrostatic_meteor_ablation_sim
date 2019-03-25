/* Gather number density and flux in one loop.
   Modeled after gather_domain.cc. and wgather.cc
   This routine DOES NOT scale the arrays, since 
   that's a single mutliplication that the user can 
   do outside this routine. */

// -->Consider just passing in particle_dist &pic 
//    and getting x,y,z as in gather_weight.cc

#include "eppic.h"

#undef __FUNCT__
#define __FUNCT__ "gather_den_flux"

#if NDIM == 1
void gather_den_flux(FArrayND &den, 
		     FArrayND &xflux, 
		     // particle_dist &pic,
		     PTYPEAVec &x,
		     PTYPEAVec &xwt, 
		     const int np,
		     const int nxmin)
{

  // Local variables 
  int i,ixl,ixh;
  FTYPE wxh,wxl;
  FTYPE v;
  // const int np=pic.np;
  // PTYPEAVec x=pic.x;

  // Zero the density and flux arrays 
  den = 0.0;
  xflux = 0.0;
  long long np = x.size();

  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    // Define the four nearest grid points: 
    if (x(i) >= 0) { // This excludes absent particles which should always be <0
      ixl = (int) x(i);
      ixh = ixl + 1;
      
      // and the corresponding linear weighting factors: 
      wxh = x(i) - ixl;
      wxl = 1. - wxh;
      
      // Update density
      den(ixh) += wxh;
      den(ixl) += wxl; 

      // Update X flux
      v = xwt(i);
      xflux(ixh,iyh) += wxh*v;
      xflux(ixh,iyl) += wxl*v; 

    } else if (x(i) >= nxmin) {
      ixh = 0;
      wxh = x(i) + 1.0;

      // Update density
      den(ixh) += wxh;

      // Update X flux
      v = xwt(i);
      xflux(ixh,iyh) += wxh*v;
    }
  } // end for (i = 0; i < np; ++i) 
} // 1-D block
#endif

#if NDIM == 2
void gather_den_flux(FArrayND &den, 
		     FArrayND &xflux, FArrayND &yflux, 
		     // particle_dist &pic,
		     PTYPEAVec &x, PTYPEAVec &y,
		     PTYPEAVec &xwt, PTYPEAVec &ywt,
		     const int np,
		     const int nxmin)
{

  // Local variables 
  int i,ixl,ixh,iyl,iyh;
  FTYPE wxh,wxl,wyh,wyl;
  FTYPE whh,whl,wlh,wll;
  FTYPE v;
  int ny=den.size(1);
  // const int np=pic.np;
  // PTYPEAVec *data[NDIM]={INDICIES(&pic.x,&pic.y,&pic.z)};
  // PTYPEAVec x=pic.x,y=pic.y;

  // Zero the density and flux arrays 
  den = 0.0;
  xflux = 0.0;
  yflux = 0.0;
  int ny = den.size(1);
  long long np = x.size();

  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    if (x(i) >= 0) { // This excludes absent particles which should always be <0

      // Define the four nearest grid points: 
      ixl = (int) x(i);
      ixh = ixl + 1;
      
      iyl = (int) y(i);
      iyh = iyl + 1;
      if (iyh == ny) iyh = 0;
      
      // and the corresponding linear weighting factors: 
      wxh = x(i) - ixl;
      wxl = 1. - wxh;
      wyh = y(i) - iyl;
      wyl = 1. - wyh;
      
      whh = wxh * wyh;
      whl = wxh * wyl;
      wlh = wxl * wyh;
      wll = wxl * wyl;

      // Update density
      den(ixh,iyh) += whh;
      den(ixh,iyl) += whl; 
      den(ixl,iyh) += wlh; 
      den(ixl,iyl) += wll; 

      // Update X flux
      v = xwt(i);
      xflux(ixh,iyh) += whh*v;
      xflux(ixh,iyl) += whl*v; 
      xflux(ixl,iyh) += wlh*v;
      xflux(ixl,iyl) += wll*v; 

      // Update Y flux
      v = ywt(i);
      yflux(ixh,iyh) += whh*v;
      yflux(ixh,iyl) += whl*v; 
      yflux(ixl,iyh) += wlh*v; 
      yflux(ixl,iyl) += wll*v; 

    } else if (x(i) >= nxmin) {
      ixh = 0;
      iyl = (int) y(i);
      iyh = iyl + 1;
      if (iyh == ny) iyh = 0;

      wxh = x(i) + 1.0;
      wyh = y(i) - iyl;
      wyl = 1. - wyh;

      whh = wxh * wyh;
      whl = wxh * wyl;

      // Add this particle's contribution to den: 
      den(ixh,iyh) += whh;
      den(ixh,iyl) += whl; 
      v = xwt(i);
      xflux(ixh,iyh) += whh*v;
      xflux(ixh,iyl) += whl*v; 
      v = ywt(i);
      yflux(ixh,iyh) += whh*v;
      yflux(ixh,iyl) += whl*v; 
    }

  } // end for (i = 0; i < np; ++i) 

} // 2-D block
#endif

#if NDIM == 3
void gather_den_flux(FArrayND &den, 
		     FArrayND &xflux, FArrayND &yflux, FArrayND &zflux, 
		     // particle_dist &pic,
		     PTYPEAVec &x, PTYPEAVec &y, PTYPEAVec &z,
		     PTYPEAVec &xwt, PTYPEAVec &ywt, PTYPEAVec &zwt,
		     const int np_dummy,
		     const int nxmin)
{

  // Local variables 
  int i,ixl,ixh,iyl,iyh,izl,izh;
  FTYPE wxh,wxl,wyh,wyl,wzh,wzl;
  FTYPE whh,whl,wlh,wll;
  FTYPE v,w;
  int ny = den.size(1);
  int nz = den.size(2);
  // const int np=pic.np;
  // PTYPEAVec *data[NDIM]={INDICIES(&pic.x,&pic.y,&pic.z)};
  // PTYPEAVec x=pic.x,y=pic.y,z=pic.z;

  // Zero the density and flux arrays 
  den = 0.0;
  xflux = 0.0;
  yflux = 0.0;
  zflux = 0.0;
  long long np = x.size();

  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    if (x(i) >= 0) { // This excludes absent particles which should always be <0

    // Define the four nearest grid points: 
    ixl = (int) x(i);
    ixh = ixl + 1;

    iyl = (int) y(i);
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;

    izl = (int) z(i);
    izh = izl + 1;
    if (izh == nz) izh = 0;

    // and the corresponding linear weighting factors: 
    wxh = x(i) - ixl;
    wxl = 1. - wxh;

    wyh = y(i) - iyl;
    wyl = 1. - wyh;

    wzh = z(i) - izl;
    wzl = 1. - wzh;
    
    // Update density
    w = wxh * wyh * wzh;
    den(ixh,iyh,izh) += w;
    w = wxh * wyh * wzl;
    den(ixh,iyh,izl) += w;
    w = wxh * wyl * wzh;
    den(ixh,iyl,izh) += w;
    w = wxh * wyl * wzl;
    den(ixh,iyl,izl) += w;
    w = wxl * wyh * wzh;
    den(ixl,iyh,izh) += w;
    w = wxl * wyh * wzl;
    den(ixl,iyh,izl) += w;
    w = wxl * wyl * wzh;
    den(ixl,iyl,izh) += w;
    w = wxl * wyl * wzl;
    den(ixl,iyl,izl) += w;

    // Update X flux
    v = xwt(i);
    w = wxh * wyh * wzh;
    xflux(ixh,iyh,izh) += w*v;
    w = wxh * wyh * wzl;
    xflux(ixh,iyh,izl) += w*v;
    w = wxh * wyl * wzh;
    xflux(ixh,iyl,izh) += w*v;
    w = wxh * wyl * wzl;
    xflux(ixh,iyl,izl) += w*v;
    w = wxl * wyh * wzh;
    xflux(ixl,iyh,izh) += w*v;
    w = wxl * wyh * wzl;
    xflux(ixl,iyh,izl) += w*v;
    w = wxl * wyl * wzh;
    xflux(ixl,iyl,izh) += w*v;
    w = wxl * wyl * wzl;
    xflux(ixl,iyl,izl) += w*v;

    // Update Y flux
    v = ywt(i);
    w = wxh * wyh * wzh;
    yflux(ixh,iyh,izh) += w*v;
    w = wxh * wyh * wzl;
    yflux(ixh,iyh,izl) += w*v;
    w = wxh * wyl * wzh;
    yflux(ixh,iyl,izh) += w*v;
    w = wxh * wyl * wzl;
    yflux(ixh,iyl,izl) += w*v;
    w = wxl * wyh * wzh;
    yflux(ixl,iyh,izh) += w*v;
    w = wxl * wyh * wzl;
    yflux(ixl,iyh,izl) += w*v;
    w = wxl * wyl * wzh;
    yflux(ixl,iyl,izh) += w*v;
    w = wxl * wyl * wzl;
    yflux(ixl,iyl,izl) += w*v;

    // Update Z flux
    v = zwt(i);
    w = wxh * wyh * wzh;
    zflux(ixh,iyh,izh) += w*v;
    w = wxh * wyh * wzl;
    zflux(ixh,iyh,izl) += w*v;
    w = wxh * wyl * wzh;
    zflux(ixh,iyl,izh) += w*v;
    w = wxh * wyl * wzl;
    zflux(ixh,iyl,izl) += w*v;
    w = wxl * wyh * wzh;
    zflux(ixl,iyh,izh) += w*v;
    w = wxl * wyh * wzl;
    zflux(ixl,iyh,izl) += w*v;
    w = wxl * wyl * wzh;
    zflux(ixl,iyl,izh) += w*v;
    w = wxl * wyl * wzl;
    zflux(ixl,iyl,izl) += w*v;

    } else if (x(i) >= nxmin) {
      ixh = 0;

      iyl = (int) y(i);
      iyh = iyl + 1;
      if (iyh == ny) iyh = 0;

      izl = (int) z(i);
      izh = izl + 1;
      if (izh == nz) izh = 0;

      // and the corresponding linear weighting factors: 
      wxh = x(i) + 1.0;

      wyh = y(i) - iyl;
      wyl = 1. - wyh;

      wzh = z(i) - izl;
      wzl = 1. - wzh;
    
      // Update density
      w = wxh * wyh * wzh;
      den(ixh,iyh,izh) += w;
      w = wxh * wyh * wzl;
      den(ixh,iyh,izl) += w;
      w = wxh * wyl * wzh;
      den(ixh,iyl,izh) += w;
      w = wxh * wyl * wzl;
      den(ixh,iyl,izl) += w;

      // Update X flux
      v = xwt(i);
      w = wxh * wyh * wzh;
      xflux(ixh,iyh,izh) += w*v;
      w = wxh * wyh * wzl;
      xflux(ixh,iyh,izl) += w*v;
      w = wxh * wyl * wzh;
      xflux(ixh,iyl,izh) += w*v;
      w = wxh * wyl * wzl;
      xflux(ixh,iyl,izl) += w*v;

      // Update Y flux
      v = ywt(i);
      w = wxh * wyh * wzh;
      yflux(ixh,iyh,izh) += w*v;
      w = wxh * wyh * wzl;
      yflux(ixh,iyh,izl) += w*v;
      w = wxh * wyl * wzh;
      yflux(ixh,iyl,izh) += w*v;
      w = wxh * wyl * wzl;
      yflux(ixh,iyl,izl) += w*v;

      // Update Z flux
      v = zwt(i);
      w = wxh * wyh * wzh;
      zflux(ixh,iyh,izh) += w*v;
      w = wxh * wyh * wzl;
      zflux(ixh,iyh,izl) += w*v;
      w = wxh * wyl * wzh;
      zflux(ixh,iyl,izh) += w*v;
      w = wxh * wyl * wzl;
      zflux(ixh,iyl,izl) += w*v;
    }

  } // end for (i = 0; i < np; ++i) 

} // 3-D block
#endif
