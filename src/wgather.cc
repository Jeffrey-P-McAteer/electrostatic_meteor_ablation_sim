/* 3 Routines below for 
   Weighted gathering of particle positions where the weight is passed in.
   Includes overloaded functions - 1D, 2D and 3D!
*/
#include "eppic.h"

#if NDIM == 1
// 1-D wgather
void wgather(FArrayND &den, PTYPEAVec &x, PTYPEAVec &wt, FTYPE scale)
{
  
  // Local variables 
  int i,ixl,ixh;
  FTYPE wxh,wxl;
  FTYPE v;
  int np=x.size();

  const FTYPE nscale = (nx)/((FTYPE) np) * scale;
  
  // Clear the density array 
  den=0;
    
  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    // Calculate the weights assuming Gaussian distributions: 

    // Define the four nearest grid points: 
    ixl = (int) x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;

    // and the corresponding linear weighting factors: 
    wxh = x(i) - ixl;
    wxl = 1. - wxh;

    // Add this particle's contribution to den: 
    v=wt(i);
    den(ixh) += wxh * v;
    den(ixl) += wxl * v;

  } // end for (i = 0; i < np; ++i) 

  // Express in physical units: 
  den *= nscale;

  
  return;

}
#endif

#if NDIM == 2
// 2-D wgather
void wgather(FArrayND &den, PTYPEAVec &x, PTYPEAVec &y, PTYPEAVec &wt, 
	      FTYPE scale)
{
  // Local variables 
  
  int i,ixl,ixh,iyl,iyh;
  FTYPE wxh,wxl,wyh,wyl;
  FTYPE whh,whl,wlh,wll;
  FTYPE v;
  int np=x.length();

  const FTYPE nscale = (nx*ny)/((FTYPE) np) * scale;
  
  // Clear the density array 
  den=0;
    
  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    // Calculate the weights assuming Gaussian distributions: 

    // Define the four nearest grid points: 
    ixl = (int) x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;
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

    // Add this particle's contribution to den: 
    v=wt(i);
    den(ixh,iyh) += whh * v;
    den(ixh,iyl) += whl * v;
    den(ixl,iyh) += wlh * v;
    den(ixl,iyl) += wll * v;

  } // end for (i = 0; i < np; ++i) 

  // Express in physical units: 
  den *= nscale;

  
  return;

}
#endif

#if NDIM == 3
// 3-D wgather
void wgather(FArrayND &den, PTYPEAVec &x, PTYPEAVec &y, PTYPEAVec &z, 
	     PTYPEAVec &wt, FTYPE scale)
{
  
  // Local variables 
  int ixl,ixh,iyl,iyh,izl,izh;
  FTYPE wxh,wxl,wyh,wyl,wzh,wzl;
  FTYPE w, v;
  const int np=x.length();

  const FTYPE nscale=(nx*ny*nz)/((FTYPE) np) * scale;

  // Zero the density array 
  den=0;

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {

    // Define the four nearest grid points: 
    ixl = (int) x(i);
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;

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
    
    v=wt(i);
    w = wxh * wyh * wzh;
    den(ixh,iyh,izh) += w * v;
    w = wxh * wyh * wzl;
    den(ixh,iyh,izl) += w * v;
    w = wxh * wyl * wzh;
    den(ixh,iyl,izh) += w * v;
    w = wxh * wyl * wzl;
    den(ixh,iyl,izl) += w * v;

    w = wxl * wyh * wzh;
    den(ixl,iyh,izh) += w * v;
    w = wxl * wyh * wzl;
    den(ixl,iyh,izl) += w * v;
    w = wxl * wyl * wzh;
    den(ixl,iyl,izh) += w * v;
    w = wxl * wyl * wzl;
    den(ixl,iyl,izl) += w * v;

  } // end for (i = 0; i < np; ++i) 

  // Express in physical units: 
  den *= nscale;

  
} // End gather 

#endif  


