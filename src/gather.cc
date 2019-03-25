/* 3 Routines below - One for each dimension

   Calculate ne on the grid by gathering particle positions and
   storing them in density.

   Assumes position points are entirely inside den and 
   peroiodic boundary conditions

*/

// 1-D Gather 
#include "eppic.h"
#if NDIM == 1
void gather(FArrayND &den, PTYPEAVec &x, const int np, FTYPE n0)
{

  // Local variables 
  int i,ixl,ixh;
  FTYPE wxh,wxl;
  FTYPE whh,whl,wlh,wll;
  //  Note den array size is padded for fftw.  Use the global nx, ny, nz
  //  const FTYPE nx=den.xsize(), ny=den.ysize(); 
  const FTYPE nscale = (nx)/((FTYPE) np) * n0;

  // Zero the density array 
  den=0;

  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    // Define the four nearest grid points: 
    if (x(i) >= 0) { // This excludes absent particles which should always be <0
      ixl = (int) x(i);
      ixh = ixl + 1;
      if (ixh == nx) ixh = 0;
      
      // and the corresponding linear weighting factors: 
      wxh = x(i) - ixl;
      wxl = 1. - wxh;
      
      // Add this particle's contribution to den: 
      den(ixh) += wxh;
      den(ixl) += wxl; 
    }
    
  } // end for (i = 0; i < np; ++i) 
  
  // Express in physical units: 
  den *= nscale;
  
} // End 1-D gather 

#endif


// 2-D Gather
#if NDIM == 2
void gather(FArrayND &den, PTYPEAVec &x,   const int np, PTYPEAVec &y, FTYPE n0)
{

  // Local variables 
  int i,ixl,ixh,iyl,iyh;
  FTYPE wxh,wxl,wyh,wyl;
  FTYPE whh,whl,wlh,wll;
  //  Note den array size is padded for fftw.  Use the global nx, ny, nz
  //  const FTYPE nx=den.xsize(), ny=den.ysize(); 
  const FTYPE nscale = (nx*ny)/((FTYPE) np) * n0;

  // Zero the density array 
  den=0;

  //  For each particle ... 
  for (i = 0; i < np; ++i) {

    if (x(i) >= 0) { // This excludes absent particles which should always be <0

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
    den(ixh,iyh) += whh;
    den(ixh,iyl) += whl; 
    den(ixl,iyh) += wlh; 
    den(ixl,iyl) += wll; 

    }
  } // end for (i = 0; i < np; ++i) 

  // Express in physical units: 
  den *= nscale;



} // End gather 

#endif


/* 3-D Gather
*/

#if NDIM == 3

void gather(FArrayND &den, PTYPEAVec &x, PTYPEAVec &y, PTYPEAVec &z, const int np,
	      FTYPE scale)
{

  // Local variables 
  int ixl,ixh,iyl,iyh,izl,izh;
  FTYPE wxh,wxl,wyh,wyl,wzh,wzl;
  FTYPE w;

  //  Note den array size is padded for fftw. Use the global nx, ny, nz
  //  const int nx=den.xsize(), ny=den.ysize(), nz=den.zsize();
  const FTYPE nscale=(nx*ny*nz)/((FTYPE) np) * scale;

  // Zero the density array 
  den=0;

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {

    if (x(i) >= 0) { // This excludes absent particles which should always be <0

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

    }

  } // end for (i = 0; i < np; ++i) 

  // Express in physical units: 
  den *= nscale;


} // End gather 

  

#endif
