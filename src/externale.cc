#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

// Impose an external Electric field
void externale(field &Efield)
{

  /*Add a DC Ey field (if the E field is know at this point) */
  if (Efield.Ey.size() == Efield.phi_rho.size()) Efield.Ey += Ey0_external;

  /* Add some random noise - Note 3D vpush routines will not see this */
  if (Ermag>0) {
    int ix;
#if NDIM==2
    int iy;
#endif
#if NDIM < 3
    static int iran = -15-mpi_rank;
#endif
    for (ix=0 ; ix<nx ; ++ix) {

#if NDIM == 1
      Efield.Ex(ix) += Ermag*(ran3(&iran)-0.5)*2.;
#endif
#if NDIM == 2
      for (iy=0 ; iy<ny ; ++iy) {
	if (nx>1) Efield.Ex(ix,iy) += Ermag*(ran3(&iran)-0.5)*2.;
	if (ny>1) Efield.Ey(ix,iy) += Ermag*(ran3(&iran)-0.5)*2.; 
      } 
#endif
    }
  }

  /* Add a particlular mode */
  /*  FTYPE Ewmag=0.0;
  kx=2*PI/(nx*dx);
  //  ky=2*PI/(ny*dy);
  w= sqrt( (n0d[0]*Sqr(qd[0])) / (eps*md[0]) );
  FTYPE slope = Ewmag*(it+1)/min(nEext*1.,256.);
  if (slope > Ewmag) slope=Ewmag;

  for (ix=0 ; ix<nx ; ++ix) 
    for (iy=0 ; iy<ny ; ++iy) {
    if (kx >0.)  Ex(ix,iy) += slope * sin(ix*kx*dx-w*it*dt);
    //    if (ky >0.) Ey(ix,iy) = slope * sin(iy*ky*dy-w*it*dt);
    }
    return;
    */    

}

