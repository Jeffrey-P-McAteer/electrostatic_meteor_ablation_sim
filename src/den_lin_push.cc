/* Advance the linearized fluid densities: */
#include "eppic.h"
void den_lin_push(fluid &adv,fluid &old, fluid &cur, FTYPE alpha)
{
  int ix, ixp1, ixm1;
  int iy, iyp1, iym1;
  FTYPE div_v;
  
  /* Assume no 0 order flow */
  ixm1=nx-1;
  for (ix=0; ix<nx; ix++) {
    ixp1=(ix+1)%nx;
    iym1=ny-1;
#if NDIM==1
    div_v =cur.vx(ixp1)-cur.vx(ixm1);
#elif NDIM==2      
    for (iy=0; iy<ny; iy++) {
      iyp1=(iy+1)%ny;
      div_v =(cur.vx(ixp1,iy)-cur.vx(ixm1,iy));
      div_v+=(cur.vy(ix,iyp1)-cur.vy(ix,iym1));
      adv.den(ix,iy)=old.den(ix,iy)-div_v/2.;
#else
    for (iy=0; iy<ny; iy++) {
      iyp1=(iy+1)%ny;
      int iz, izp1, izm1;
      izm1=nz-1;
      for (iz=0; iz<nz; iz++) {
	izp1=(iz+1)%nz;
	div_v =(cur.vx(ixp1,iy,iz)-cur.vx(ixm1,iy,iz));
	div_v+=(cur.vy(ix,iyp1,iz)-cur.vy(ix,iym1,iz));
	div_v+=(cur.vy(ix,iy,izp1)-cur.vy(ix,iy,izm1));
	adv.den(ix,iy,iz)=old.den(ix,iy,iz)-div_v/2.;
	izm1=iz;
      }
    }
#endif
#if NDIM==2      
      iym1=iy;
    }
#endif
    ixm1=ix;
  } 
  
}

