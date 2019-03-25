// Take the gradient of array phi and place it in Ex, Ey and Ez.
#include "eppic.h"
//#include <iostream>

void gradient(FArrayND &phi, INDICIES(FArrayND &Ex, FArrayND &Ey, FArrayND &Ez), 
	      INDICIES(FTYPE dx, FTYPE dy, FTYPE dz), FTYPE scaler)
{
  
  FTYPE dx_scale=scaler/(2*dx);
  FTYPE dy_scale=scaler/(2*dy);
#if NDIM==3
  FTYPE dz_scale=scaler/(2*dz);
#endif
  
  int ixm1=(nx-1);
  for (int ix=0; ix<nx; ix++) {
    int ixp1=(ix+1)%nx;
    int iym1=ny-1;
    for (int iy=0; iy<ny; iy++) {
      int iyp1=(iy+1)%ny;
#if NDIM==3
      int izm1=nz-1;
#endif
      for (int iz=0; iz<nz; iz++) {
        
        Ex(INDICIES(ix, iy, iz)) = (phi(INDICIES(ixp1,iy,iz))-phi(INDICIES(ixm1,iy,iz)))*dx_scale;
          
        /*if (boundary_type==INJECT & ix==0) Ex(INDICIES(0, iy, iz)) = phi(INDICIES(1,iy,iz))*dx_scale;
        if (boundary_type==INJECT & ix==nx-1) Ex(INDICIES(nx-1, iy, iz)) = -phi(INDICIES(nx-2,iy,iz))*dx_scale;
        if (boundary_type==INJECT & (ix=!0 & ix=!nx-1)) Ex(INDICIES(ix, iy, iz)) = (phi(INDICIES(ixp1,iy,iz))-phi(INDICIES(ixm1,iy,iz)))*dx_scale;
        
        if (boundary_type==PERIODIC & (ix=!0 & ix=!nx-1)) Ex(INDICIES(ix, iy, iz)) = (phi(INDICIES(ixp1,iy,iz))-phi(INDICIES(ixm1,iy,iz)))*dx_scale;
        */

#if NDIM > 1	
	if (Ey.length()>0)
	  Ey(INDICIES(ix, iy, iz)) = (phi(INDICIES(ix,iyp1,iz))-
				      phi(INDICIES(ix,iym1,iz)))*dy_scale;
#if NDIM==3
	int izp1 = (iz+1)%nz;
	if (Ez.length()>0)
	  Ez(INDICIES(ix, iy, iz)) = (phi(INDICIES(ix,iy,izp1))-
				      phi(INDICIES(ix,iy,izm1)))*dz_scale;
	izm1=iz;
#endif
#endif
      }
      iym1=iy;
    }
    ixm1=ix;
  }
  
  //X boundary correction of the E field. THIS IS WEIRD! WHY ny*boundary_type??? 
  for (int iy=0; iy<ny*boundary_type[0]; iy++) {
    for (int iz=0; iz<nz*boundary_type[0]; iz++) {
      Ex(INDICIES(0, iy, iz)) += phi(INDICIES(nx-1,iy,iz))*dx_scale;
      Ex(INDICIES(nx-1, iy, iz)) -= phi(INDICIES(1,iy,iz))*dx_scale;
    }
  }
  //cout<<"phi   Ex"<<endl;
  //for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) cout<<phi(INDICIES(ix,iy,0))<<" "<<Ex(INDICIES(ix,iy,0))<<endl;  
  
}

