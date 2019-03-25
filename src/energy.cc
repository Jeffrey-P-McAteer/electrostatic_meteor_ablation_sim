/* Calculate the perturbed energy contained in the simulation*/

#include <math.h>
#include <stdio.h> 
#include <string.h> 
#include "eppic.h"
#include <iostream>

FTYPE energy(field &Efield, particle_dist *pic, FTYPE &Wf, 
	      FTYPEAVec &Wp)
{

  /* Total Field energy: (Divide by volume to obtain energy density)
     (This could be made to represent the inter-grid areas better .) */
  Wf = 0.;
  int nx2=Efield.phi_rho.xsize()-Efield.phi_rho.xstart();
  int ixm1=-1;
  if (ixm1 < Efield.phi_rho.xstart()) ixm1=nx-1; // This assumes it's periodic (awkward)
  for (int ix=0; ix<nx; ix++) {
    int ixp1=(ix+1)%nx2;
    int iym1=ny-1;
    for (int iy=0; iy<ny; iy++) { 
      int iyp1=(iy+1)%ny;
#if NDIM==3
      int izm1=nz-1;
#endif
      for (int iz=0; iz<nz; iz++) {
	Wf += Sqr((Efield.phi_rho(INDICIES(ixp1,iy,iz)) -
		   Efield.phi_rho(INDICIES(ixm1,iy,iz)))/(2*dx));
#if NDIM > 1
	Wf += Sqr((Efield.phi_rho(INDICIES(ix,iyp1,iz)) -
		   Efield.phi_rho(INDICIES(ix,iym1,iz)))/(2*dy));
#endif
#if NDIM== 3
	int izp1=(iz+1)%nz;
	Wf += Sqr((Efield.phi_rho(INDICIES(ix,iy,izp1)) -
		   Efield.phi_rho(INDICIES(ix,iy,izm1)))/(2*dz));
	izm1=iz;
#endif
      }
      iym1=iy;
    }
    ixm1=ix;
  }
  if (nx > 1) Wf *= dx;
  if (ny > 1) Wf *= dy;
  if (nz > 1) Wf *= dz;
  Wf *= eps/2.;

  /* Particle energy */

  FTYPE Wp_total = 0.;
  for (int id=0; id<ndist; ++id) {
    Wp[id]=0;
    if (method[id] >= 0 ) {
      if (method[id] == 0) {
	for (int i=0;i<pic[id].np;i++) {
	  Wp[id] += Sqr(pic[id].vx(i)*dx/dt);
	  //cout<<id<<" "<<pic[id].vx(i)<<endl;
	  if (vel_dim[id] >= 2) Wp[id] += Sqr(pic[id].vy(i)*dy/dt);
	  //cout<<id<<" "<<Sqr(pic[id].vx(i)*dx/dt)<<" "<<Sqr(pic[id].vy(i)*dy/dt)<<" "<<Wp[id]<<endl;
	  if (vel_dim[id] == 3) Wp[id] += Sqr(pic[id].vz(i)*dz/dt);
	}
	//cout<<id<<Wp[id]<< " "<<Wf<<endl;
	/*
	cout << "\nProc: " << mpi_rank << " id: " 
	     << id << " Wp[id]: " << Wp[id] << "\n";
      
	if (isnan(Wp[id])) {
	  cout << "nabsent = " << pic[id].nabsent << "\n";
	  for (int iabsent=0;iabsent<pic[id].nabsent;iabsent++) {
	    cout << " iabsent=" << iabsent 
		 << " absent_idx=" << pic[id].absent(iabsent)
		 << "\n";
	  }

	  Wp[id]=0;
	  for (int i=0;i<pic[id].np;i++) {
	    Wp[id]=0;
	    Wp[id] += Sqr(pic[id].vx(i)*dx/dt);
	    //cout<<i<<" "<<pic[id].vx(i)<<endl;
	    if (vel_dim[id] >= 2) Wp[id] += Sqr(pic[id].vy(i)*dy/dt);
	    //cout<<i<<" "<<Sqr(pic[id].vx(i)*dx/dt)<<" "<<Sqr(pic[id].vy(i)*dy/dt)<<" "<<Wp[id]<<endl;
	    if (vel_dim[id] == 3) Wp[id] += Sqr(pic[id].vz(i)*dz/dt);
	    if (isnan(Wp[id])) {
	      cout << "Problem particle = " << i 
		   << " x= " << pic[id].x(i) 
		   << " y= " << pic[id].y(i) 
		   << " vx= " << pic[id].vx(i) 
		   << " vy= " << pic[id].vy(i) 
		   << "\n";
	    }
	  }
	  
	  }
	*/
	// This normalization means we are calculating the total particle 
	// energy:
	if (pic[id].np > 0) Wp[id] *= 0.5 * pic[id].m * (pic[id].n0avg/pic[id].np);
	if (nx > 1) Wp[id] *= nx*dx;
	if (ny > 1) Wp[id] *= ny*dy;
	if (nz > 1) Wp[id] *= nz*dz;
      }
      
      Wp_total += Wp[id];
      
    } else {
      // fluid
    }
    
  }

  return Wp_total;
  
}

	
