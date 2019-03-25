// advance density and velocities one step
// does not work in 1D yet

#include <stdio.h>
#include <cmath>
#include <math.h>
#include "eppic.h"

void fluid_den_meNon0_step(fluid &fspecie, field &Efield, fluid &f0, fluid &f1,
			   fluid &f2,fluid &f_old, bool store_f, FTYPE dt2, 
			   FTYPE dt3)
{

  //f0  is used to evaluate fluid components on the RHS of the equations
  //f1  is the target fluid component where the result go (cannot be the same structure as f0)
  //f2  is the fluid components where the time step initiates 

  // An external, Ey0_external, driver needs scaling for rapid use
  static const FTYPE Ex0dx=Ex0_external*(-2.0*dx);
  static const FTYPE Ey0dy=Ey0_external*(-2.0*dy);


  static const FTYPE dx2 = 1.0/(2.0*dx);
  static const FTYPE dy2 = 1.0/(2.0*dy);
  static const FTYPE dz2 = 1.0/(2.0*dz);

  static const FTYPE dx4 = 1.0 / (pow(dx, 4) );
  static const FTYPE dy4 = 1.0 / (pow(dy, 4) );
  static const FTYPE dz4 = 1.0 / (pow(dz, 4) );
  static const int xrange = nx+denx_guard_size;
#ifdef HYPER_DIFF
  static const int xm2shift = 2*denx_guard_size + nx - 2;
#endif

  // here I add the vx, vy stuff 
  FTYPE vt2 =fspecie.gamma*fspecie.T0/fspecie.m;
	//      not needed
	// 	static const FTYPE diffzv=t2/fspecie.nu; 
  FTYPE q_m =fspecie.q/fspecie.m;///fspecie.nu; This was dropped, not sure why it's here
#ifdef CHECK
  FTYPE vx_max=0;
  FTYPE vy_max=0;
  FTYPE vgrid_max=dx;
  if (ndim >= 2) vgrid_max=max(vgrid_max,dy);
  if (ndim >= 3) vgrid_max=max(vgrid_max,dz);
  vgrid_max /= dt;
#endif

  int ixm1 = (2*denx_guard_size + nx - 1)%(nx+denx_guard_size)
    -denx_guard_size;
  for (int ix=0; ix < nx; ix++) {
    int ixp1 = (ix + 1)%(xrange);
#ifdef HYPER_DIFF
    int ixp2 = (ix + 2)%(xrange);
    int ixm2 = (ix + xm2shift)%(xrange)
      -denx_guard_size;
#endif


    int iym1 = ny - 1;
    for (int iy=0; iy<ny; iy++) {
      int iyp1 = (iy + 1)%(ny);
#ifdef HYPER_DIFF
      int iyp2 = (iy + 2)%(ny);
      int iym2 = (iy - 2 + ny)%ny;
#endif
#if NDIM == 3	

      int izm1 = nz - 1;
      for (int iz=0; iz<nz; iz++) {
	int izp1 = (iz + 1)%nz;
	// calculate efields from phi_rho
#endif	  
	FTYPE Ex = (Efield.phi_rho(INDICIES(ixp1,iy,iz))
		    - Efield.phi_rho(INDICIES(ixm1,iy,iz))+Ex0dx)*(-dx2);
	FTYPE Ey = (Efield.phi_rho(INDICIES(ix,iyp1,iz)) 
		    - Efield.phi_rho(INDICIES(ix,iym1,iz))+Ey0dy)*(-dy2);      
	FTYPE Ez = (Efield.phi_rho(INDICIES(ix,iy,izp1)) - 
		    Efield.phi_rho(INDICIES(ix,iy,izm1)))*(-dz2);
	
	FTYPE dxnvx= (+f0.den(INDICIES(ixp1,iy,iz))*f0.vx(INDICIES(ixp1,iy,iz))
		      -f0.den(INDICIES(ixm1,iy,iz))*f0.vx(INDICIES(ixm1,iy,iz))
		      ) * dx2;
	FTYPE dynvy= (+ f0.den(INDICIES(ix,iyp1,iz))*f0.vy(INDICIES(ix,iyp1,iz)) 
		      - f0.den(INDICIES(ix,iym1,iz))*f0.vy(INDICIES(ix,iym1,iz)) 
		      ) * dy2;

	FTYPE fne = -dxnvx - dynvy;
#ifdef HYPER_DIFF
	
	FTYPE d4ne = (+dx4 * (+ 1 * f0.den(INDICIES(ixp2,iy,iz))
			      - 4 * f0.den(INDICIES(ixp1,iy,iz))
			      + 6 * f0.den(INDICIES(ix  ,iy,iz))
			      - 4 * f0.den(INDICIES(ixm1,iy,iz))
			      + 1 * f0.den(INDICIES(ixm2,iy,iz)) )
		      +dy4 * (+ 1 * f0.den(INDICIES(ix,iyp2,iz))
			      - 4 * f0.den(INDICIES(ix,iyp1,iz))
			      + 6 * f0.den(INDICIES(ix,iy  ,iz))
			      - 4 * f0.den(INDICIES(ix,iym1,iz))
			      + 1 * f0.den(INDICIES(ix,iym2,iz)) ) );
	
	fne -= fspecie.diffc * d4ne;

#endif

	FTYPE dnedx = (f0.den(INDICIES(ixp1,iy,iz))-f0.den(INDICIES(ixm1,iy,iz)))*dx2;
	FTYPE dnedy = (f0.den(INDICIES(ix,iyp1,iz))-f0.den(INDICIES(ix,iym1,iz)))*dy2;

	FTYPE dvxdx = (f0.vx(INDICIES(ixp1,iy,iz))-f0.vx(INDICIES(ixm1,iy,iz)))*dx2;
	FTYPE dvxdy = (f0.vx(INDICIES(ix,iyp1,iz))-f0.vx(INDICIES(ix,iym1,iz)))*dy2;
	FTYPE dvydx = (f0.vy(INDICIES(ixp1,iy,iz))-f0.vy(INDICIES(ixm1,iy,iz)))*dx2;
	FTYPE dvydy = (f0.vy(INDICIES(ix,iyp1,iz))-f0.vy(INDICIES(ix,iym1,iz)))*dy2;

	// these terms ignore the effect of B_x,y_, which I think are assumed to be 0
	FTYPE fvx =
	  //          + q_m*Efield.Ex(INDICIES(ix,iy,iz))
          + q_m*Ex
          + q_m*f0.vy(INDICIES(ix,iy,iz))*Bz
          - (f0.vx(INDICIES(ix,iy,iz))*dvxdx + f0.vy(INDICIES(ix,iy,iz))*dvxdy) 
          - vt2*dnedx/f0.den(INDICIES(ix,iy,iz))
          - fspecie.nu*f0.vx(INDICIES(ix,iy,iz));
	  
	FTYPE fvy =
	  //          + q_m*Efield.Ey(INDICIES(ix,iy,iz))
	  + q_m*Ey
          - q_m*f0.vx(INDICIES(ix,iy,iz))*Bz
          - (f0.vx(INDICIES(ix,iy,iz))*dvydx + f0.vy(INDICIES(ix,iy,iz))*dvydy)
          - vt2*dnedy/f0.den(INDICIES(ix,iy,iz))
          - fspecie.nu*f0.vy(INDICIES(ix,iy,iz));


#if NDIM == 3
	FTYPE dznvz= (+f0.den(INDICIES(ix,iy,izp1))*f0.vz(INDICIES(ix,iy,izp1))
		      -f0.den(INDICIES(ix,iy,izm1))*f0.vz(INDICIES(ix,iy,izm1))
		      ) * dz2;

#ifdef HYPER_DIFF

	int izp2 = (iz + 2)%nz;
	int izm2 = (iz - 2 + nz)%nz;
	FTYPE d4zne = dz4 * (+ 1 * f0.den(INDICIES(ix,iy,izp2))
			     - 4 * f0.den(INDICIES(ix,iy,izp1))
			     + 6 * f0.den(INDICIES(ix,iy,iz  ))
			     - 4 * f0.den(INDICIES(ix,iy,izm1))
			     + 1 * f0.den(INDICIES(ix,iy,izm2)) );

	fne += - dznvz - fspecie.diffc * d4zne;
#else
	fne += - dznvz;
#endif

	// velocity stuff
	FTYPE dvxdz = (f0.vx(INDICIES(ix,iy,izp1))-f0.vx(INDICIES(ix,iy,izm1)))*dz2;
	FTYPE dvydz = (f0.vy(INDICIES(ix,iy,izp1))-f0.vy(INDICIES(ix,iy,izm1)))*dz2;

	fvx -= f0.vz(INDICIES(ix,iy,iz))*dvxdz;
	fvy += q_m*f0.vz(INDICIES(ix,iy,iz))*Bx-f0.vz(INDICIES(ix,iy,iz))*dvydz;

	FTYPE dnedz=( f0.den(INDICIES(ix,iy,izp1))
		    -f0.den(INDICIES(ix,iy,izm1)) )*dz2;


	FTYPE dvzdx = (f0.vz(INDICIES(ixp1,iy,iz))-f0.vz(INDICIES(ixm1,iy,iz)))*dx2;
	FTYPE dvzdy = (f0.vz(INDICIES(ix,iyp1,iz))-f0.vz(INDICIES(ix,iym1,iz)))*dy2;
	FTYPE dvzdz = (f0.vz(INDICIES(ix,iy,izp1))-f0.vz(INDICIES(ix,iy,izm1)))*dz2;

	FTYPE fvz = 
	  //	  + q_m * Efield.Ez(INDICIES(ix,iy,iz))
	  + q_m * Ez
	  - vt2*dnedz/f0.den(INDICIES(ix,iy,iz))
	  - f0.vx(INDICIES(ix,iy,iz))*dvzdx
	  - f0.vy(INDICIES(ix,iy,iz))*dvzdy
	  - f0.vz(INDICIES(ix,iy,iz))*dvzdz;

	f1.vz(INDICIES(ix,iy,iz)) = 
	  f2.vz(INDICIES(ix,iy,iz)) 
	  + dt2*fvz
	  + dt3*f_old.vz(INDICIES(ix,iy,iz));

	if (store_f == true ) f_old.vz(INDICIES(ix,iy,iz)) = fvz;

#endif
#ifdef CHECK

// 	if (std::isnan(fvx)) { // we have a nan
// 	  cout << "You have a nan!\n"; // for debuging really
// 	} else if ( std::isinf(fvx)) {
// 	  cout << "You have an inf!\n"; // for debuging really
// 	}

	// Setup Courant condition
	vx_max=max(vx_max,fabs(fspecie.vx(INDICIES(ix,iy,iz))));
	vy_max=max(vy_max,fabs(fspecie.vy(INDICIES(ix,iy,iz))));
#endif

	f1.vx(INDICIES(ix,iy,iz)) = 
	  f2.vx(INDICIES(ix,iy,iz)) 
	  + dt2*fvx
	  + dt3*f_old.vx(INDICIES(ix,iy,iz));

	f1.vy(INDICIES(ix,iy,iz)) = 
	  f2.vy(INDICIES(ix,iy,iz)) 
	  + dt2*fvy
	  + dt3*f_old.vy(INDICIES(ix,iy,iz));

	f1.den(INDICIES(ix,iy,iz)) = 
	  f2.den(INDICIES(ix,iy,iz))
	  + dt2*fne 
	  + dt3*f_old.den(INDICIES(ix,iy,iz)) ;
	
	/* Store the rhs of the current time step as the old value
	   to be used in the next call */
	if (store_f == true ) {
	  f_old.den(INDICIES(ix,iy,iz)) = fne;
	  f_old.vx(INDICIES(ix,iy,iz)) = fvx;
	  f_old.vy(INDICIES(ix,iy,iz)) = fvy;
	}
	
#if NDIM == 3
	izm1 = iz;
      }
#endif
      
      iym1 = iy;
    }
    
    ixm1 = ix;
  }
#ifdef CHECK
  FTYPE vmax=sqrt(vx_max*vx_max + vy_max*vy_max);
  if (vmax > 2*vgrid_max) 
    terminate(-1,"Courant violation: fluid velocities too large");

#endif  

}

 
