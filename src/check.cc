// Check the input parameters to make sure no obvious errors were made 

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "eppic.h"
#include "eppic-mpi.h"

void check(void)
{
  int id;
  FTYPE kxt, kyt, kzt, debye, wpe, Vest;

  //dls///////////////////////////////////////////////
  //Check domain decomposition
  if (mpi_np%nsubdomains != 0) {
    char string[130];
    sprintf(string,"%s%d%s%d%s\n","Error: number of processors (mpi_np = ", mpi_np,
	    ") not evenly divisible by number of subdomains (nsubdomains = ",
	    nsubdomains,")\n");
    terminate(-1,string);
  }
  //end dls///////////////////////////////////////////

  //Make sure the number of processors is consistent with the PETSc decomposition
  if (efield_algorithm == 2 && nx < mpi_np/nsubdomains) {
    char string[130];
    sprintf(string,"The PETSc solver requires nx (%d) >= mpi_np (%d) / nsubdomains (%d)\n",
	    nx,mpi_np,nsubdomains);
    terminate(-1,string);
  }

  //Check and fix dimensional problems
  if (ndim != ndim_space) {
    char string[130];
    sprintf(string,"Error: ndim_space=%d does not match compiled ndim=%d\n",
	    ndim_space,ndim);
    terminate(-1,string);
  }

  if (ndim <= 2) {
    if (nz != 1) {
      if (mpi_rank == 0) printf("\n WARNING: setting nz=1, dz=1\n\n");
      nz=1;
      dz=1;
    }
    if (ndim == 1)
      if (ny != 1) {
	if (mpi_rank == 0) printf("\n WARNING: setting ny=1, dy=1\n\n");
	ny=1;
	dy=1;
    }
  }

  // Echo some parameters related to Efield solver
  if (mpi_rank==0) cout << "Efield Algorithm = " << efield_algorithm << "\n";

  // Check Debye length conformity 
  if (mpi_rank == 0) {
    printf("\n\n;Checking a few parameters:\n");
  
    float k_b;
    if (eps!=1) k_b=KB; else k_b=1;
    for (id=0; id<ndist; ++id) {
      printf(";Distribution %d:\n", id);
      
      kxt=Sqr(vxthd[id])*md[id];
      printf(";\tX Temperature: %g (K)\n", kxt/k_b);
      
      kyt=Sqr(vythd[id])*md[id];
      printf(";\tY Temperature: %g (K)\n", kyt/k_b);

      kzt=Sqr(vzthd[id])*md[id];
      printf(";\tZ Temperature: %g (K)\n", kzt/k_b);
      
      printf(";\tn0d[%d] = %g\n",id,n0d[id]);
      if (unscale_density){
        // n0d is #real particles/sim particles, npd is # sim particles
        debye=sqrt( eps*kxt / (n0d[id]*npd[id]/(nx*ny*nz*dx*dy*dz)*Sqr(qd[id])));
      }
      else{
        debye=sqrt( eps*kxt / (n0d[id]*Sqr(qd[id])) );
      }
      printf(";\tX Debye Length: %g (m)\n", debye);
      if (dx > debye) printf(";WARNING: dx (%g) > debye length %s",dx,
			     "... Possible numerical instability\n");

      if (unscale_density){
        // n0d is the total number of initialized particles
        debye=sqrt( eps*kyt / (n0d[id]*npd[id]/(nx*ny*nz*dx*dy*dz)*Sqr(qd[id])));
      }
      else{
        debye=sqrt( eps*kyt / (n0d[id]*Sqr(qd[id])) );
      }
      printf(";\tY Debye Length: %g (m)\n", debye);
      if (dy > debye) printf(";WARNING: dy (%g) > debye length %s",dy,
			     "... Possible numerical instability\n");
      if (unscale_density){
        // n0d is the total number of initialized particles
        debye=sqrt( eps*kzt / (n0d[id]*npd[id]/(nx*ny*nz*dx*dy*dz)*Sqr(qd[id])));
      }
      else{
        debye=sqrt( eps*kzt / (n0d[id]*Sqr(qd[id])) );
      }
      printf(";\tZ Debye Length: %g (m)\n", debye);
      if (dz > debye) printf(";WARNING: dz (%g) > debye length %s",dz,
			     "... Possible numerical instability\n");
      
      if (unscale_density){
        // n0d is the total number of initialized particles
        wpe=sqrt( (n0d[id]*npd[id]/(nx*ny*nz*dx*dy*dz)*Sqr(qd[id])) / (eps*md[id]) );
      }
      else{
        wpe=sqrt( (n0d[id]*Sqr(qd[id])) / (eps*md[id]) );
      }
      printf(";\tPlasma frequency: %g (rad/s) %g (Hz)\n", wpe, wpe/(2*PI));

      FTYPE wce=qd[id]*fabs(Bz)/md[id];
      if (wce == 0 && Bx != 0) wce=qd[id]*fabs(Bx)/md[id];
      if (wce !=0 ) {
	printf(";\tCyclotron frequency: %g\n", wce);
	if (Bz*Bz > 0.0) {
	  printf(";\tEx0/Bz0: %f\n", fabs(Ex0_external/Bz));
	  printf(";\tEy0/Bz0: %f\n", fabs(Ey0_external/Bz));
	}
	if (Bx*Bx > 0.0) printf(";\tEy0/Bx0: %f\n", fabs(Ey0_external/Bx));
      }

      printf(";\tTime steps per plasma oscillation: %g\n",  1./(dt*wpe/(2*PI)));
      if (dt*wpe/(2*PI) > .25) 
	printf(";WARNING: dt (%g) cannot resolve plasma freqency\n",dt);

      if(method[id] < 0) {//fluid
	// test Courant condition
	FTYPE vgrid_max=dx;
	if (ndim >= 2) vgrid_max=max(vgrid_max,dy);
	if (ndim >= 3) vgrid_max=max(vgrid_max,dz);
	vgrid_max/=2.*dt;
	if (efield_algorithm == 2) {
	  cout << ";\tCourant violation not tested due to Efield Algorithm";
	} else {
	  if (vxthd[id] > vgrid_max) {
	    cout << "Vgrid_max = " << vgrid_max << "\n";
	    cout << "Vxthd = " << vxthd[id] << "\n";
	    terminate(-1,"Courant violation: thermal vx (vxth) too large");
	  }
	  if (ndim >= 2) {
	    if (vythd[id] > vgrid_max) {
	    cout << "Vgrid_max = " << vgrid_max << "\n";
	    cout << "Vxthd = " << vxthd[id] << "\n";

	      terminate(-1,"Courant violation: thermal vy (vyth) too large");
	    }
	  }
	  if (ndim >= 3) {
	    if (vzthd[id] > vgrid_max) {
	    cout << "Vgrid_max = " << vgrid_max << "\n";
	    cout << "Vxthd = " << vxthd[id] << "\n";

	      terminate(-1,"Courant violation: thermal vz (vzth) too large");
	    }
	  }
	  if (fabs(Bz) > 0) {
	    if (fabs(Ex0_external/Bz) > vgrid_max)
	      terminate(-1,
			"Courant violation: Vdrift (Ex0_external/Bz) too large");
	    if (fabs(Ey0_external/Bz) > vgrid_max)
	      terminate(-1,
			"Courant violation: Vdrift (Ey0_external/Bz) too large");
	  } 
	  if (fabs(Bx) > 0) {
	    if (fabs(Ey0_external/Bx) > vgrid_max)
	      terminate(-1,
			"Courant violation: Vdrift (Ey0_external/Bx) too large");

	  }
	}
      } else {

	FTYPE picPerCell = 	
	  (FTYPE((long)npd[id]*mpi_np)/FTYPE((long)nx*(long)ny)/FTYPE((long)nz*nsubdomains));
        if (unscale_density){
          // n0d is #real particles/sim particles, npd is # sim particles
          // Number of real particles in the entire domain
          printf(";\tPeak Particles in simulation volume: %g\n", 
                 FTYPE(n0d[id]*npd[id]*nsubdomains));
          // Number of real particles/number of simulation particles
          printf(";\tPeak Particle/PIC: %g\n", (FTYPE(n0d[id])));
        }
        else{
          // n0d is the real plasma density, npd is # sim particles
          printf(";\tPeak Particles in simulation volume: %g\n",
                 (FTYPE(n0d[id])*FTYPE(nx*nsubdomains)*dx*ny*dy*nz*dz));
          printf(";\tPeak Particle/PIC: %g\n",
                 (FTYPE(n0d[id])*FTYPE(nx*nsubdomains)*dx*ny*dy*nz*dz)/(FTYPE(npd[id]*mpi_np)));
        }
	printf(";\tSmallest Density Difference +/-(%%): %g\n",
	       100*((picPerCell+1)/picPerCell-1));
	printf(";\tPeak PIC/cell: %g\n",picPerCell);
	printf(";\tPeak Number of pic particles/proc: %d\n",
	       npd[id]);
      }


      
      printf("\n");
      
      
    }
  }
  

  // Check to insure that the initial velocities do not cause a
  // particle to cross multiple grid spacings per time step 
    for (id=0; id<ndist; ++id) {
      if (method[id] >= 0) {
	Vest=(fabs(vx0d[id])+2*fabs(vxthd[id]))*dt/dx;
	if (Vest > 2.0 ) {
	  char string[130];
	  sprintf(string,"%s %g %s\n",
		  "Time Step to large: Particles will cross more than", 
		  Vest, " cells per iteration" );
	  terminate(-1,string);
	}
	Vest=(fabs(vy0d[id])+2*fabs(vythd[id]))*dt/dy;
	if (Vest > 2.0 ) {
	  char string[130];
	  sprintf(string,"%s %g %s\n",
		  "Time Step to large: Particles will cross more than", 
		  Vest, " cells per iteration" );
	  terminate(-1,string);
	}
      }
    }

    if (mpi_rank == 0) printf("\n");

    if (nout <= 0) 
      terminate(-1,"Error: nout must be a positive integer greater than 0");
    for (int ifdist=0; ifdist<fndist; ifdist++) {
      if (fnout[ifdist] <= 0)
	terminate(-1,
		  "Error: fnout must be a positive integer greater than 0");
    }

    
    // Check subcycling parameters
#ifndef CHECK    
    for (id=0; id < ndist; ++id) {
      if (iwrite%subcycle[id] !=0) 
	terminate(-1,"Error: iwrite not an integral multiple of the subcycling variable");
      if (nout%subcycle[id] !=0)
	terminate(-1,"Error: nout not an integral multiple of the subcycling variable");
    }
#endif
    
    // EFIELD INJECT TRIDIAG PARAMS
    if (MAX_TRIDIAG_SOLVE > 0) {
      if (MAX_TRIDIAG_SOLVE%2 != 0) {
	terminate(-1,"Error: max_tridiag_solve must be even or -1");
      }
    }
}
