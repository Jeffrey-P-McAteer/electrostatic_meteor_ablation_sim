#include <stdio.h> 
#include <string> 
#include <math.h>
#include <iostream>
#include "eppic.h"
#include "eppic-mpi.h"

void output_energy(FILE *fcons, particle_dist *pic, field &Efield){
  // Calculate and output energies 
  //Local energy variables
  FTYPE Wfield;
  FTYPE netnp = 0;
  FTYPE energy(field &Efield, particle_dist *pic, FTYPE &Wfield, 
	       FTYPEAVec &Wpart);
  FTYPEAVec Wpart(ndist); // Particle Energy DENSITY

  energy(Efield, pic, Wfield, Wpart);
  
  // fluid enregy -- to be done later

  if (isnan(Wpart.sum()+Wfield)) {
    
    cout << "Proc(" << mpi_rank << "): "
	 << "Wpart.sum = " << Wpart.sum()
	 << " Wfield = " << Wfield << "\n";
    
    if (isnan(Wpart.sum())) {
      for (int id=0; id<ndist; ++id) {
	if (method[id] >= 0 ) {
	  if (method[id] == 0) {
	    for (int i=0;i<pic[id].np;i++) {
	      FTYPE picE = Sqr(pic[id].vx(i)*dx/dt);
	      if (vel_dim[id] >= 2) picE += Sqr(pic[id].vy(i)*dy/dt);
	      if (vel_dim[id] == 3) picE += Sqr(pic[id].vz(i)*dz/dt);
	      if (isnan(picE)) {
		cout << "particle " 
		     << i
		     << " has NAN energy: "
		     << endl
		     << " x = " << pic[id].x(i) 
#if NDIM>1		 
		     << " y = " << pic[id].y(i)
#endif
#if NDIM>2     
		     << " z = " << pic[id].z(i)
#endif		  
		     << endl;
		cout << " vx = " << pic[id].vx(i);
		if (vel_dim[id] > 1)
		  cout << " vy = " << pic[id].vy(i);
		if (vel_dim[id] > 2)
		  cout << " vz = " << pic[id].vz(i);
		cout << endl;
	      }	      
	    }
	  }
	}
      }
    }
    terminate(-1,"Error: Infinite energy");
  }

  
  for (int i =0; i<ndist; i++)
    {
      netnp += pic[i].np - pic[i].nabsent;
    }
#ifdef USE_MPI
    {
      FTYPEAVec Wpart_all(ndist);  
      Wpart_all=0;
      int mpi_err=MPI_Reduce(Wpart.address(),Wpart_all.address(), 
			    ndist,MPI_FTYPE,MPI_SUM,0,MPI_COMM_WORLD);
      if ( mpi_err != MPI_SUCCESS) 
	mpi_error(mpi_rank, mpi_err,"Reduce call in output");
     Wpart = Wpart_all;

     // Output the number of particles
     FTYPE netnp_global;
     FTYPE netnp_subdomain;
     mpi_err=MPI_Reduce(&netnp, &netnp_subdomain,
                        1, MPI_FTYPE, MPI_SUM,0,subdomain.internal_comm);
     // Get average number of particles inside subdomain processors
     netnp = netnp_subdomain/subdomain.np;
     // Add up across entire global gomain
     mpi_err=MPI_Reduce(&netnp, &netnp_global,
                        1, MPI_FTYPE, MPI_SUM,0,subdomain.neighbor_comm);
     /*
     for (int i =0; i<ndist; i++)
       {
	 netnp += pic[i].np_all;
       }
     mpi_err=MPI_Reduce(&netnp, &netnp_global,
			    1, MPI_FTYPE, MPI_SUM,0,MPI_COMM_WORLD);

     */
     netnp = netnp_global;
    }
#endif

#ifdef USE_DOMAINS
    {
      FTYPE Wfield_all;
      int mpi_err=MPI_Reduce( &Wfield, &Wfield_all, 
			     1, MPI_FTYPE,MPI_SUM,0,subdomain.neighbor_comm);
      if ( mpi_err != MPI_SUCCESS) 
	mpi_error(mpi_rank, mpi_err,"Reduce call in output_energy");
      Wfield = Wfield_all;
    }
#endif    


  
  if (mpi_rank == 0) {
    
    printf(" %10.5g %10.5g %10.5g", Wpart.sum(), Wfield, netnp);
    fprintf(fcons," %12.5g %12.5g\n", Wpart.sum(), Wfield);
  }    
}
