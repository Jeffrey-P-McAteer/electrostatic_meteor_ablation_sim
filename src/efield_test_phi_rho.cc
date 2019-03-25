
#include "eppic-efield.h"
#include "eppic-mpi.h"
#include "eppic-system.h"
#include <cmath>
#include "eppic-io.h"

void efield_test_phi_rho(field &Efield, FArrayND rho, eppic_system sys, int it)
{

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  FTYPE totalErr=0;
  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

  FTYPE dx2=sys.dx*sys.dx;
  FTYPE dy2=sys.dy*sys.dy;
  FTYPE dz2=sys.dz*sys.dz;

  int xstart=0;
  if (subdomain.id_number==0)
    xstart=0;
  int xend=sys.nx;
  if (subdomain.id_number==nsubdomains-1)
    xend=sys.nx;
  
  FTYPE rho_dc=0;
  if (Efield.zero_dc) {
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<ny;iy++) 
	for(int iz=0;iz<nz;iz++) 
	  rho_dc+=rho(INDICIES(ix,iy,iz));
    rho_dc/=(nx*ny*nz);
    
    if (nsubdomains>1)  {
      FTYPE rho_dc_tmp=rho_dc;
      MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
		    subdomain.neighbor_comm);
      rho_dc/=nsubdomains;
    }
  }
  


  
  for (int ix=xstart;ix<xend;ix++)
    for (int iy=0;iy<sys.ny;iy++)
      for (int iz=0;iz<sys.nz;iz++) {
	int iyp1=(iy+1)%sys.ny;
	int iym1=(iy-1+sys.ny)%sys.ny;
#if NDIM==3
	int izp1=(iz+1)%sys.nz;
	int izm1=(iz-1+sys.nz)%sys.nz;
#endif
	FTYPE rhoxyz = -sys.eps*((phi(INDICIES(ix+1,iy  ,iz  ))+
				  phi(INDICIES(ix-1,iy  ,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dx2
				 +
				 (phi(INDICIES(ix  ,iyp1,iz  ))+
				  phi(INDICIES(ix  ,iym1,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dy2
				 +
				 (phi(INDICIES(ix  ,iy  ,izp1))+
				  phi(INDICIES(ix  ,iy  ,izm1))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dz2
				 );
	FTYPE err = 
	  abs((rho(INDICIES(ix,iy,iz))-rhoxyz-rho_dc)/
	      (rho(INDICIES(ix,iy,iz))-rho_dc));
	totalErr+=err;
      }
  if (totalErr/(sys.ny*sys.nz*(sys.nx-1)) > 1e-8) {
    if (subdomain.rank==0) {
      if (mpi_rank==0)cout << "domain\tavg err\ttotal err\ttime" << endl;
      cout 	 << subdomain.id_number << "\t"
		 << totalErr/(sys.ny*sys.nz*(sys.nx-1))<< "\t"
		 << totalErr<< "\t"
		 << it << "\t"
		 << endl;

      for (int ix=xstart;ix<xend;ix++)
	for (int iy=0;iy<sys.ny;iy++)
	  for (int iz=0;iz<sys.nz;iz++) {
	    int iyp1=(iy+1)%sys.ny;
	    int iym1=(iy-1+sys.ny)%sys.ny;
#if NDIM==3
	    int izp1=(iz+1)%sys.nz;
	    int izm1=(iz-1+sys.nz)%sys.nz;
#endif
	    FTYPE rhoxyz = -sys.eps*((phi(INDICIES(ix+1,iy  ,iz  ))+
				  phi(INDICIES(ix-1,iy  ,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dx2
				 +
				 (phi(INDICIES(ix  ,iyp1,iz  ))+
				  phi(INDICIES(ix  ,iym1,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dy2
				 +
				 (phi(INDICIES(ix  ,iy  ,izp1))+
				  phi(INDICIES(ix  ,iy  ,izm1))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dz2
				 );
	    FTYPE err = 
	      abs((rho(INDICIES(ix,iy,iz))-rhoxyz-rho_dc)/
		  (rho(INDICIES(ix,iy,iz))-rho_dc));
	    if (err > 1e-8) 
	      cout 	 << ix+sys.nx*subdomain.id_number << "\t"
			 << iy << "\t"
			 << iz << "\t"
			 << rho(INDICIES(ix,iy,iz)) << "\t"
			 << rhoxyz << "\t"
			 << rho_dc << "\t"
			 << rho(INDICIES(ix,iy,iz)) - rhoxyz-rho_dc<< "\t"
			 << endl;
	  }

      terminate(-1,"\n");    
    }

  }
}
