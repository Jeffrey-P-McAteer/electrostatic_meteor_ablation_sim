// A test program for the efield routine found in the libeppic.a library
// First we initialize all the variables need to call efield.
// Then we generate an initial periodic density on a grid.
// We calculate the dc component to rho, so that a comparison can be made.
// We use efield to calculate phi (Efield.phi_rho).
// Then we numerically differentiate phi to compare with rho. 
// We compare the two results and pass if they are the same within a tolerance.


using namespace std;
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "ArrayNd.h"
#include "eppic.h"
#include "eppic-efield.h"
#include "eppic-mpi.h"



clock_t first_time=0;

#include "global_defs.h"
int main(int argc, char* argv[])
{
  // initialize functions called
  void efield(field &Efield, FArrayND &rho);

  FArrayND rho;
  FArrayND rho_start;
  field Efield;

  // setup MPI
  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if (mpi_rank == 0) {
    MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
  }

  // initialize variables
  // system defining parameters
#ifdef DEBUG  

  // initialize variables
  // system defining parameters
  nx = 48;
  ny = 24;
#if NDIM==2
  nz = 1;
#else 
  nz = 24;
#endif


  int Nit=1; /* Number of times to run efield (for scaling only) */


  if (mpi_np > 1) {
    if (mpi_np > 4) 
      nsubdomains = 4;
    else 
      nsubdomains = mpi_np;
      nx/=nsubdomains;
  }

  if (mpi_rank==0)
    cout << "NSUBDOMAINS = " << nsubdomains << endl;

  //int boundarys[2]={neumann,open}; // natural, ie open
  bc_type boundarys[2]={periodic,periodic}; // natural, ie open
#else


  if (mpi_rank==0) {
    if (argc<7) {
      cout << "Not enough arguments" <<endl;
      cout << "NX NY NIT NSUBDOMAINS LHS_BOUNDARY RHS_BOUNDARY" << endl;
      MPI_Finalize();
      return 1;
    }
  }
  // initialize variables
  // SYSTEM defining parameters
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  ny = atoi(argv[3]);
  int Nit=atoi(argv[4]); /* Number of times to run efield (for scaling) */
  
  nsubdomains = atoi(argv[5]);  
  bc_type boundarys[2]={static_cast<bc_type>(atoi(argv[6])),
			static_cast<bc_type>(atoi(argv[7]))}; // natural, ie open
  
#endif
  dx=.9,dy=.9,dz=.9;
  eps=1;//8.8542e-12;

  FTYPE boundary_vals[4] = {0,0,0,0};
  init_field(Efield,rho,nx,ny,nz,dx,dy,dx,eps,boundarys,boundary_vals);
  domain_decomp();

  // set rho to test case
  // take care of NDIM specific loops
  rho=0.0;
  for(int ix=0;ix<nx;ix++) {
    for(int iy=0; iy<ny;iy++) {
      for(int iz=0; iz<nz;iz++) {
      rho(ix,iy) = 
	cos(2*PI*1*(ix+mpi_rank*nx)/FTYPE(nx*nsubdomains))*
	sin(2*PI*1*iy/FTYPE(ny))*
	sin(2*PI*1*iz/FTYPE(nz));
	//sin(ix/FTYPE(nx))*sin(iy/FTYPE(ny));
      }
    }
  }

  //rho dc part calculation
  FTYPE rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) rho_dc+=rho(ix,iy);
  rho_dc/=(nx*ny);
  FTYPE rho_dc_tmp=rho_dc;
  MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
		subdomain.neighbor_comm);
  rho_dc/=nsubdomains;
  
  rho_start=rho;
  for (int it=0;it<Nit;it++){
    rho = rho_start;
    //get phi
    Efield.phi_rho = 0.0;
    efield(Efield,rho);
    
    if (mpi_rank == 0) 
      if (it==0)
	first_time = efield_time;


    /* binary write to file */

    if (it==0) {
      void write_domains(FArrayND_ranged &array,MPI_Comm across_comm,
			 char *path_name,int nghost_pts[2],
			 bool non_periodic_x);

      FArrayND_ranged rho_ranged;

      if (subdomain.rank==0) {
	int non_periodic = 1;
	if ((Efield.boundary[0]==periodic)||(Efield.boundary[1]==periodic))
	  non_periodic = 0;
	int nghost_phi[2]={1,2};
	write_domains(Efield.phi_rho,
		      MPI_COMM_WORLD,"efield_test_phi*.out",
		      nghost_phi,non_periodic);
	
	int nghost_rho[2]={0,1};
	//int nstart_rho[]={INDICIES(0,0,0)};
	//rho_ranged= FArrayND_ranged(nstart_rho,nghost_rho);
	rho_ranged=rho_start;
	write_domains(rho_ranged,
		      MPI_COMM_WORLD,"efield_test_rho*.out",
		      nghost_rho,0);
      }
    }

#ifdef DEBUG
    int failed=0;

    //differentiate phi
    char errormsg[90];
    int ixm1=-1;
    int ixp1=0;
    int iz=0;
    FTYPE errtot=0.;
    if (subdomain.id_number==0) {
    sprintf(errormsg,"%10s %4s %4s %16s %16s %16s %16s\n",
	    " ","ix","iy","deviation","relative error","rho_calc","rho_orig");
    cout << errormsg;
    }
    for (int ix=0; ix<nx; ix++) {
      ixp1=ix+1;
      int iym1=ny-1;
      for (int iy=0; iy<ny; iy++) {
	int izm1=nz-1;
	for (int iz=0; iz<nz; iz++) {

	FTYPE rho_test,rho_now;
	
	rho_now = rho_start(INDICIES(ix,iy,iz));
	rho_test=0;
	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dx));
	rho_test += Efield.phi_rho(INDICIES(ix+1,iy,iz))/Sqr(dx);
	rho_test += Efield.phi_rho(INDICIES(ixm1,iy,iz))/Sqr(dx);
	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dy));
	rho_test += Efield.phi_rho(INDICIES(ix,(iy+1)%ny,iz))/Sqr(dy);
	rho_test += Efield.phi_rho(INDICIES(ix,iym1,iz))/Sqr(dy);
#if NDIM == 3
	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dz));
	rho_test += Efield.phi_rho(INDICIES(ix,iy,(iz+1)%nz))/Sqr(dz);
	rho_test += Efield.phi_rho(INDICIES(ix,iy,izm1))/Sqr(dz);
#endif	
	rho_test *= -1.0*eps;
	rho_test += rho_dc;
	

	
	if (fabs(rho_now) > 1e-6) {
	  if (fabs((rho_test-rho_now)/rho_now)>1e-6){
	    failed++;
	    if (subdomain.rank == 0) {
	      sprintf(errormsg,"%10s %4d %4d %16g %16g %16g %16g\n",
		      "error is",
		      ix+subdomain.id_number*nx,iy,
		      fabs((rho_test-rho_now)),
		      fabs((rho_test-rho_now)/rho_now)
		      ,rho_test,rho_now);
	      cout << errormsg;
	    }
	  }
	}
	else
	  if (fabs((rho_test-rho_now))>1e-6){
	    failed++;
	    sprintf(errormsg,"%10s %4d %4d %16g %16g %16g %16g\n",
		    "error is",
		    ix+subdomain.id_number*nx,iy,
		    fabs((rho_test-rho_now)),
		    fabs((rho_test-rho_now)/rho_now)
		    ,rho_test,rho_now);
	    cout << errormsg;
	  }
	
	
	errtot +=fabs(rho_test-rho_now);

	izm1=iz;
	}
	
	iym1=iy;
      }
      ixm1=ix;
    }
    if (failed) {
      MPI_Finalize();
      return 1;
    }
#endif
  }
#ifndef DEBUG
  /* total time, average time, total less first, average less first */
  if (mpi_rank == 0) cout <<  nx << "\t";
  if (mpi_rank == 0) cout <<  ny << "\t";
  if (mpi_rank == 0) cout <<  mpi_np << "\t";
  if (mpi_rank == 0) cout <<  nsubdomains << "\t";
  if (mpi_rank == 0) cout <<  Nit << "\t";
  if (mpi_rank == 0) cout <<  boundarys[0] << "\t";
  if (mpi_rank == 0) cout <<  boundarys[1] << "\t";
  if (mpi_rank == 0) cout <<  efield_time << "\t";
  if (mpi_rank == 0) cout <<  efield_time/Nit<< "\t";
  if (mpi_rank == 0) cout <<  efield_time - first_time<< "\t";
  if (mpi_rank == 0) cout << (efield_time - first_time)/(Nit-1) << endl;
#endif

  MPI_Finalize();
  return 0;
}




// int main(int argc, char* argv[]){
//   cout << "3D test not implemented!" << endl;
//   return 1;

// }



