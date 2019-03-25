  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     A test program for pic_current routine. 

     Notes:
     Currently the frame work is setup to test this routine. The PIC 
     distribution is too ideal to consider this test trust worthy, and 
     it should be adjusted. See my_init_particles (found below) to 
     control how the distributions are setup.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

using namespace std;
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "ArrayNd.h"
#include "eppic.h"
#include "eppic-mpi.h"
#include "global_defs.h"

// pic specific 
int ndist_ion=1;
FTYPE qion=QI;

FTYPE (*vFunc)(INDICIES(int ix, int iy, int iz));

int main(int argc, char* argv[])
{
  nx = 40;
  ny = 40;
  nz = 1;
  
  int failed = 0;
  particle_dist *pic;
  FArrayND_ranged J;
  FArrayND Jtmp;
  const int nsize_J[]={INDICIES(nx+1+1,ny,nz)}; 
  const int nsize_Jtmp[]={INDICIES(nx,ny,nz)}; 
  int nstart_J[]={INDICIES(-1,0,0)};
  int testdim=0; /* configure using command-line */
  int vFuncType=-1;

  /* function declarations */
  void my_init_particles(particle_dist *&pic,
			 FTYPE (*vFunc)(INDICIES(int ix, int iy, int iz)));
  void pic_flux(particle_dist &pic,int dim,FArrayND &fluxDim,
		FTYPE scaler, int n_avg);

  void domain_decomp();
  int test(FTYPE J,INDICIES(int ix, int iy, int iz));
  // possible velocity functions
  FTYPE vel_one(INDICIES(int ix, int iy, int iz));
  FTYPE vel_sine(INDICIES(int ix, int iy, int iz));

  // setup MPI
  int mpi_err,numprocs;
  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

  MPI_Init(&argc,&argv);
  mpi_err=MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"rank call in main");
  mpi_err=MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"size call in main");
  if (mpi_rank == 0) {
    mpi_err=MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
    //printf("\nMPI STARTING %d PROCESSORS\n MPI Error Handler: %d\n\n", 
    //mpi_np, mpi_err_hand);
  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize arguments
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize variables
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ndist += ndist_ion;
  method=intAVec(ndist) = 0;  
  method[0]=-1; // make first dist electrons
  J = FArrayND_ranged(nstart_J, nsize_J);
  J=0.0;
  Jtmp = FArrayND(nsize_J);
  Jtmp = 0.0;

  if (mpi_np > 1) {
    nsubdomains = mpi_np;
    nx/=nsubdomains;
  }

  domain_decomp();
  

  switch(vFuncType)
    {
    default:
      vFunc = &vel_one;
    }

  vFunc = &vel_sine;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize pic with ideal density and velocity
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  my_init_particles(pic,vFunc);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Call pic_current
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (int id=0;id<ndist;id++) {
    if (method[id]>=0) {
      pic_flux(pic[id],testdim,Jtmp,pic[id].q,1);
      for(int ix=0;ix<nx;ix++)
	for (int iy=0;iy<ny;iy++)
	  for (int iz=0;iz<nz;iz++)
	    J(INDICIES(ix,iy,iz)) += Jtmp(INDICIES(ix,iy,iz));
      }
  }

  if (nsubdomains>1) {
    void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
    const int nx_guard[] = {1,1};
    pass_guards(J,nx_guard);
  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Test value of current returned
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  // test velocity at each mesh point
  //cout << "....... proc: " << mpi_rank << " Starting mesh point test\n";
  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      for (int iz=0;iz<nz;iz++){
	failed+= test(J(INDICIES(ix,iy,iz)),INDICIES(ix,iy,iz));
      }
    }
  }

  // test the boundaries
  //cout << "....... proc: " << mpi_rank << " Starting boundary point test\n";
  if (nsubdomains>1) {
    int ix,iz=0;
    
    for(int iy=0;iy<ny;iy++){
      ix=-1;
      failed+= test(J(INDICIES(ix,iy,iz)),INDICIES(ix,iy,iz));
      ix=nx;
      failed+= test(J(INDICIES(ix,iy,iz)),INDICIES(ix,iy,iz));
    }
    
  }
  
  MPI_Finalize();
  //cout << "....... proc: " << mpi_rank << " Failed: " << failed << "\n";
  if (failed) return 1;
  else return 0;
}


void my_init_particles(particle_dist *&pic,
		       FTYPE (*vFunc)(INDICIES(int ix, int iy, int iz))) 
{
  
  /* initializing particles using code striped from init_particles.cc */
  pic=new particle_dist [ndist];
  for (int id=0;id<ndist;id++) {
    if (method[id]>=0) {
      /* using cpp defaults defined in eppic.h, or infile.cc */
      pic[id].m=MI; 
      pic[id].q=qion;
      pic[id].n0avg=1;//1.e6;
      pic[id].np=nx*ny*nz;
      int npad = (int)(pic[id].np*1.1); 
      pic[id].x=PTYPEAVec(npad) = 0.;
      if (ndim >= 2) pic[id].y=PTYPEAVec(npad) =0.; 
      if (ndim == 3) pic[id].z=PTYPEAVec(npad) = 0.;
      
      if (nsubdomains > 1 || boundary_type == INJECT) {
	pic[id].absent=intAVec((int) (pic[id].np*.1))=-1;
	pic[id].nabsent=0;
      }
      pic[id].vx=PTYPEAVec(npad) = 0.;
      if (ndim >= 2) pic[id].vy=PTYPEAVec(npad) = 0.;
      if (ndim == 3) pic[id].vz=PTYPEAVec(npad) = 0.;
      // Should be one particle per cell, so place evenly
      int part_id=0;
      for (int ix=0;ix<nx;ix++){
	for(int iy=0;iy<ny;iy++) {
	  for(int iz=0;iz<nz;iz++) {
	    pic[id].x(part_id) = ix;//ix+dx/2.;
	    if (ndim >=2 ) pic[id].y(part_id) = iy;//iy+dy/2.;
	    if (ndim ==3 ) pic[id].z(part_id) = iz;//iz+dz/2.;
	    part_id++;
	  }
	}
      }
      for (int part_id=0;part_id<pic[id].np;part_id++){
	pic[id].vx[part_id] = vFunc(INDICIES(static_cast<int>(pic[id].x[part_id]),static_cast<int>(pic[id].y[part_id]),static_cast<int>(pic[id].z[part_id])));
      }
    }
  }
}



FTYPE vel_one(INDICIES(int ix, int iy, int iz))
{
  return FTYPE(1.0);
}

FTYPE vel_sine(INDICIES(int ix, int iy, int iz))
{
  int ix_local = ix-subdomain.id_number*nx;
  if (ix_local<0) ix_local += nx*nsubdomains;
  ix_local = ix_local%(nx*nsubdomains);
  FTYPE vel = sin(ix_local/(2*nx*dx*nsubdomains));
  if (ndim>=2) vel*=sin(iy/(2*ny*dy));
  return vel;
}


int test(FTYPE J,INDICIES(int ix, int iy, int iz))
{

  if (fabs(J-qion*vFunc(INDICIES(ix,iy,iz)))>1.e-8) {
    ix+= nx*subdomain.id_number;
    cout << "Rank = " 
	 << mpi_rank
	 << "FAILED AT " 
	 << ix << ","
	 << iy << ","
	 << " value = " 
	 << J
	 << " value v = " 
	 << qion*vFunc(INDICIES(ix,iy,iz))
	 << " diff = " 
	 << (J-qion*vFunc(INDICIES(ix,iy,iz)))
	 << "\n";
    return 1;
  }
  return 0;
}

