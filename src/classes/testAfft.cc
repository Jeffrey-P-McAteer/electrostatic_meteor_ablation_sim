// Global variables, functions, types and structures

#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <locale.h>
#include <signal.h>
#include <math.h>

#include <fstream.h>
#include <string>
#include <sstream>

#undef NDIM
#define NDIM 3
#define FTYPE double


#include <iostream>
#include <sstream>

#ifdef CHECK
#define DEBUG
#else
#ifndef NOCHECK
#define NOCHECK
#endif
#endif

#if NDIM == 1
#define INDICIES(ix,iy,iz) ix
#elif NDIM == 2
#define INDICIES(ix,iy,iz) ix, iy
#else
#define INDICIES(ix,iy,iz) ix, iy, iz
#endif

//Parallel processing variables will be known globally:
int mpi_rank=0, mpi_np=1;

inline double Sqr(double x) {return x*x;}
inline float Sqr(float x) {return x*x;}

#ifdef USE_MPI
#include "ArrayNd_fftwmpi.h"
typedef ArrayNd_fftwmpi<FTYPE, NDIM> Arrayfft ; 
#else
#include "ArrayNd_fftw2.h"
typedef ArrayNd_fftw2<FTYPE, NDIM> Arrayfft;
#endif

#include "ArrayNd_ranged.h"


template <class T> std::string ConvertToString(const T& val)
{
std::stringstream ost;
ost << val;
return ost.str();
}

#ifdef USE_MPI
void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg) {
  int stringlen=MPI_MAX_ERROR_STRING;
  char string[MPI_MAX_ERROR_STRING];
  printf("ERROR in processor %d: MPI error %d from %s\n",
	 mpi_rank,mpi_err,mpi_msg);
  fprintf(stderr,"ERROR in processor %d: MPI error %d from %s\n",
	  mpi_rank,mpi_err,mpi_msg);
  MPI_Error_string(mpi_err, string, &stringlen);
  printf("                             : %s\n", string);

  char message[]="ArrayNd_fftwmpi termination";
  
  printf("\n%s\n\nPPIC3D process %d terminating all %d processes\n",
	 message, mpi_rank, mpi_np);
  fprintf(stderr,"\n%s\n\nPPIC3D process %d terminating all %d processes\n",
	  message, mpi_rank, mpi_np);
  
#ifdef CHECK
  // Raise a trap signal to allow a debugger to catch this problem.
  raise(SIGTRAP);
#endif
  
  MPI_Abort(MPI_COMM_WORLD,-1);
  MPI_Finalize();

  exit(-1);

}
#endif

main(int argc, char* argv[])
{
  
  int mpi_err;

#ifdef USE_MPI
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
    printf("\nMPI STARTING %d PROCESSORS\n MPI Error Handler: %d\n\n", 
	   mpi_np, mpi_err_hand);
  }

  MPI_Comm MPI_COMM_DOMAINS=MPI_COMM_WORLD;
#endif



  const index_type nx=8, ny=4, nz=4;
  const index_type nsize[3]={nx*mpi_np, ny, nz};

  // Define an array:

#ifdef USE_MPI
  Arrayfft phi(MPI_COMM_DOMAINS,nsize);
#else
  Arrayfft phi(nsize);
#endif


 /* initialize data to a simple wave function: */
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz) {
	phi(INDICIES(ix, iy, iz))
	  = sin((ix + phi.x_start)/(nx/(2.*M_PI)));//*cos(iy/(nx/2/M_PI));
      }

  std::string filename="phi_initial"+ ConvertToString(mpi_rank) + ".txt";
  phi.output(filename.c_str());
	
  /* Now, compute the forward transform: */
  phi.transform();

  filename = "phi_trans"+ ConvertToString(mpi_rank) + ".txt";
  phi.output(filename.c_str());

  //  cout << "Transform Complete...printing:" <<endl ;

  /* transformed phi will be multiplied by a ksqrinv array: 
     generate this here: */ 

  Arrayfft ksqrinv(
#ifdef USE_MPI
MPI_COMM_DOMAINS, 
#endif
nsize);

  // ksqrinv is just the fft of the finite difference opperator:
  double dx=1., dy=1., dz=1.;

  cout << "setting ksqrinv" << mpi_rank <<" "<< phi.x_start <<" "<< nx <<" "<< nsize[0] <<" "<< (phi.x_start + nx == nsize[0]) <<endl;
  if (phi.x_start + nx == nsize[0] && nx > 0) // Only the last processor:
    ksqrinv(INDICIES(nx-1,0,0))+= 1/Sqr(dx);
  if (phi.x_start == 0 && nx > 0) {  // Only the first processor
    ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dx);
    ksqrinv(INDICIES(1%nx,0,0))+= 1/Sqr(dx);
#if NDIM >1
    ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dy);
    ksqrinv(INDICIES(0,1%ny,0))+= 1/Sqr(dy);
    ksqrinv(INDICIES(0,ny-1,0))+= 1/Sqr(dy);
#endif
#if NDIM == 3
    ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dz);
    ksqrinv(INDICIES(0,0,1%nz))+= 1/Sqr(dz);
    ksqrinv(INDICIES(0,0,nz-1))+= 1/Sqr(dz);
#endif
  }
  filename = "ksqrinv0_" + ConvertToString(mpi_rank) + ".txt";
  ksqrinv.output(filename.c_str());

  ksqrinv.transform();

  filename="ksqrinv_trans"+ ConvertToString(mpi_rank) + ".txt";
  ksqrinv.output(filename.c_str());

  // ksqrinv is now a k^2 like operator. Invert, element by element
  // and make the complex = real part for rapid multiplication later.
  for (int ik=2; ik<ksqrinv.length(); ik += 2) {
    ksqrinv(ik) =1./ksqrinv(ik);
    ksqrinv(ik+1)=ksqrinv(ik);
  }
  
  // Incorporate the epsilon coefficient and the n_elements factor into
  // ksqrinv:
  
  const double eps = 1.0;
  ksqrinv /= -eps * (mpi_np*nx)*ny;
  if (NDIM == 3) ksqrinv /= nz;

  filename="ksqrinv"+ ConvertToString(mpi_rank) + ".txt";
  ksqrinv.output(filename.c_str());

  /* By multiplying fft of the charge density (stored in rho) by ksqrinv, 
     we obtain the fft of phi */

  phi *= ksqrinv;

  /* Now, compute the inverse transform to obtain phi: */
  phi.invtransform();

  filename="phi_final"+ ConvertToString(mpi_rank) + ".txt";
  phi.output(filename.c_str());

  //Test copy constructor
  ksqrinv=phi;
  //  cout << ksqrinv;
    

  // Let's now test the ranged array type:

  const int astart[]={-1,0,0};
  const int arange[]={nx+2,ny,nz};
  ArrayNd_ranged<double, NDIM> a(astart,arange,1), b;

  // This does not work because it goes into the destructor???
  int asize[NDIM][2]={{-1,nx},{0,ny-1},{0,nz-1}};
  ArrayNd_ranged<double, NDIM> aa(asize,1), bb;
  
  b=a;

  a=2.;
  cout << " a(-1,0,0)=" << a(-1,0,0) 
       << "\t a(nx,ny-1,nz-1)="<< a(nx,ny-1,nz-1) 
       << "\t a(static_cast<int>(nx),0,0)=" << a(static_cast<int>(nx),0,0) 
       << endl << flush;
  a /= 0.5;
  a(-1,0,0)=1;
  cout << " a(-1,0,0)="<< a(-1,0,0) 
       << " a(nx,ny-1,nz-1)=" << a(nx,ny-1,nz-1) 
       << " a(static_cast<int>(nx),0,0)=" << a(static_cast<int>(nx),0,0) 
       << endl << flush;

  a *= 3.33333;

  cout <<"a " <<  a <<endl << flush;


  if (a(-1,0,0) != 1.0*3.33333) 
    cout << a(-1,0,0) <<  "Error: a(-1,0,0) != 1.0*3.33333:  ArrayND_ranged not working\n";

  if (a(nx,ny-1,nz-1) != 2.0/0.5*3.33333) 
    cout << a(nx,ny-1,nz-1) << "Error: a(nx,ny-1,nz-1) != 2.0/0.5*3.33333:  ArrayND_ranged not working\n";


  aa=2.;
  aa-=3;
  bb=aa;
  bb(-1,1,1)=2;
  cout << "bb " << bb;

  //  bb(-2,1,1)=2;
  


#ifdef USE_MPI  
  MPI_Finalize();
#endif
  {// Watch deconstructor
  int asize2[NDIM][2]={{-1,nx},{-3,ny-1},{-4,nz-1}};
    ArrayNd_ranged<double, NDIM> a2(asize2,1);
    a2 *= 4.;
    a2(-1,ny-1,-4)=6.0;
    a2(-1,-3,-4)=-1.0;
    a2(nx,ny-1,nz-1)=8.0;
    cout <<  a2;
      
    
    }
  //  exit(-1);

}
