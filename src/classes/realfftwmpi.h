// realfftwmpi<Type, A_NDIM>:  A simple A_NDIM array class - 
//                       designed to be entirely inline.

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)

#ifndef realfftwmpi_H 
#define realfftwmpi_H

#include "ArrayNd_ranged.h"

// We do not want MPI C++ bindings
#define MPI_NO_CPPBIND
#define MPICH_SKIP_MPICXX
#define MPICH_IGNORE_CXX_SEEK

#include <rfftw_mpi.h>

const index_type nx_guard_default0[]={0,0,0};


template <class TYPE,int A_NDIM>
class realfftwmpi 
{

 public:

  ArrayNd_ranged<TYPE,A_NDIM> data;
  rfftwnd_mpi_plan plan, iplan;
  fftw_complex *cdata;

  int x_start; // x start position in the global array
  int ny_transpose, y_start_transpose; // Local transposed values
  int local_size; //Total local size
  int local_nx; // nx on the local processor

 protected:
    FTYPE *workspace;

 public:
  //Constructors:
    //Default Constructor:

  realfftwmpi(): 
    ny_transpose(0), y_start_transpose(0),
    local_size(0), plan(0), iplan(0), cdata(0) {};

  //Copy Constructor:
  /*  realfftwmpi(const realfftwmpi &in) 
    data=in.data;
    x_start=in.x_start; 
    ny_transpose=0;
    y_start_transpose=0;
    cdata=in.cdata;
    plan=in.plan;
    iplan=in.plan;
    local_size=in.local_size;
{}
*/

  //Constructor with size and domain defined:
  realfftwmpi(const MPI_Comm MPI_COMM_DOMAINS, const index_type n_in[], 
	      const index_type nx_guard[]=nx_guard_default0, TYPE x=0) 
    {

      // An fftw plan must be established to determine the sizes in each domain:
      int j;

      //      cout << mpi_rank << " A_NDIM= " << A_NDIM << " sizes: ";
      //      for (int i=0;i<A_NDIM;i++) cout << n_in[i] << " " ;
      //      cout << "\n";
    int fftw_plan_type = FFTW_MEASURE;
#ifdef DEBUG 
   fftw_plan_type = FFTW_ESTIMATE;
#endif
      if (A_NDIM == 3) {

	  plan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1], n_in[2],
					 FFTW_REAL_TO_COMPLEX, );

	  iplan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1], n_in[2],
					  FFTW_COMPLEX_TO_REAL, fftw_plan_type);

      } else if(A_NDIM == 2) {

	  plan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
					 FFTW_REAL_TO_COMPLEX,fftw_plan_type);

	  iplan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
					  FFTW_COMPLEX_TO_REAL,	fftw_plan_type);
      } else {
	cerr << "realfftwmpi only 2D or 3D implemented"<<"\n" ;
	exit(-1);
      }

      //      cout << mpi_rank << " plan done" <<"\n" ;

      //This routine returns the sizes in each domain:

      int fftw_local_size=0;
      rfftwnd_mpi_local_sizes(plan, &local_nx, &x_start,
			      &ny_transpose, &y_start_transpose,
			      &fftw_local_size);

      //      cout << mpi_rank << " local_nx=" <<  local_nx << " x_start=" <<  x_start 
      //      	   << " ny_transpose=" << ny_transpose << " y_start_transpose=" 
      //	   <<  y_start_transpose << " fftw_local_size=" << fftw_local_size << "\n";
      

      // Note: the local size is padded in the first and last dimensions.
      // It is padded in the last dim to make space for the complex transform
      // and in the first dim for guard cells.
      // A 3D array will have ranges of array(0-nx_guard[0]..nx-1+nx_guard[1],0..ny-1,-0.nz+1).

      index_type nsize[A_NDIM];
      nsize[0]=local_nx+nx_guard[0]+nx_guard[1];
      for (index_type i=1; i<A_NDIM-1; i++) nsize[i]=n_in[i];
      //Note: This matrix pads the final dimension of the array as required 
      //      by the fftw real to complex transforms so that for 3-D, 
      //      nz = 2*(n_in[2]/2+1)
      nsize[A_NDIM-1]=2*(n_in[A_NDIM-1]/2+1);

      index_type memory_after_0=local_nx+nx_guard[1];
      int nynz=nsize[1];
      if (A_NDIM == 3) nynz *= nsize[A_NDIM-1];
      memory_after_0 *= nynz;

      // Set the starting index for the array:
      int  nstart[]={0,0,0};
      nstart[0]=-nx_guard[0];

      if (memory_after_0 >= fftw_local_size) {
	workspace=0;
	local_size=memory_after_0+nx_guard[0]*nynz;
//	data=ArrayNd_ranged<TYPE,A_NDIM>(nstart,nsize,x);
	data.build_ArrayNd(nstart, nsize, x);
      } else {
//	cout << "\nWarning from realfftwmpi: fftwmpi requesting more memory on proc "
//	     << mpi_rank << " than fits.  Padding end of array\n" ;
	local_size=fftw_local_size+nx_guard[0]*nynz;
	workspace=new TYPE[local_size];
	data.build_ArrayNd(nstart, nsize, x, workspace);
      } 

      //      cout << mpi_rank << " local_size=" << local_size << " nsize=("; 
      //      for (index_type i=0; i<A_NDIM; i++) cout << " "<<nsize[i]; 
      //      cout << ") nstart[0]=(";
      //      for (index_type i=0; i<A_NDIM; i++) cout << " "<< nstart[i]; 
      //      cout << ")\n";
      
    };

  ~realfftwmpi()
    {
      if (workspace !=0){
	delete [] workspace;
	workspace=0;
      }
    }

  //  Access and setting -- Speed is everything
  TYPE& operator() (int n0){ return data(n0);}
  TYPE& operator[] (int n0){ return data[n0];}
  TYPE& operator() (int n0, int n1){return data(n0,n1);}
  TYPE& operator() (int n0, int n1, int n2){return data(n0,n1,n2);}

  //Fourier Transform Array:
  void transform(TYPE* work=0) {
    
      bool work_defined = false;
      if (work == 0) {
	  work=new TYPE[local_size];
	  work_defined = true;
      }

    int index_zero[3]={0,0,0};

    /* Now, compute the forward transform: */
    TYPE* start_add=data.address(index_zero);
    rfftwnd_mpi(plan, 1, start_add, work, FFTW_NORMAL_ORDER);

    /* the data is now complex, so typecast a pointer: */
     cdata = (fftw_complex*) start_add;   

     if (work_defined) delete [] work;
     
  }

  //Inverse Fourier Transform Array:
  void invtransform(TYPE *work=0) {
      
      bool work_defined = false;
      if (work == 0) {
	  work=new TYPE[local_size];
	  work_defined = true;
      }
    int index_zero[3]={0,0,0};
    /* Now, compute the forward transform: */
    rfftwnd_mpi(iplan, 1,  data.address(index_zero), work, FFTW_NORMAL_ORDER);
    
    if (work_defined) delete [] work;
    
  }

};

#endif
