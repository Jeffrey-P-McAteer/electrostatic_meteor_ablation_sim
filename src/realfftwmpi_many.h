// realfftwmpi_many<Type, A_NDIM>:  A simple A_NDIM array class - 
//                       designed to be entirely inline.

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)

#ifndef realfftwmpi_many_H 
#define realfftwmpi_many_H

#include "ArrayNd_ranged.h"

// We do not want MPI C++ bindings
#define MPI_NO_CPPBIND
#define MPICH_SKIP_MPICXX
#define MPICH_IGNORE_CXX_SEEK

#include "eppic_fftw_mpi.h"
#include <complex>

template <class TYPE,int A_NDIM>
class realfftwmpi_many 
{

 public:

  ArrayNd_ranged<TYPE,A_NDIM+1> data;
  rfftwnd_mpi_plan plan, iplan;
  rfftwnd_plan plan_nompi, iplan_nompi;
  ArrayNd_ranged<std::complex<TYPE>,A_NDIM+1> cdata;
  ArrayNd_ranged<TYPE,A_NDIM+1> data_transposed; //Container for transposed data

  int x_start; // x start position in the global array
  int ny_transpose, y_start_transpose; // Local transposed values
  int local_size; //Total local size
  int local_nx; // nx on the local processor
  int howmany; // how many transforms to perform 
  index_type nsize[4];
  MPI_Comm comm;

  enum order {NOT_TRANSPOSED, TRANSPOSED};

 protected:
  order ordering;
  TYPE *workspace;

 public:
  //Constructors:
    //Default Constructor:

  realfftwmpi_many(): 
    ny_transpose(0), y_start_transpose(0),
      local_size(0), plan(0), iplan(0), plan_nompi(0),iplan_nompi(0),
      cdata(0), ordering(NOT_TRANSPOSED) {};

  //Copy Constructor:

  //Constructor with size and domain defined:
  realfftwmpi_many(const MPI_Comm MPI_COMM_DOMAINS, const index_type n_in[]
	      , int ntransform=1, TYPE x=0) :  ordering(NOT_TRANSPOSED)
{
  
  // An fftw plan must be established to determine the sizes in each domain:

	howmany=ntransform;
	comm = MPI_COMM_DOMAINS;

	int fftw_plan_type = FFTW_MEASURE;
#ifdef DEBUG 
	fftw_plan_type = FFTW_ESTIMATE;
#endif

	if (A_NDIM == 3) {
	  
	  plan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, 
					 n_in[0], n_in[1], n_in[2],
					 FFTW_REAL_TO_COMPLEX, fftw_plan_type);
	  
	  iplan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS,
					  n_in[0], n_in[1], n_in[2],
					  FFTW_COMPLEX_TO_REAL, fftw_plan_type);
	  
	} else if(A_NDIM == 2) {

	  plan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
					 FFTW_REAL_TO_COMPLEX, fftw_plan_type);
	  iplan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
					  FFTW_COMPLEX_TO_REAL, fftw_plan_type);
	} else {
	  /// non_mpi -- do it inplace, to be consistant with mpi setup
	  fftw_plan_type = fftw_plan_type|FFTW_IN_PLACE;
	  plan_nompi=
	    rfftwnd_create_plan(1,&n_in[0],FFTW_FORWARD,fftw_plan_type);
	  iplan_nompi=
	    rfftwnd_create_plan(1,&n_in[0],FFTW_BACKWARD,fftw_plan_type);
	}


	int fftw_local_size=0;
	if (A_NDIM>1) {
	  //This routine returns the sizes in each domain:
	  rfftwnd_mpi_local_sizes(plan, &local_nx, &x_start,
				  &ny_transpose, &y_start_transpose,
				  &fftw_local_size);
	} else {
	  // In 1D without MPI calls, we don't have a routine to do this, 
	  // so we just define values
	  local_nx = n_in[0];
	  fftw_local_size = 2*(local_nx/2+1);
	  ny_transpose=local_nx;
	  x_start=0;
	  y_start_transpose=0;
	}
	  
	// Note: the local size is padded in the first and last dimensions.
	// Yann: I don't think this last statement is actually true
	// It is padded in the last dim to make space for the complex transform
	// Yann: I agree with the last point, but...
	// and in the first dim for guard cells.
	// Yann: I don't agree with this one
	// A 3D array will have ranges of 
	//  array(0-nx_guard[0]..nx-1+nx_guard[1],0..ny-1,-0.nz+1).
	// Yann: Yes, if nx_guard = 0
	
	nsize[0]=local_nx;
	for (index_type i=1; i<A_NDIM-1; i++) nsize[i]=n_in[i];
	//Note: This matrix pads the final dimension of the array as required 
	//      by the fftw real to complex transforms so that for 3-D, 
	//      nz = 2*(n_in[2]/2+1)
	nsize[A_NDIM-1]=2*(n_in[A_NDIM-1]/2+1);
	nsize[A_NDIM] = howmany;
	if (A_NDIM <= 2) nsize[3]=0;
	if (A_NDIM <= 1) nsize[2]=0;
	
	// Generate the desired size:
	// Yann: basically, you need a working array that an handle either
	//       local_nx*ny*[2*(nz/2+1)]  -or-
	//       ny_transpose*nx*[2*(nz/2+1)]
	//       For the mpi routines, each process will possibly have an
	//       ny_transpose /= ny
	//       This is obviously not a conern in the non-mpi call/routine.
	//       The following is from the original realfftwmpi_many class, it
	//       works, so I'm keeping it...

	index_type memory_after_0=nsize[0];
	int nynz=1;
	if (A_NDIM > 1 ) nynz=nsize[1];
	if (A_NDIM == 3) nynz *= nsize[A_NDIM-1];
	memory_after_0 *= nynz;

	// Set the starting index for the array:
	int  nstart[]={0,0,0};
	
	if (A_NDIM >1 ) nstart[0]=x_start;

	if (allocatable()) {
	  if (memory_after_0 >= fftw_local_size) { 
	    //Our array supplies sufficient memory:
	    workspace=0;
	    local_size=memory_after_0*howmany;
	    data=ArrayNd_ranged<TYPE,A_NDIM+1>(nstart,nsize,x);
	  } else { 
	    // We need a bigger worksize than the array supplies:
	    local_size=fftw_local_size*howmany;
	    workspace=new TYPE[local_size];
	    data.build_ArrayNd_ranged(nstart, nsize, x, workspace);
	  } 

	
	  // transposed array
	  // Generally {0,0,0}
	  int nstart_transpose[]={y_start_transpose,0,0,0}; 
	
	  // Dimension of 3D real array
	  index_type nsize_transpose[]={ny_transpose,n_in[0],nsize[2],howmany}; 
	  // Dim of 2D real array
	  if (A_NDIM == 2) nsize_transpose[1] *= 2;
	  // Dim of 1D real array
	  if (A_NDIM == 1) nsize_transpose[0] = 2*(local_nx/2+1);
	  nsize_transpose[A_NDIM] = howmany;
	  int index_zero[4]={x_start,0,0,0};
	  TYPE* start_add=data.address(index_zero);
	  if (nsize_transpose[0]>0)
	    data_transposed.build_ArrayNd_ranged(nstart_transpose, 
						 nsize_transpose, 
						 x, start_add);

	  //Dimension of complex array
	  index_type ncsize_transpose[]={ny_transpose,n_in[0],nsize[2]/2,howmany};
	  if (A_NDIM==1) ncsize_transpose[0] = local_nx/2+1;
	  ncsize_transpose[A_NDIM] = howmany; // catches non-3D case
	  if (ncsize_transpose[0]>0) {
	    cdata.build_ArrayNd_ranged(nstart_transpose, ncsize_transpose, 
				       std::complex<TYPE>(0.,0.), 
				       (std::complex<TYPE> *)(start_add));
	  }
	}
      };

    ~realfftwmpi_many()
      {
	if (workspace !=0){
	  delete [] workspace;
	  workspace=0;
	}
      }

//  Access and setting -- Speed is everything
// if you pass only 1 index, assumes you want the 0th transform...
TYPE& operator() (int n0){ return data(n0,0);}
TYPE& operator[] (int n0){ return data(n0,0);}
// otherwise just access the value you'd expect
TYPE& operator() (int n0, int n1){
  int index[A_NDIM+1];
  index[0] = n0;
  index[1] = n1;
  if (A_NDIM==2)
    index[A_NDIM] = 0;

  if (ordering==TRANSPOSED) return data_transposed(index);
  else return data(index);
}
TYPE& operator() (int n0, int n1, int n2){
  int index[A_NDIM+1]=0;
  index[0] = n0;
  index[1] = n1;
  index[2] = n2;
  if (A_NDIM==3)
    index[A_NDIM] = 0;
  if (ordering==TRANSPOSED) return data_transposed(index);
  else return data(index);
}

//Fourier Transform Array:
int transform(TYPE* work=0, int nwork=0) {
      
      bool work_defined = false;
      if (A_NDIM>1) { // need to define work, otherwise in place
	if (work==0) {
	  work=new TYPE[local_size];
	  work_defined = true;
	} else if (nwork < local_size) {
	  // just expand work
	  nwork=local_size;
	  delete [] work;
	  work=new TYPE[local_size];
	}
      } 
	      
      /* Now, compute the forward transform: */
      int index_zero[4]={x_start,0,0,0};
      TYPE* start_add;
      if (allocatable()) {
	start_add=data.address(index_zero);
      }
      if (A_NDIM>1) 
	rfftwnd_mpi(plan, howmany, start_add, work, FFTW_TRANSPOSED_ORDER);
      else  {
	rfftwnd_real_to_complex(plan_nompi,howmany,
			    start_add,
			    howmany,1
			    ,NULL,
			    howmany,1);
      }
      
      ordering=TRANSPOSED;
    
      if (work_defined) {
	  delete [] work;
	  nwork=0;
      }
     
      return nwork; // Return the size of the work array that remains
    }

    //Inverse Fourier Transform Array:
    int invtransform(TYPE* work=0, int nwork = 0) {
      
      bool work_defined = false;
      if (A_NDIM>1) { // need to define work, otherwise in place
	if (work==0) {
	  work=new TYPE[local_size];
	  work_defined = true;
	} else if (nwork < local_size) {
	  // just expand work
	  nwork=local_size;
	  delete [] work;
	  work=new TYPE[local_size];
	}
      }
      int index_zero[4]={x_start,0,0,0};
      TYPE* start_add;
      if (allocatable()) {
	start_add=data.address(index_zero);
      }

      /* Now, compute the forward transform: */
      if (A_NDIM>1) 
	rfftwnd_mpi(iplan, howmany,  start_add, work, 
		    FFTW_TRANSPOSED_ORDER);
      else {
	rfftwnd_complex_to_real(iplan_nompi,howmany,
				(fftw_complex *)cdata.address(index_zero),
				howmany,1,
				NULL,
				howmany,1);
      }

      ordering=NOT_TRANSPOSED;

      if (work_defined) {
	  delete [] work;
	  nwork=0;
      }


      return nwork; // Return the size of the work array that remains
    
    }

    inline bool allocatable() {
      return local_nx > 0;
    }
    
};

#endif
