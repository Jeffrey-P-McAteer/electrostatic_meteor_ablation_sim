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

#include "eppic_fftw_mpi.h"
#include<complex>

template <class TYPE,int A_NDIM>
class realfftwmpi 
{

 public:

  ArrayNd_ranged<TYPE,A_NDIM> data;
  rfftwnd_mpi_plan plan, iplan;
  ArrayNd_ranged<std::complex<TYPE>,A_NDIM> cdata;
  ArrayNd_ranged<TYPE,A_NDIM> data_transposed; //Container for transposed data

  int x_start; // x start position in the global array
  int ny_transpose, y_start_transpose; // Local transposed values
  int local_size; //Total local size
  int local_nx; // nx on the local processor

  enum order {NOT_TRANSPOSED, TRANSPOSED};
  order ordering;

  FTYPE *workspace;

 protected:

 public:

  //Constructors:
    //Default Constructor:

  realfftwmpi(): 
    ny_transpose(0), y_start_transpose(0),
    local_size(0), plan(0), iplan(0), cdata(0), ordering(NOT_TRANSPOSED) {};

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
    realfftwmpi(const MPI_Comm MPI_COMM_DOMAINS, const index_type n_in[], TYPE x=0) :
      ordering(NOT_TRANSPOSED)
      {

	// An fftw plan must be established to determine the sizes in each domain:
	
#ifdef DEBUG 
	const int fftw_plan_type = FFTW_ESTIMATE;
#else
	const int fftw_plan_type = FFTW_MEASURE;
#endif

	if (A_NDIM == 3) {
	  
	  plan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1], n_in[2],
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
	  cerr << "realfftwmpi only 2D or 3D implemented"<<"\n" ;
	  exit(-1);
	}

	make_workspace(n_in, x);
      };

    ~realfftwmpi()
      {
	del_workspace();
      }

    //  Access and setting -- Speed is everything
    TYPE& operator() (int n0){ return data(n0);}
    TYPE& operator[] (int n0){ return data[n0];}
    TYPE& operator() (int n0, int n1){
      if (ordering==TRANSPOSED) return data_transposed(n0,n1);
      else return data(n0,n1);
    }
    TYPE& operator() (int n0, int n1, int n2){
      if (ordering==TRANSPOSED) return data_transposed(n0,n1,n2);
      else return data(n0,n1,n2);
    }

    // To enable one to create and delete the (large) storage of this array, 
    // but keep everything else, we add del_workspace and make_workspace functions.
    void del_workspace() {
      if (workspace !=0) {
	delete [] workspace;
	workspace=0;
      }
    }
    
    //Make and organize the work space
    void make_workspace(const index_type n_in[], TYPE x) {

      //This routine returns the sizes in each domain:
      
      int fftw_local_size=0;
      rfftwnd_mpi_local_sizes(plan, &local_nx, &x_start,
			      &ny_transpose, &y_start_transpose,
			      &fftw_local_size);
      
      // Note: the local size is padded in the first and last dimensions.
      // It is padded in the last dim to make space for the complex transform
      // and in the first dim for guard cells.
      // A 3D array will have ranges of array(0-nx_guard[0]..nx-1+nx_guard[1],0..ny-1,-0.nz+1).
      
      index_type nsize[]={0,0,0};
      nsize[0]=local_nx;
      for (index_type i=1; i<A_NDIM-1; i++) nsize[i]=n_in[i];
	//Note: This matrix pads the final dimension of the array as required 
	//      by the fftw real to complex transforms so that for 3-D, 
	//      nz = 2*(n_in[2]/2+1)
      nsize[A_NDIM-1]=2*(n_in[A_NDIM-1]/2+1);

	// Generate the desired size:
	index_type memory_after_0=local_nx;
	int nynz=nsize[1];
	if (A_NDIM == 3) nynz *= nsize[A_NDIM-1];
	memory_after_0 *= nynz;

	// Set the starting index for the array:
	int  nstart[]={0,0,0};

	if (memory_after_0 >= fftw_local_size) { //Our array supplies sufficient memory:
	  workspace=0;
	  local_size=memory_after_0;
	  data=ArrayNd_ranged<TYPE,A_NDIM>(nstart,nsize,x);
	} else { // We need a bigger worksize than the array supplies:
	  local_size=fftw_local_size;
	  // delete [] workspace;
	  workspace=new TYPE[local_size];
	  data.build_ArrayNd_ranged(nstart, nsize, x, workspace);
	} 

	/*transposed array*/
	int nstart_transpose[]={nstart[1],nstart[0],nstart[2]}; // Generally {0,0,0}
	// Dimension of 3D real array
	index_type nsize_transpose[]={ny_transpose,n_in[0],nsize[2]}; 
	if (A_NDIM == 2) nsize_transpose[1] *= 2; // Dim of 2D real array
	if (A_NDIM == 1) nsize_transpose[0] =(nstart[0])*2; // Dim of 1D real array
	int index_zero[3]={0,0,0};
	TYPE* start_add=data.address(index_zero);
	if (nsize_transpose[0]>0)
	  data_transposed.build_ArrayNd_ranged(nstart_transpose, nsize_transpose, x, start_add);

	//Dimension of complex array
	index_type ncsize_transpose[]={ny_transpose,n_in[0],nsize[2]/2};
	if (ncsize_transpose[0]>0)
	  cdata.build_ArrayNd_ranged(nstart_transpose, ncsize_transpose, 
			      std::complex<TYPE>(0.,0.), 
			      (std::complex<TYPE> *)(start_add));
    }      


    //Fourier Transform Array:
    int transform(TYPE* work=0, int nwork=0) {
      
      bool work_defined = false;
      if (work==0) {
	work=new TYPE[local_size];
	work_defined = true;
      } else if (nwork < local_size) {
	  // just expand work
	  nwork=local_size;
	  delete [] work;
	  work=new TYPE[local_size];
      }
	      
      /* Now, compute the forward transform: */
      int index_zero[3]={0,0,0};
      TYPE* start_add=data.address(index_zero);
      rfftwnd_mpi(plan, 1, start_add, work, FFTW_TRANSPOSED_ORDER);
      
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
      if (work==0) {
	work=new TYPE[local_size];
	work_defined = true;
      } else if (nwork < local_size) {
	  // just expand work
	  nwork=local_size;
	  delete [] work;
	  work=new TYPE[local_size];
      }

      int index_zero[3]={0,0,0};
      /* Now, compute the forward transform: */
      rfftwnd_mpi(iplan, 1,  data.address(index_zero), work, FFTW_TRANSPOSED_ORDER);
    
      ordering=NOT_TRANSPOSED;

      if (work_defined) {
	  delete [] work;
	  nwork=0;
      }


      return nwork; // Return the size of the work array that remains
    
    }

};

#endif
