// ArrayNd_fftwmpi<Type, A_NDIM>:  A simple A_NDIM array class - 
//                       designed to be entirely inline.

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)

#ifndef ArrayNd_fftwmpi_H 
#define ArrayNd_fftwmpi_H

#include "ArrayNd.h"

#include <rfftw_mpi.h>

template <class TYPE,int A_NDIM>
class ArrayNd_fftwmpi : public ArrayNd<TYPE,A_NDIM>
{
  // All values inherited from ArrayND are local!

 protected:
  // Because of the templating, c++ requires the following:
  using ArrayNd<TYPE,A_NDIM>::nsize; 
  using ArrayNd<TYPE,A_NDIM>::nelem;
  using ArrayNd<TYPE,A_NDIM>::data; 
  using ArrayNd<TYPE,A_NDIM>::check_dimensions;

 public:

  int x_start; // x start position in the global array
  int ny_transpose, y_start_transpose; // Local transposed values
  int local_size; //Total local size

  rfftwnd_mpi_plan plan, iplan;
  fftw_complex *cdata;

  //Constructors:
  
  //Default Constructor

  ArrayNd_fftwmpi() : ArrayNd<TYPE,A_NDIM>(), x_start(0),  
    ny_transpose(0), y_start_transpose(0),
    local_size(0), plan(0), iplan(0), cdata(0) {}

  //Copy Constructor:
  ArrayNd_fftwmpi(const ArrayNd_fftwmpi &in) : ArrayNd<TYPE,A_NDIM>(in) {
    x_start=0; 
    ny_transpose=0;
    y_start_transpose=0;
    cdata=in.cdata;
    plan=in.plan;
    iplan=in.plan;
    local_size=in.local_size;
  }

  //Constructor with size and domain defined:
  ArrayNd_fftwmpi(const MPI_Comm MPI_COMM_DOMAINS, const index_type n_in[], TYPE x=0) 
    : ArrayNd<TYPE,A_NDIM>()
    {

      // An fftw plan must be established to determine the sizes in each domain:
      int j;

      cout << mpi_rank << " sizes: "<< n_in[0] << " " << n_in[1] << " " << n_in[2] <<"\n" ;

    int fftw_plan_type = FFTW_MEASURE;
#ifdef DEBUG 
   fftw_plan_type = FFTW_ESTIMATE;
#endif

      if (A_NDIM == 3) {

	  plan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1], n_in[2],
				       FFTW_REAL_TO_COMPLEX,  fftw_plan_type);

	  iplan = rfftw3d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1], n_in[2],
					  FFTW_COMPLEX_TO_REAL,  fftw_plan_type);

      } else if(A_NDIM == 2) {

	plan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
				       FFTW_REAL_TO_COMPLEX, fftw_plan_type);

	iplan = rfftw2d_mpi_create_plan(MPI_COMM_DOMAINS, n_in[0], n_in[1],
					FFTW_COMPLEX_TO_REAL, fftw_plan_type);
      } else {
	ArrayNd<TYPE,A_NDIM>::check_error(0==0, "ArrayNd_fftwmpi only 2D or 3D implemented");
      }

      cout << mpi_rank << " plan done" <<"\n" ;

      //This routine returns the sizes in each domain:
      int local_nx;

      rfftwnd_mpi_local_sizes(plan, &local_nx, &x_start,
			      &ny_transpose, &y_start_transpose,
			      &local_size);

      // Note: the local size is padded in the last dimesion 
      // to make space for the complex transform
 
      //      cout << mpi_rank << " local_nx=" <<  local_nx << " x_start=" <<  x_start 
      //   << " ny_transpose=" << ny_transpose << " y_start_transpose=" 
      //   <<  y_start_transpose << " local_size=" << local_size << "\n";
      
      nsize[0]=local_nx;
      for (index_type i=1; i<A_NDIM; i++) nsize[i]=n_in[i];

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)
      nsize[A_NDIM-1]=2*(nsize[A_NDIM-1]/2+1);
      check_dimensions();

      // nelem contains the number of elements between sucessive elements
      nelem[A_NDIM-1]=nsize[A_NDIM-1];
      for (int id = A_NDIM-2; id>=0; id--)  nelem[id]=nsize[id]*nelem[id+1];

      delete [] data;
      data=new TYPE[local_size];

      for (index_type ip=0;ip<nelem[0];ip++) data[ip]=x;

    }

  //  Assignment - I should not have to do this here.
  const ArrayNd_fftwmpi& operator= (TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i]=x;
    return *this;
  };


  //Fourier Transform Array:
  void transform() {
    
    TYPE work[local_size];
    
    /* Now, compute the forward transform: */
    rfftwnd_mpi(plan, 1, data, work, FFTW_NORMAL_ORDER);

    /* the data is now complex, so typecast a pointer: */
     cdata = (fftw_complex*) data;   
     
  }

  //Inverse Fourier Transform Array:
  void invtransform() {
    
    TYPE work[local_size];
    
    /* Now, compute the forward transform: */
    rfftwnd_mpi(iplan, 1, data, work, FFTW_NORMAL_ORDER);

  }

};




#endif
