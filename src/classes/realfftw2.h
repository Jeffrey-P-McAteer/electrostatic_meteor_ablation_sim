// realfftw2<Type, A_NDIM>:  A simple A_NDIM array class - 
//                       designed to be entirely inline.

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)

#ifndef realfftw2_H 
#define realfftw2_H

// We use the following class type:
#include "ArrayNd.h"

#ifdef EPPIC_FFTW_USE_D_PREFIX 
#include <drfftw.h>
#else
#include <rfftw.h>
#endif

template <class TYPE,int A_NDIM>
class realfftw2 
{
 public:

  // All value are public

  ArrayNd<TYPE,A_NDIM> data;
  rfftwnd_plan plan, iplan;
  // fftw_complex *cdata;
  ArrayNd<std::complex<TYPE>,A_NDIM> cdata;

  
  // These values are to be compatible with fftw_mpi:
  static const int x_start=0; // x start position in the global array
  static const int ny_transpose=0, y_start_transpose=0; // Local transposed values
  int local_size; //Total local size

  //Constructors:
  
  //Default Constructor

  realfftw2()  {
    cdata=0;
    plan=0;
    iplan=0;
    local_size=0;
  }
  
  //Copy Constructor:
  realfftw2(const realfftw2 &in) {
    data=in.data;
    //    cdata = (fftw_complex*) data.address();   
    cdata=build_ArrayNd(data.size(), std::complex<TYPE>(0.,0.), 
			(std::complex<TYPE> *)(data.address()));
    plan=in.plan;
    iplan=in.plan;
    local_size=in.local_size;
  }

  //Constructor with size and domain defined:
  realfftw2(const index_type n_in[], TYPE x=0)
    {
      // An fftw plan must be established to determine the sizes in each domain:
      int j;

      // Type change for n_in
      int n_in2[A_NDIM];
      for (index_type i=0; i<A_NDIM; i++) n_in2[i]=n_in[i];

    int fftw_plan_type = FFTW_MEASURE | FFTW_IN_PLACE;
#ifdef DEBUG 
   fftw_plan_type = FFTW_ESTIMATE | FFTW_IN_PLACE;
#endif
      plan = rfftwnd_create_plan(A_NDIM, n_in2, FFTW_FORWARD, fftw_plan_type);
      iplan = rfftwnd_create_plan(A_NDIM, n_in2, FFTW_BACKWARD, fftw_plan_type);

      // Note: the local size is padded in the last dimension 
      // to make space for the complex transform
      n_in2[A_NDIM-1]=2*(n_in[A_NDIM-1]/2+1);

      data=ArrayNd<TYPE,A_NDIM>(n_in2,x);

      //For compatability with the mpi version set the following:
      local_size=data.size();
    }

  //  Access and setting -- Speed is everything
  TYPE& operator() (int n0){ return data(n0);}
  TYPE& operator[] (int n0){ return data[n0];}
  TYPE& operator() (int n0, int n1){return data(n0,n1);}
  TYPE& operator() (int n0, int n1, int n2){return data(n0,n1,n2);}
  TYPE& operator() (int n0, int n1, int n2, int n3){return data(n0,n1,n2,n3);}
  TYPE& operator() (int n0, int n1, int n2, int n3, int n4){return data(n0,n1,n2,n3,n4);}
  TYPE& operator() (int n0, int n1, int n2, int n3, int n4, int n5){return data(n0,n1,n2,n3,n4,n5);}

  //Fourier Transform Array:
  void transform() {
    
    /* Now, compute the forward transform: */

    rfftwnd_one_real_to_complex(plan, data.address(), cdata.address());

    /* the data is now complex, so typecast a pointer: */
    //     cdata = (fftw_complex*) data.address();
     
  }

  //Inverse Fourier Transform Array:
  void invtransform() {
    
    /* Now, compute the forward transform: */
    rfftwnd_one_complex_to_real(iplan, cdata, data.address());


  }

};




#endif
