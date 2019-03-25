// ArrayNd_fftw2<Type, A_NDIM>:  A simple A_NDIM array class - 
//                       designed to be entirely inline.

//Note: This matrix pads the final dimension of the array as required by 
//      fftw real to complex transforms so that for 3-D,
//      nz = 2*(n_in[2]/2+1)

#ifndef ArrayNd_fftw2_H 
#define ArrayNd_fftw2_H

#include "ArrayNd.h"

#include <rfftw.h>

template <class TYPE,int A_NDIM>
class ArrayNd_fftw2 : public ArrayNd<TYPE,A_NDIM>
{
  // All values inherited from ArrayND are local!
 protected:
  using ArrayNd<TYPE,A_NDIM>::nsize; 
  using ArrayNd<TYPE,A_NDIM>::nelem;
  using ArrayNd<TYPE,A_NDIM>::data; 
  using ArrayNd<TYPE,A_NDIM>::check_dimensions;

 public:
  
  // These values are to be compatible with fftw_mpi:
  static const int x_start=0; // x start position in the global array
  static const int ny_transpose=0, y_start_transpose=0; // Local transposed values
  int local_size; //Total local size

  rfftwnd_plan plan, iplan;
  fftw_complex *cdata;

  //Constructors:
  
  //Default Constructor

  ArrayNd_fftw2() : ArrayNd<TYPE,A_NDIM>() {
    cdata=0;
    plan=0;
    iplan=0;
    local_size=0;
}

  //Copy Constructor:
  ArrayNd_fftw2(const ArrayNd_fftw2 &in) : ArrayNd<TYPE,A_NDIM>(in) {
    cdata=in.cdata;
    plan=in.plan;
    iplan=in.plan;
    local_size=in.local_size;
  }

  //Constructor with size and domain defined:
  ArrayNd_fftw2(const index_type n_in[], TYPE x=0) : ArrayNd<TYPE,A_NDIM>()
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

      /* I should be able to do the following but I cannot
      index_type n_in3[A_NDIM];
      for (index_type i=0; i<A_NDIM; i++) n_in3[i]=n_in2[i];

      ArrayNd<TYPE,A_NDIM>::ArrayNd<TYPE,A_NDIM> (*this) (n_in3, x);
      so..  I do the following */

      for (index_type i=0; i<A_NDIM; i++) nsize[i]=n_in2[i];
      check_dimensions();

      // nelem contains the number of elements between sucessive elements
      nelem[A_NDIM-1]=nsize[A_NDIM-1];
      for (int id = A_NDIM-2; id>=0; id--)  nelem[id]=nsize[id]*nelem[id+1];

      data=new TYPE[nelem[0]];
      for (index_type ip=0;ip<nelem[0];ip++) data[ip]=x;


      //For compatability with the mpi version set the following:
      local_size=nelem[0];
    }

  const ArrayNd_fftw2& operator= (ArrayNd<TYPE,A_NDIM> &in) {
    // If they're the same size copy the data except for the extra space
    // If they're different sizes, reinitialize:
    int same=true;
    for (index_type i=0; i<A_NDIM-1; i++) if (nsize[i] != in.size(i)) same=false;
    if (nsize[A_NDIM-1] != 2*(in.size(A_NDIM-1)/2+1)) same=false;
    
    if (same!=true) {
      index_type old_nelem=nelem[0];
      for (index_type i=0; i<A_NDIM; i++) nsize[i]=in.size(i);
      nsize[A_NDIM-1]=2*(nsize[A_NDIM-1]/2+1);
      // nelem contains the number of elements between sucessive elements
      nelem[A_NDIM-1]=nsize[A_NDIM-1];
      for (int id = A_NDIM-2; id>=0; id--)  nelem[id]=nsize[id]*nelem[id+1];

      delete [] data;
      if  (nelem[0] > 0) data=new TYPE[nelem[0]];
      else data=0;
    }

    // Copy old to new:
    int id2=0;
    for (index_type id=0;id<nelem[0];id++) 
      if (id%nsize[A_NDIM-1] < in.size(A_NDIM-1) ) {
	data[id]= in(id2);
	id2++;
      }
    else data[id]=0;
    
    return *this;
  };


  //Fourier Transform Array:
  void transform() {
    
    /* Now, compute the forward transform: */

    rfftwnd_one_real_to_complex(plan, data, cdata);

    /* the data is now complex, so typecast a pointer: */
     cdata = (fftw_complex*) data;   
     
  }

  //Inverse Fourier Transform Array:
  void invtransform() {
    
    /* Now, compute the forward transform: */
    rfftwnd_one_complex_to_real(iplan, cdata, data);


  }

};




#endif
