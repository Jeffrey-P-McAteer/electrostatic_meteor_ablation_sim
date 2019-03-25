
#include "eppic-types.h"

void transpose2d(ArrayNd<FTYPE,2> &in_array) {
  
//   transposed_array.build_ArrayNd(in_array.size(1),in_array.size(0),
// 				 in_array.address(0,0);
  
  if (in_array.size(0) == in_array.size(1)) {
    // square matrix is a special, fast case
    FTYPE value_holder=0;
    for (int ix=0;ix<in_array.size(0);ix++) {
      for (int iy=ix+1;iy<in_array.size(1);iy++) {
	value_holder = in_array(ix,iy);
	in_array(ix,iy) = in_array(iy,ix);
	in_array(iy,ix) = value_holder;
      }
    }
  } else {
    // rectangular, much more expensive
    // there are faster routines, but I'm not doing them now.
    ArrayNd<FTYPE,2> transposed_array;
    int dimensions[2]={in_array.size(1),in_array.size(0)};
    transposed_array.build_ArrayNd(dimensions);
    for (int ix=0;ix<in_array.size(0);ix++) {
      for (int iy=0;iy<in_array.size(1);iy++) {
	transposed_array(iy,ix) = in_array(ix,iy);
      }
    }
    in_array = transposed_array;
  }
  
}
