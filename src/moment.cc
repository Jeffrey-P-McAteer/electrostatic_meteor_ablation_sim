// Routine to calculate the n first moments of x where n is determined by 
// the size of moments and x is a vector.
// On return:
//    moment(0) gives the mean value of x called <x>
//    moment(1) gives <(x-<x>)^2)
//    moment(2) gives <(x-<x>)^3)
//    ...
//    moment(n) gives <(x-<x>)^n)
  //The binomial theorem is applied to efficiently calculate the moments
  //using a single pass algorithm

#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

int moment(PTYPEAVec &x, PTYPEAVec &moment, PTYPEAVec &flag_array, int np=0)
{
  
  // For each member of the distribution:
  if (np==0) np=x.size();

  int nmom=moment.size();
  ArrayNd<double,1> sum(nmom+1);
  sum=0.0;

  if (nmom<1) {
      
      
      return -2; // Error
  }
  //    ArrayNd<double,1> sum(nmom+1);
  bool use_flag=false;
  if (flag_array.size() >= np) use_flag=true;
  
  // Calculate the sums of the products:
  for (int i = 0; i < np ; ++i) {
      if (use_flag && flag_array(i)>0) {
	  PTYPE prod=1;
	  for (int imom=0;imom<nmom;imom++) {
	  prod *= x(i);
	  sum(imom) += prod;
	  }
	  // Store the number of particles (not flagged)
	  sum(nmom) += 1;
      }
  }

  // Add across processors if necessary
#ifdef USE_MPI
    {
      
      ArrayNd<double,1> sum_all(nmom+1);  
      int mpi_err=MPI_Reduce((void*) &(sum[0]),(void*) &(sum_all[0]), 
			     nmom+1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if ( mpi_err != MPI_SUCCESS) 
	mpi_error(mpi_rank, mpi_err,"Reduce call in output");
      sum = sum_all;
    }
#endif    
  if (sum(nmom) >0) sum /= sum(nmom); else sum=0;

  //Use the binomial theorem to turn these into moments
  moment(0)=sum(0);
  for (int imom=2;imom<nmom+1;imom++) {
    moment(imom-1) = pow(-sum(0),imom);
    PTYPE coeff=1;
    PTYPE fact=1;
    for (int iseq=1;iseq<=imom;iseq++) {
      coeff *= imom-(iseq-1);
      fact *= iseq;
      moment(imom-1) += (coeff/fact)*sum(iseq-1)*pow(-sum(0),imom-iseq);
    }
  }

  
  
  return 0; //Successful completion!
}
  

