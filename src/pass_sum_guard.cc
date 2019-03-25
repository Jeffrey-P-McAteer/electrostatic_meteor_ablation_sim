  // Passe the guard cell data in rho(nx,*,*) to the east and sums it to
  // rho(0,*,*)

#include "eppic.h"
#include "eppic-mpi.h"

void pass_sum_guard(FArrayND &rho, int nx_local)
{
  
  MPI_Request request1;
  MPI_Status  status1;
  //  int size_send=rho.length()/rho.size(0);
  long size_send=rho.length()/rho.size(0);
  
  // Pass the last row:
  const int index0[]={nx_local,0,0};
  if (proc_east < mpi_np) {
    MPI_Isend(rho.address(index0), size_send, MPI_FTYPE, proc_east, 
	      1, MPI_COMM_WORLD, &request1);
  }

  if (proc_west >= 0) {
    FTYPEAVec receive(size_send);
    MPI_Recv(receive.address(), size_send, MPI_FTYPE, proc_west, 
	     1, MPI_COMM_WORLD, &status1);

    // Add receive to den(0,*,*)
    // if statement added by Glenn 3/1/2018 to fix OPEN boundary bug
    if (!(boundary_type[0] != PERIODIC && subdomain.id_number == 0)){
      int ic=0,ny2=1,nz2=1;
      if (NDIM > 1) ny2=rho.size(1);
      if (NDIM > 2) nz2=rho.size(2);
      for (int iy=0; iy<ny2; iy++)
        for (int iz=0; iz<nz2; iz++)
          rho(INDICIES(0,iy,iz))+=receive(ic++);
    }
  }
  if (proc_east < mpi_np) {
    MPI_Wait(&request1, &status1);
  }
}	  

void pass_sum_guard(ArrayNd<float, NDIM> &rho, int nx_local)
{
  
  MPI_Request request1;
  MPI_Status  status1;
  //  int size_send=rho.length()/rho.size(0);
  long size_send=rho.length()/rho.size(0);  

  // Pass the last row:
  const int index0[]={nx_local,0,0};
  if (proc_east < mpi_np) {
    MPI_Isend(rho.address(index0), size_send, MPI_FLOAT, proc_east, 
	      1, MPI_COMM_WORLD, &request1);
  }
  if (proc_west >= 0) {
    ArrayNd<float, 1> receive(size_send);
    MPI_Recv(receive.address(), size_send, MPI_FLOAT, proc_west, 
	     1, MPI_COMM_WORLD, &status1);
   

  /*  if (boundary_type[0] == INJECT && proc_west <0) {
      This is the left injection row and the density needs supplementing.
       We will just add n0/2.  This makes the injection between left and right assymmetric.
       NOTE: A guard cell needs to be added to rho and the density(-1,*,*)
       would then contain the correct amount.  This requires substantial revisions 
       throughout the code
       ...Not done */

    // Add receive to den(0,*,*)
    // if statement added by Glenn 3/1/2018 to fix OPEN boundary bug
    // Don't add to den(0,*,*) if 0th subdomain and non-periodic particle x boundary
    if (!(boundary_type[0] != PERIODIC && subdomain.id_number == 0)){
      int ic=0,ny2=1,nz2=1;
      if (NDIM > 1) ny2=rho.size(1);
      if (NDIM > 2) nz2=rho.size(2);
      for (int iy=0; iy<ny2; iy++)
        for (int iz=0; iz<nz2; iz++)
          rho(INDICIES(0,iy,iz))+=receive(ic++);
    }
  }
  if (proc_east < mpi_np) {
    MPI_Wait(&request1, &status1);
  }

}	  
