#include "eppic.h"
#ifdef USE_DOMAINS
#include "eppic-mpi.h"


void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[])
{
  //This routine passes guard cell(s) from an array on one processor 
  //to another on another processor.  It only passes the last dimension(s).
  // nx_guard[0] indicates how many guard cells to be passed on the LHS
  // nx_guard[1] indicates how many guard cells to be passed on the RHS

  MPI_Request request1, request2;
  MPI_Status  status1,  status2;

  int size_send=in.length()/in.size(0);
  int this_nx = in.xsize() - nx_guard[0] - nx_guard[1];
  // Pass RHS to the East:
  // added if statement by Glenn 3/1/2018 to fix OPEN boundary bug
  // Send to east processor if not periodic particle x boundary
  // and 
  //if (proc_east < mpi_np && 
  //    !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //      && subdomain.id_number == nsubdomains-1)) {
  // Pass in(nx-1,*,*) to east prox
  if (proc_east < mpi_np){
    const int index0[]={INDICIES(this_nx-nx_guard[0],0,0)};
    MPI_Isend(in.address(index0), size_send*nx_guard[0], MPI_FTYPE, 
              proc_east, 2, MPI_COMM_WORLD, &request2);
  }
  
  // Recieve this data from the West and place in the LHS:
  //if (proc_west >= 0 &&
  //    !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //      && subdomain.id_number == 0)) {
  if (proc_west >=0){
    const int index1[]={INDICIES(-((int)nx_guard[0]),0,0)};
    MPI_Recv(in.address(index1), size_send*nx_guard[0], MPI_FTYPE, 
             proc_west, 2, MPI_COMM_WORLD, &status2);
  }

  // Pass LHS to the West:
  // added if statement by Glenn 3/1/2018 to fix OPEN boundary bug
  //if (proc_west >= 0 &&
  //   !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //     && subdomain.id_number == 0)) {
  if (proc_west >= 0){
    const int index2[]={INDICIES(0,0,0)};
    MPI_Isend(in.address(index2), size_send*nx_guard[1], MPI_FTYPE, proc_west, 
              1, MPI_COMM_WORLD, &request1);
  }
  
  // Recieve this data from the East and place in the RHS:
  //  if (proc_east < mpi_np &&
  //   !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //     && subdomain.id_number == nsubdomains-1)) {
  if (proc_east < mpi_np){
    const int index3[]={INDICIES(this_nx,0,0)};
    MPI_Recv(in.address(index3), size_send*nx_guard[1], MPI_FTYPE, proc_east, 
             1, MPI_COMM_WORLD, &status1);
  }
  
  // make sure non-blocking sends have completed.
  //if (proc_west >= 0) MPI_Wait(&request2, &status2);
  //if (proc_east < mpi_np) MPI_Wait(&request1, &status1);
  //  if (proc_east < mpi_np && 
  //   !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //     && subdomain.id_number == nsubdomains-1)) {
  if (proc_east < mpi_np){
    MPI_Wait(&request2, &status2);
  }
  
  //if (proc_west >= 0 &&
  //   !(boundary_type[0] != PERIODIC && field_boundary_type[0][0] != periodic
  //     && subdomain.id_number == 0)) {
  if (proc_west >= 0){
    MPI_Wait(&request1, &status1);
  }

}
  
#endif
