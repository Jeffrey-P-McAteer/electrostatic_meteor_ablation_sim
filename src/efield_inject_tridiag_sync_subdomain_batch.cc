

#include "efield_inject_tridiag.h" /* for field, types */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */

/// In 2D case, ny is subdomain-decomposed, this routine passes the local
///  portion of the working arrays to the other processors.
    
///  This is used if there are multiple processors per subdomain. 

///  \param rfftw_many working_array -- solutions on the grid
///  \param rfftw_many boundary_working_array -- solutions at the boundary
///  \param field      Efield -- used for frequnetly used values
///
///  \todo rename variables. The names were chosen when it was thought
///        that both 2D and 3D varients needed this code

void  efield_inject_tridiag_sync_subdomain_batch(rfftw_many &working_array,
					   rfftw_many &boundary_working_array,
					   field &Efield)
{

  if (subdomain.np==1) return;
#if NDIM==3
  return;
#else
  MPI_Barrier(MPI_COMM_WORLD);
  int nx = Efield.nx;
  int ny = Efield.ny;


  int nNy=ny/2+1;

  int dnslow = nNy/subdomain.np;
  // old way: static_cast<int>(ceil((nNy)/subdomain.np));
  int slow_start = dnslow*subdomain.rank;
  int slow_end = slow_start+dnslow;
  if (slow_end+dnslow > nNy) slow_end = nNy; 
  int sendcount=2*nx*(slow_end-slow_start);

  ArrayNd<FTYPE,2> send_buf = ArrayNd<FTYPE,2>(2*(slow_end-slow_start),nx);

  for(int iy=slow_start*2;iy<slow_end*2;iy++)
    for (int ix=0;ix<nx;ix++)
      send_buf(iy-slow_start*2,ix) = working_array.data(iy,ix);
  
  int recvcounts[subdomain.np];
  int displs[subdomain.np];

  for (int ip=0;ip<subdomain.np;ip++) {
    int ip_start=dnslow*ip;
    int ip_end = ip_start+dnslow;
    if (ip_end +dnslow > nNy) ip_end = nNy; 
    recvcounts[ip]=2*nx*(ip_end-ip_start);
    displs[ip]=2*ip_start*(nx);
  }

  int zero_addr[2]={0,0};

  MPI_Allgatherv(send_buf.address(zero_addr),
		 sendcount,MPI_FTYPE,
		 working_array.data.address(zero_addr),
		 recvcounts,displs,MPI_FTYPE,
		 subdomain.internal_comm);


  for (int ip=0;ip<subdomain.np;ip++) {
    recvcounts[ip]=(recvcounts[ip]*3)/nx;
    displs[ip]=(displs[ip]*3)/nx;
  }

  send_buf = ArrayNd<FTYPE,2>(2*(slow_end-slow_start),3);

  for(int iy=slow_start*2;iy<slow_end*2;iy++)
    for (int ix=0;ix<3;ix++)
      send_buf(iy-slow_start*2,ix) = 
	boundary_working_array.data(iy,ix);
  
  sendcount=(sendcount*3)/nx;
  MPI_Allgatherv(send_buf.address(zero_addr),
		 sendcount,MPI_FTYPE,
		 boundary_working_array.data.address(zero_addr),
		 recvcounts,displs,MPI_FTYPE,
		 subdomain.internal_comm);


#endif
}
