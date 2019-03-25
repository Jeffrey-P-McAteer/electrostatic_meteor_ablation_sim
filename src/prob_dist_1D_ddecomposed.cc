

#include "eppic-types.h"
#include "eppic-mpi.h"

/** Calculated the max value and max local average of a 1D function
    across domains. 
**/

void prob_dist_1D_ddecomposed(FTYPEAVec &local_dist,FTYPE &global_max,
			      FTYPE &global_avg_max, FTYPE &right_area, FTYPE &left_area) {
  int nx = local_dist.length();
  FTYPE local_max = local_dist.max();
  FTYPE local_avg = local_dist.sum();
  if (nsubdomains>1) {
    
    FTYPEAVec domain_avgs = FTYPEAVec(nsubdomains);

    /*
    MPI_Allreduce(&local_avg,&global_avg_max,1,MPI_FTYPE,
		  MPI_MAX,subdomain.neighbor_comm);
    */
    MPI_Allgather(&local_avg,1,MPI_FTYPE,
		  domain_avgs.address(0),1,MPI_FTYPE,
		  subdomain.neighbor_comm);

    FTYPE total_area = 0;

    total_area=domain_avgs.sum();
    left_area = 0;
    for (int idomain = 0; idomain < subdomain.id_number; idomain++) 
      left_area+=domain_avgs(idomain);
    left_area/=total_area;
    right_area = left_area+local_avg/total_area;
    global_avg_max = domain_avgs.max()/nx;
    /* //old way, when x=-1 was an acceptable position
    if (subdomain.id_number==0) {
      if (global_avg_max==domain_avgs(0))
	global_avg_max/=nx;
      else
	global_avg_max/=(nx-1);
    }
    else {
      if (global_avg_max==domain_avgs(0))
	global_avg_max/=(nx+1);
      else
	global_avg_max/=nx;
    }
    */
    
    MPI_Allreduce(&local_max,&global_max,1,MPI_FTYPE,
		  MPI_MAX,subdomain.neighbor_comm);
    MPI_Request send_request,recv_request;
    MPI_Status  status;

    // get boundary cell for correctly laying out particles    
    if (subdomain.id_number>0) 
      MPI_Isend(&local_dist(0),1, MPI_FTYPE,
		subdomain.id_number-1, 2,
		subdomain.neighbor_comm, &send_request);
    if (subdomain.id_number<nsubdomains-1)
      MPI_Irecv(&local_dist(nx-1), 1, MPI_FTYPE,
		subdomain.id_number+1, 2,
		subdomain.neighbor_comm, &recv_request);
    if (subdomain.id_number>0)
      MPI_Wait(&send_request, &status);
    if (subdomain.id_number<nsubdomains-1)
      MPI_Wait(&recv_request, &status);

  } else {
    global_max = local_max;
    global_avg_max = local_avg;
  }



}
