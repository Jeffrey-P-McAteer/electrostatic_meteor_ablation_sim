/**
   ----------------------------------------------------

   pass particles across domains in the x direction
   use non-blocking sends
   
   This routine also initiates particle injections (REMOVED BY GLENN 4/12/2018)
   Glenn: Injection is now done in another funciton
   
   Meers Oppenheim & Doug Sondak
   
   \bug There is a memory leak at line 183, build_ArrayNd
   It likely has to do with ext_workspace not being initialized
   correctly, to the value false.

   ---------------------------------------------------- 
**/

#include <cmath>
#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-velsample.h"
#include <fstream>   // file I/O
#include <cstdio> 
#include <csignal>

void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                         bool add_absent=false);
  
void pass_particles(int id, particle_dist  *adv, int nabsent_start, int vel_dim)
{
  
  int n_per_particle=ndim+vel_dim;
  int max_absent=(int) (adv[id].x.size()*(part_pad[id]-1)/part_pad[id]);
  PTYPE fnx = (PTYPE) nx;
  
  PTYPE smallest_float_below(PTYPE x);
  static PTYPE fnx_minus_epsilon=smallest_float_below(fnx);
  static int multi_cross=0;
  
  // Max_absent is sometimes large because of creation so use the lesser of:
  int max_pass=std::min(max_absent*n_per_particle,
			(adv[id].nabsent-nabsent_start)*n_per_particle+1);
  
  PTYPEAVec send_east(max_pass);
  PTYPEAVec send_west(max_pass);
  MPI_Request request1, request2;
  MPI_Status  status1,  status2;
  
  // Loop through absent array
  int nsend_east=0;
  int nsend_west=0;
  
  int nabsent_end = adv[id].nabsent;

  // Pack the array for passing east:
  //for (int i=nabsent_start;i<nabsent_end;i++) {
  for (int j=nabsent_start;j<nabsent_end;j++) {
    //  printf("nabsent_end: %d, nabsent_start: %d", nabsent_end, nabsent_start);
    int i = nabsent_end - j + nabsent_start-1;
    int isend = adv[id].absent(i);
    
    // See if we can decrement np (if the absent particle index is the last element) Added by Glenn 2018/01/11
    if (isend == adv[id].np-1){
      adv[id].np--;
      adv[id].absent(i) = -1;
      adv[id].nabsent--;
    }

    if (adv[id].x(isend) >=  fnx) {
      // Test for the array being too small:
      if (nsend_east+n_per_particle > send_east.size()) {
	cout << "Error: insufficient space in the send_east array in distribution " 
	     << id << " rank " << mpi_rank << endl 
	     << "Need to increase part_pad" << endl;
     	terminate(-1,"");
      }
      
      // We will adjust for the boundary condition here: (periodic only now)
      PTYPE x=adv[id].x(isend)-fnx; // This is done for rounding error reasons.
      
      // Test to see if the particle has crossed more than 1 domain.  
      if (x >= fnx || x < (PTYPE)0. ) {
	// If this happens a lot or the particle is really moving, terminate
	if (multi_cross++ >= adv[id].np*0.01){
	  terminate(-1,"Too many particles cross multiple domains");
	}
	if  (x >= 2*fnx || x < -1*fnx ) {
	  terminate(-1,"Particle cross more than 2 domains");
	}
	// This is just a rare event, so I will issue a warning and
	// adjust the particles position and velocity, 
	// Reduce the velocity so the particle just reaches the far side 
	// instead of overshooting:
	PTYPE x_orig=x+fnx-adv[id].vx(isend)*subcycle[id];
	PTYPE x_new=fnx_minus_epsilon;
	PTYPE vx_new=( (x_new+fnx) -x_orig);
	cout << "Warning: "
	     << "Particle jumps more than 1 domain "
	     << "- adjusting to prevent. Dist: " 
	     << id << " rank " << mpi_rank << endl 
	     << "Old position and velocity " 
	     << x <<"  " << adv[id].vx(isend) <<endl
	     << "New position and velocity " 
	     << x_new <<"  " << vx_new <<endl;
	adv[id].vx(isend)=vx_new;
	x=x_new;
      }
      
      // only populate send_east if there is a processor to the east willing to accept
      if (!(boundary_type[0] != PERIODIC && subdomain.id_number == nsubdomains-1)){
        send_east(nsend_east++) = x;
        if(ndim > 1) send_east(nsend_east++) = adv[id].y(isend);
        if(ndim > 2) send_east(nsend_east++) = adv[id].z(isend);
        send_east(nsend_east++) = adv[id].vx(isend);
        if(vel_dim > 1) send_east(nsend_east++) = adv[id].vy(isend);
        if(vel_dim > 2) send_east(nsend_east++) = adv[id].vz(isend);
      }

      // Set the existing particle to an indicator value, 
      // outside all possible values.
      adv[id].x(isend)=fabsent2;
      if (NDIM >= 2) adv[id].y(isend)=fabsent2;
      if (NDIM >= 3) adv[id].z(isend)=fabsent2;
      // Insure that it does not move (also removes the energy)
      adv[id].vx(isend)=0.0;
      if (vel_dim >= 2) adv[id].vy(isend)=0.0;
      if (vel_dim == 3) adv[id].vz(isend)=0.0;
    } // end if (adv[id].x(isend) >=  fnx) (exited right side)
    else if (adv[id].x(isend) > fabsent2) {  //Tests for absent particle
      // Pack the array for passing west:
      // Test for the array being too small:
      if (nsend_west+n_per_particle > send_west.size()) {
	cout << "Error: "
	     << "insufficient space in the send_west array in distribution "  
	     << id << " rank " << mpi_rank << endl 
	     << "Need to increase part_pad" << endl;
	terminate(-1,"");
      }
      
      // We will adjust for the boundary condition here: (periodic only now)
      PTYPE x=adv[id].x(isend) + fnx; // This is done for rounding error reasons.
      if (x == fnx) x=fnx_minus_epsilon;
      
      // Test to see if the particle has crossed more than 1 domain.  
      if (x >= fnx || x < (PTYPE)0. ) {
	// If this happens a lot or the particle is really moving, terminate
	if (multi_cross++ >= adv[id].np*0.01)
	  terminate(-1,"Too many particles cross multiple domains");
	if  (x >= 2*fnx || x < -1*fnx ) 
	  terminate(-1,"Particle cross more than 2 domains");
	// This is just a rare event, so I will issue a warning and
	// adjust the particles position and velocity, 
	// Reduce the velocity so the particle just reaches the far side 
	// instead of overshooting:
	PTYPE x_orig=x-fnx-adv[id].vx(isend)*subcycle[id];
	PTYPE x_new=0;
	PTYPE vx_new=((x_new-fnx)-x_orig);
	cout << "Warning: "
	     << "Particle jumps more than 1 domain "
	     << "- adjusting to prevent. Dist: " 
	     << id << " rank " << mpi_rank << endl 
	     << "Old position and velocity " 
	     << x <<"  " << adv[id].vx(isend) <<endl
	     << "New position and velocity " 
	     << x_new <<"  " << vx_new <<endl;
	adv[id].vx(isend)=vx_new;
	x=x_new;
      }
      // only populate send_west if there is a processor to the west willing to accept
      if (!(boundary_type[0] != PERIODIC && subdomain.id_number == 0)){
        send_west(nsend_west++) = x;
        if(ndim > 1) send_west(nsend_west++) = adv[id].y(isend);
        if(ndim > 2) send_west(nsend_west++) = adv[id].z(isend);
        send_west(nsend_west++) = adv[id].vx(isend);
        if(vel_dim > 1) send_west(nsend_west++) = adv[id].vy(isend);
        if(vel_dim > 2) send_west(nsend_west++) = adv[id].vz(isend);
      }

      // Set the existing particle to an indicator value, 
      // outside all possible values.
      adv[id].x(isend)=fabsent2;
      if (NDIM >= 2) adv[id].y(isend)=fabsent2;
      if (NDIM >= 3) adv[id].z(isend)=fabsent2;
      // Ensure that it does not move in x direction:
      adv[id].vx(isend)=0.0;
      if(vel_dim >= 2) adv[id].vy(isend)=0.0;
      if(vel_dim == 3) adv[id].vz(isend)=0.0;
    }
  }      
  
  // send data to east and west processors
  // message-passing
  
  if (proc_east < mpi_np)
    MPI_Isend(send_east.address(), nsend_east, MPI_PTYPE, proc_east,
		2, MPI_COMM_WORLD, &request1);

  int nwest=0;
  PTYPEAVec recv_west;
  if (proc_west >= 0) {
    MPI_Probe(proc_west, 2, MPI_COMM_WORLD, &status1);
    MPI_Get_count(&status1, MPI_PTYPE, &nwest);
    int nwest_tmp[1]={nwest+1};
    recv_west.build_ArrayNd(nwest_tmp,0.0);
    MPI_Recv(recv_west.address(), nwest, MPI_PTYPE, proc_west,
	     2, MPI_COMM_WORLD, &status1);
  }

  MPI_Barrier(subdomain.neighbor_comm);

  if (proc_west >= 0) 
    MPI_Isend(send_west.address(), nsend_west, MPI_PTYPE, proc_west,
	      1, MPI_COMM_WORLD, &request2);

  int neast = 0;
  PTYPEAVec recv_east;
  if (proc_east < mpi_np){
    MPI_Probe(proc_east, 1, MPI_COMM_WORLD, &status2);
    MPI_Get_count(&status2, MPI_PTYPE, &neast);
    int neast_tmp[1]={neast+1};
    recv_east.build_ArrayNd(neast_tmp,0.0);
    MPI_Recv(recv_east.address(), neast, MPI_PTYPE, proc_east, 1, 
	     MPI_COMM_WORLD, &status2);
  }

  // make sure non-blocking sends have completed
  
  if (proc_west >= 0) MPI_Wait(&request2, &status2);
  if (proc_east < mpi_np) MPI_Wait(&request1, &status1);


  fill_particle_holes(recv_east, neast, adv[id], vel_dim);
  fill_particle_holes(recv_west, nwest, adv[id], vel_dim);
}
