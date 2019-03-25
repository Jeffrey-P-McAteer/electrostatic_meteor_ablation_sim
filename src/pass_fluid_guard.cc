//----------------------------------------------------

// Copy fluid guard cells from neighboring domains
// use non-blocking sends

//----------------------------------------------------

#include "eppic.h"
#include "eppic-mpi.h"

void pass_fluid_guard(fluid &fspecie,int id)
{
  
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, 
		   const int nx_guard[]);
  //  int n_per_particle=ndim+vel_dim;
  //  int max_absent=(int) (adv.np*(part_pad[id]-1)/part_pad[id]);
  //  PTYPE fnx = (PTYPE) nx;

  //  PTYPE smallest_float_below(PTYPE x);
  //  static PTYPE fnx_minus_epsilon=smallest_float_below(fnx);

  const int den_guard_size[] = {denx_guard_size,denx_guard_size};
  const int vel_guard_size[] = {velx_guard_size,velx_guard_size};

  pass_guards(fspecie.den,den_guard_size);
  pass_guards(fspecie.vx,vel_guard_size);
  if (vel_dim[id] >= 2) pass_guards(fspecie.vy,vel_guard_size);
  if (vel_dim[id] == 3) pass_guards(fspecie.vz,vel_guard_size);
  
 }


  /*
  fluid send_west,send_east,recv_east,recv_west;
  // densities
  send_east.den = FArrayND_ranged(allstart,denx_guard_size);
  recv_west.den = FArrayND_ranged(allstart,denx_guard_size);
  send_west.den = FArrayND_ranged(allstart,denx_guard_size);
  recv_east.den = FArrayND_ranged(allstart,denx_guard_size);

  send_east.vx =  FArrayND_ranged(allstart,velx_guard_size);
  recv_west.vx =  FArrayND_ranged(allstart,velx_guard_size);
  send_west.vx =  FArrayND_ranged(allstart,velx_guard_size);
  recv_east.vx =  FArrayND_ranged(allstart,velx_guard_size);

  send_east.vy =  FArrayND_ranged(allstart,velx_guard_size);
  recv_west.vy =  FArrayND_ranged(allstart,velx_guard_size);
  send_west.vy =  FArrayND_ranged(allstart,velx_guard_size);
  recv_east.vy =  FArrayND_ranged(allstart,velx_guard_size);

  send_east.vz =  FArrayND_ranged(allstart,velx_guard_size);
  recv_west.vz =  FArrayND_ranged(allstart,velx_guard_size);
  send_west.vz =  FArrayND_ranged(allstart,velx_guard_size);
  recv_east.vz =  FArrayND_ranged(allstart,velx_guard_size);


  // den from west to east
  for(ix=0;ix<denx_guard_size;ix++) {
#if NDIM > 1
    for(iy=0;iy<ny;iy++) {
#endif
#if NDIM > 2
      for(iz=0;iz<nz;iz++) {
#endif
	send_west.den(INDICIES(ix,iy,iz)) = fspecie.den(INDICIES(ix,iy,iz))
#if NDIM > 2
      } //for iz
#endif
#if NDIM > 1
    } // for iy
#endif
  } // for ix

  // den from east to west
  for(ix=nx-denx_guard_size;ix<nx;ix++) {
#if NDIM > 1
    for(iy=0;iy<ny;iy++) {
#endif
#if NDIM > 2
      for(iz=0;iz<nz;iz++) {
#endif
	send_east.den(INDICIES(ix-nx+denx_guard_size,iy,iz)) = 
	  fspecie.den(INDICIES(ix,iy,iz));
#if NDIM > 2
      } //for iz
#endif
#if NDIM > 1
    } // for iy
#endif
  } // for ix

  // v from east to west
  for(ix=0;ix<velx_guard_size;ix++) {
#if NDIM > 1
    for(iy=0;iy<ny;iy++) {
#endif
#if NDIM > 2
      for(iz=0;iz<nz;iz++) {
#endif
	send_east.vx(INDICIES(ix,iy,iz)) = fspecie.vx(INDICIES(ix,iy,iz))
	send_east.vy(INDICIES(ix,iy,iz)) = fspecie.vy(INDICIES(ix,iy,iz))
	send_east.vz(INDICIES(ix,iy,iz)) = fspecie.vz(INDICIES(ix,iy,iz))
#if NDIM > 2
      } //for iz
#endif
#if NDIM > 1
    } // for iy
#endif
  } // for ix

  // v from east to west
  for(ix=nx-velx_guard_size;ix<nx;ix++) {
#if NDIM > 1
    for(iy=0;iy<ny;iy++) {
#endif
#if NDIM > 2
      for(iz=0;iz<nz;iz++) {
#endif
	send_east.vx(INDICIES(ix-nx+velx_guard_size,iy,iz)) = 
	  fspecie.vx(INDICIES(ix,iy,iz))
	send_east.vy(INDICIES(ix-nx+velx_guard_size,iy,iz)) = 
	  fspecie.vy(INDICIES(ix,iy,iz))
	send_east.vz(INDICIES(ix-nx+velx_guard_size,iy,iz)) = 
	  fspecie.vz(INDICIES(ix,iy,iz))
#if NDIM > 2
      } //for iz
#endif
#if NDIM > 1
    } // for iy
#endif
  } // for ix


  MPI_Request request1, request2;
  MPI_Status  status1,  status2;
  
  // Loop through absent array
  int nsend_east=0;
  int nsend_west=0;

  // send data to east and west processors
  
  // message-passing

  MPI_Isend(send_east.address(), nsend_east, MPI_PTYPE, proc_east,
            2, MPI_COMM_WORLD, &request1);

  MPI_Recv(recv_west.address(), max_absent*n_per_particle, MPI_PTYPE, proc_west,
            2, MPI_COMM_WORLD, &status1);


  MPI_Isend(send_west.address(), nsend_west, MPI_PTYPE, proc_west,
            1, MPI_COMM_WORLD, &request2);

  MPI_Recv(recv_east.address(), max_absent*n_per_particle, MPI_PTYPE, proc_east,
            1, MPI_COMM_WORLD, &status2);

  // determine how much data was received from 
  // adjacent processor

  //  int neast, nwest;
  //  MPI_Get_count(&status1, MPI_PTYPE, &nwest);
  //  MPI_Get_count(&status2, MPI_PTYPE, &neast);

  
  // make sure non-blocking sends have completed

  MPI_Wait(&request1, &status1);
  MPI_Wait(&request2, &status2);
  */
