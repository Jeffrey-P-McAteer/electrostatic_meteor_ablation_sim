
/// Takes the solutions stored in the working arrays and returns them to 
/// the Efield phi member, as well as syncing domain boundary ghost cells. 
///
/// In 3D the working data is subdomain-decomposed, so it needs to be gathered
/// with in a domain. To simplify the conversion between 2D and 3D, the 
/// same array used in the 3D gather is used in the 2D part of the code. 
/// Since the 2D is not subdomain-decomposed (at this point, it is in fourier
/// space, but that's already taken care of else where), this "gathering" 
/// array, with prefix "full_", is simply defined using the same workspace as
/// the original working arrays. 

#include "efield_inject_tridiag.h" /* for field, type defs */
#include "eppic-mpi.h" /* for mpi, subdomain stuff */
#include "timerFile.h"

void efield_inject_tridiag_tophi(rfftw_many &working_array,
				 rfftw_many &boundary_working_array,
				 field &Efield)
{
  static timerFile timer;
  clock_t before_time;
  static bool firstEntry=true;
  if (firstEntry) {
    firstEntry = false;
    init_timerFile(timer,"efield_tophi_times.dat",64,0);
  }



  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  // periodic? 
  bool use_cyclic = false;
  // Commented out by Glenn 4/10/2018 because use_cyclic=true gives an error because
  // at 0 and np-1 processors there is no east and west processors respectively
  //if ((Efield.boundary[0] == periodic) || (Efield.boundary[1] == periodic)) use_cyclic = true;

  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

  // setup full working array
  before_time=times(&timer.refTimer);
#if NDIM==2
  // in 2D, this is all set, so just redefine full working array for
  // variable name consistancy from 2D to 3D

  ArrayNd_ranged<FTYPE,2> full_working_array, full_boundary_working_array;

  int array_size[2] = {working_array.data.size(0),
		       working_array.data.size(1)};
  
  full_working_array = working_array.data;
  full_boundary_working_array=boundary_working_array.data;
  
#else
  // in 3D, the data is subdomain-decomposed, so need an allgather, 
  // until there is actual subdomain-decomposition elsewhere in the code 
  // (ie PIC part of code)

  ArrayNd_ranged<FTYPE,3> full_working_array, full_boundary_working_array;
  int start_index[3] ={0,0,0};

  // get ranges on each processor

  int mpi_np_working=0,mpi_rank_working;
  MPI_Comm_size(working_array.comm, &mpi_np_working);
  MPI_Comm_rank(working_array.comm, &mpi_rank_working);

  if (mpi_np_working==1) { // then just need to copy
    
    int array_size[3] = {working_array.data.size(0),
			 working_array.data.size(1),
			 working_array.data.size(2)
    };
    
    full_working_array = working_array.data; 
    array_size[2]=boundary_working_array.data.size(3);
    full_boundary_working_array = boundary_working_array.data;

  } else {

    FTYPE* old_start_add;
    FTYPE* start_add;

    int array_size[3] = {ny,
			 working_array.nsize[1],
			 working_array.nsize[2]};

    full_working_array.build_ArrayNd_ranged(start_index,array_size);


    int xranges[mpi_np_working];// = new int[mpi_np_working];

    MPI_Allgather(&(working_array.local_nx),1,MPI_INT,
		  &xranges,1,MPI_INT,
		  working_array.comm);
    
    
    // gather data from each processor
    
    int recvcounts[mpi_np_working];
    int displs[mpi_np_working];
    
    recvcounts[0] = xranges[0]*array_size[1]*array_size[2];
    displs[0] = 0;
    for (int ip=1;ip<mpi_np_working;ip++) {
      recvcounts[ip] = xranges[ip]*array_size[1]*array_size[2]; 
      displs[ip] = recvcounts[ip-1]+displs[ip-1];
    }
    
    int old_start_index[3]={working_array.x_start,0,0};

    if (working_array.allocatable()) {
      old_start_add = working_array.data.address(old_start_index);
    }
    start_add = full_working_array.address(start_index);
    
    MPI_Allgatherv(old_start_add,
		   xranges[mpi_rank_working]*array_size[1]*array_size[2],
		   MPI_FTYPE,
		   start_add,
		   recvcounts,
		   displs,
		   MPI_FTYPE,
		   working_array.comm);
    
    
    array_size[2]=boundary_working_array.nsize[2];
    full_boundary_working_array.build_ArrayNd_ranged(start_index,array_size);
    
  // gather data on each processor

    recvcounts[0] = (recvcounts[0]*3)/working_array.nsize[2];
    for (int ip=1;ip<mpi_np_working;ip++) {
      recvcounts[ip] = (recvcounts[ip]*3)/working_array.nsize[2];
      displs[ip] = recvcounts[ip-1]+displs[ip-1];
    }
    
    if (boundary_working_array.allocatable()) {
      old_start_add = boundary_working_array.data.address(old_start_index);
    }

    MPI_Allgatherv(old_start_add,
		   xranges[mpi_rank_working]*array_size[1]*array_size[2],
		   MPI_FTYPE,
		   full_boundary_working_array.address(start_index),
		   recvcounts,displs,MPI_FTYPE,
		   boundary_working_array.comm);
    
  } // else mpi_np_working != 1
#endif
  timer.times["init"] += times(&timer.refTimer)-before_time;

  before_time=times(&timer.refTimer);
  // phi dc calculation 
  FTYPE phi_avg=0.;
  for(int ix=0;ix<nx;ix++) 
    for(int iy=0;iy<ny;iy++)
      for(int iz=0;iz<nz;iz++){
	phi_avg+=full_working_array(RFFTW_INDICIES(ix,iy,iz));
	phi(INDICIES(ix,iy,iz))=
	  full_working_array(RFFTW_INDICIES(ix,iy,iz));
      }

  if (!use_cyclic) {
    if (subdomain.id_number==0) {
      for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++) {
	  phi_avg+=full_boundary_working_array(RFFTW_INDICIES(0,iy,iz));
	  phi(INDICIES(-1,iy,iz))=
	    full_boundary_working_array(RFFTW_INDICIES(0,iy,iz));
	}
    }
    if (subdomain.id_number==nsubdomains-1) {
      for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++) {
	  phi_avg+=
	    full_boundary_working_array(RFFTW_INDICIES(1,iy,iz))
	    +full_boundary_working_array(RFFTW_INDICIES(2,iy,iz));
	  phi(INDICIES(nx,iy,iz))=
	    full_boundary_working_array(RFFTW_INDICIES(1,iy,iz));
	  phi(INDICIES(nx+1,iy,iz))=
	    full_boundary_working_array(RFFTW_INDICIES(2,iy,iz));
	}
    }
  }

  timer.times["phidcCalc"] += times(&timer.refTimer)-before_time;
  if (nsubdomains>1) {
    
    before_time=times(&timer.refTimer);
    FTYPE phi_avg_all=0;
    MPI_Allreduce(&phi_avg,&phi_avg_all,1,MPI_FTYPE,MPI_SUM,
		  subdomain.neighbor_comm);
    phi_avg = phi_avg_all;
    timer.times["phidcSync"] += times(&timer.refTimer)-before_time;

  }
  // The division of ny is for normalization of the IFFT.
  phi_avg/=ny;
#if NDIM==3
  phi_avg/=nz;
#endif

  if (use_cyclic)phi_avg/=((nx*nsubdomains));
  else phi_avg/=((nx*nsubdomains+3));

  int ix_start=0;
  int ix_end=nx;
  if (!use_cyclic){
    if(subdomain.id_number==0)
      ix_start=-1;
    if(subdomain.id_number==nsubdomains-1) 
      ix_end=nx+2;
  }

  // force phi_avg=0 here; ideally zeroing the average should not be needed
  before_time=times(&timer.refTimer);
  phi_avg=0;
  for(int ix=ix_start;ix<ix_end;ix++) 
    for(int iy=0;iy<ny;iy++) 
      for(int iz=0;iz<nz;iz++) {
	phi(INDICIES(ix,iy,iz))=(phi(INDICIES(ix,iy,iz))-phi_avg)/nz/ny;
      }
  timer.times["phiset"] += times(&timer.refTimer)-before_time;

  MPI_Request request_west_send, request_east_send;
  MPI_Request request_west_recv, request_east_recv;
  MPI_Status status_west, status_east;

  /* pass west */
  /* if first domain, set to zero, no pass */
  before_time=times(&timer.refTimer);
  int size_send=phi.length()/phi.size(0);
  const int index0[]={INDICIES(0,0,0)};
  if ((subdomain.id_number > 0) || use_cyclic)
    MPI_Isend(phi.address(index0),size_send*phix_guard_size[1],MPI_FTYPE,
	      proc_west,2,MPI_COMM_WORLD,
	      &request_west_send);

  /* if last domain, no receive */
  const int index1[]={INDICIES(nx,0,0)};
  if ((subdomain.id_number < nsubdomains-1) || use_cyclic)
    MPI_Irecv(phi.address(index1),size_send*phix_guard_size[1],MPI_FTYPE,
	      proc_east,2,MPI_COMM_WORLD,
	      &request_west_recv);

  if ((subdomain.id_number < nsubdomains-1) || use_cyclic)
    MPI_Wait(&request_west_recv,&status_west);
  if ((subdomain.id_number > 0) || use_cyclic)
    MPI_Wait(&request_west_send,&status_west);  

  /* pass east */
  /* if last domain, set to zero, no pass */
  const int index2[]={INDICIES(nx-phix_guard_size[0],0,0)};
  if ((subdomain.id_number < nsubdomains-1) || use_cyclic)
    MPI_Isend(phi.address(index2),size_send*phix_guard_size[0],MPI_FTYPE,
	      proc_east,3,MPI_COMM_WORLD,
	      &request_east_send);

  /* if first domain, no receive */
  const int index3[]={INDICIES(-(int(phix_guard_size[0])),0,0)};
  if ((subdomain.id_number > 0) || use_cyclic)
    MPI_Irecv(phi.address(index3),size_send*phix_guard_size[0],MPI_FTYPE,
	      proc_west,3,MPI_COMM_WORLD,
	      &request_east_recv);

  if ((subdomain.id_number > 0) || use_cyclic)
    MPI_Wait(&request_east_recv,&status_east);
  if ((subdomain.id_number < nsubdomains - 1) || use_cyclic)
    MPI_Wait(&request_east_send,&status_east);

  timer.times["bordercomm"] += times(&timer.refTimer)-before_time;

  flush_timerFile(timer,-1);

}
