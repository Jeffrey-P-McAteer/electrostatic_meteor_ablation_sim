#include "eppic-mpi.h"
#include "eppic-system.h"

#if HAVE_PETSC
#include "petsc.h"
#endif

#ifdef USE_DOMAINS


void domain_decomp()
{

  MPI_Group group_world, group_within_subdomain, group_across_subdomains;
  int i, j;

  /* UNDERDEVELOPMENT ...
  // Applying unequal numbers of processors in each domain requires 
  // calculating the number per domain:
  if (ndomain_mult.max() > 1) {
    // We need to sum the requested # of subdomains
    int ndomain_min=nsubdomain;
    for (int i=0; i<ndomain_index_max; i++) ndomain_min += (ndomain_mult-1);
    // Is the number of processors evenly divisible by ndomain_min:
    if (mpi_np%ndomain_min != 0) 
      terminate (-1,"The min number of domains requested not evenly divisible by mpi_np");
    // Set the number of processors in each subdomain
    int min_subdomain_np=mpi_np/ndomain_min;

    subdomain.np=min_subdomain_np;
    for (int i=0; i< ndomain_index_max, i++) 
      if m

    // Set the subdomain id num
    subdomain.id_number=mpi_rank/ndomain_min

  END NEW DEVELOPMENT SECTION */ 

  // PETSC Communicators written on 1/30/2018 by Glenn
  // MPI processors split up z, y, x but Petsc splits up x, y, z
  // (Note we don't split up domains in the z direction which is
  // why we set the z variable to 0).
  // We need to switch the PETSC_COMM_WORLD communicator for the 
  // field solver
  /*#if HAVE_PETSC
  MPI_Comm NewComm;
  int NewRank, x, y, z, x_procs, y_procs, z_procs;

  x_procs = nsubdomains;
  y_procs = mpi_np/nsubdomains;
  x = mpi_rank / x_procs;
  y = mpi_rank % x_procs;
  z = 0;

  // set new rank according to PETSc ordering:
  NewRank = z*y_procs*x_procs + y*x_procs + x;

  printf("xprocs: %d, yprocs: %d, zprocs %d\n", x_procs, y_procs, z_procs);
  printf("Old rank: %d, New Rank: %d\n", mpi_rank, NewRank);
  // create communicator with new ranks according to PETSc ordering:
  MPI_Comm_split(PETSC_COMM_WORLD, 1, NewRank, &NewComm);
    
  // override the default communicator (was MPI_COMM_WORLD as default)
  PETSC_COMM_WORLD = NewComm;
  #endif*/

  // array of communicators, one element for each subdomain
  comm_within_subdomain = ArrayNd<MPI_Comm,1>(nsubdomains);

  // array of communicators, one element for each processor on a given subdomain
  subdomain.np = mpi_np/nsubdomains;

  comm_across_subdomains = ArrayNd<MPI_Comm,1>(subdomain.np);

  // group containing all global procs
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  // create intracommunicator comm_within_subdomain[i] for each subdomain i
  ArrayNd<int,1> ranks_within(subdomain.np);
  for(i=0; i<nsubdomains; i++){

    if (SUBDOMAINS_ADJACENT) {
      for(j=0; j<subdomain.np; j++) ranks_within[j] = j*nsubdomains + i;
    } else {
      for(j=0; j<subdomain.np; j++) ranks_within[j] = i*subdomain.np + j;
    }

    MPI_Group_incl(group_world, subdomain.np,
		   ranks_within.address(), &group_within_subdomain);
    MPI_Comm_create(MPI_COMM_WORLD, group_within_subdomain,
		    (comm_within_subdomain.address(i)) );
  }

  // create intercommunicator comm_across_subdomains[j]
  // containing one processor in each subdomain
  ArrayNd<int,1> ranks_across(nsubdomains);
  for(j=0; j<subdomain.np; j++){

    if (SUBDOMAINS_ADJACENT) {
      for(i=0; i<nsubdomains; i++) ranks_across[i] = j*nsubdomains + i;
    } else {
      for(i=0; i<nsubdomains; i++) ranks_across[i] = i*subdomain.np + j;
    }

    MPI_Group_incl(group_world, nsubdomains,
		   ranks_across.address(), &group_across_subdomains);
    MPI_Comm_create(MPI_COMM_WORLD, group_across_subdomains,
		    comm_across_subdomains.address(j));
  }

  // local subdomain identifiers
  if (SUBDOMAINS_ADJACENT) {
    subdomain.id_number = mpi_rank%nsubdomains; // The subdomain id number 
    subdomain.rank = mpi_rank/nsubdomains; //The number of processors on the subdomain
  } else {
    subdomain.id_number = mpi_rank/subdomain.np; // The subdomain id number 
    subdomain.rank = mpi_rank%subdomain.np; //The number of processors on the subdomain
  }
  subdomain.root = 0; //subdomain rank of subdomain processor treated as root
  subdomain.internal_comm=comm_within_subdomain[subdomain.id_number]; // Communicator for within the subdomain and for processors
  subdomain.neighbor_comm=comm_across_subdomains[subdomain.rank]; // Communicator to all processors with the same subdomain.rank

  // identify east and west processors

  if (SUBDOMAINS_ADJACENT) {
    proc_east = (mpi_rank + 1);
    proc_west = (mpi_rank - 1);
    if (true){//(boundary_type[0] == PERIODIC || boundary_type[0] == OPEN ) {// Incorrect for OPEN but for now...
      if (subdomain.id_number == nsubdomains-1) proc_east=mpi_rank-(nsubdomains-1);
      if (subdomain.id_number == 0) proc_west=mpi_rank+(nsubdomains-1);
    }
  } else {
    if (true){//(boundary_type[0] == PERIODIC || boundary_type[0] == OPEN) {// Incorrect for OPEN but for now...
      proc_east = (mpi_rank + subdomain.np)%mpi_np;
      proc_west = (mpi_rank + mpi_np - subdomain.np)%mpi_np;
    } else {
      proc_east = (mpi_rank + subdomain.np);
      proc_west = (mpi_rank - subdomain.np);
    }
  }
}

#endif
