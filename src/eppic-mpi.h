

#ifndef EPPIC_MPI_H
#define EPPIC_MPI_H
#include "eppic-types.h"


#ifdef USE_MPI


// We do not want MPI C++ bindings
#define MPI_NO_CPPBIND
#define MPI_SKIP_MPICXX
#define MPICH_IGNORE_CXX_SEEK

//dls//////////////////////////////////////
#define NPMAX 55000
#define ABSENT -999.0
//dls end//////////////////////////////////

#include <mpi.h>




#if F_PRECISION == 2
const MPI_Datatype MPI_FTYPE = MPI_DOUBLE;
#else
const MPI_Datatype MPI_FTYPE = MPI_FLOAT;
#endif

#if P_PRECISION == 2
const MPI_Datatype MPI_PTYPE = MPI_DOUBLE;
#define PTYPE_MPI MPI_DOUBLE
#else
const MPI_Datatype MPI_PTYPE = MPI_FLOAT;
#define PTYPE_MPI MPI_FLOAT
#define PTYPE_MPI MPI_FLOAT
#endif

#define MPI_OTYPE MPI_FLOAT


// domain decomposition
extern ArrayNd<MPI_Comm,1> comm_within_subdomain;
extern ArrayNd<MPI_Comm,1> comm_across_subdomains;

typedef struct {
  int id_number, np; // The subdomain id number and number of processors.
  int rank,root;          // rank within subdomain, root of subdomain
  MPI_Comm internal_comm; // Communicator for within the subdomain and for processors
  MPI_Comm neighbor_comm; // Communicator to all processors with the same subdomain.rank
} subdomain_traits;

extern subdomain_traits subdomain;

void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);
void mpi_check(int error_code);

//set the following to make adjacent subdomains have adjacent mpi ranks,
// false will make processors within a subdomain have adajacent ranks.
#if NDIM == 3 && USE_P3DFFT ==1 || HAVE_PETSC == 1
const bool SUBDOMAINS_ADJACENT = false;
#else
const bool SUBDOMAINS_ADJACENT = true;
#endif

void domain_decomp();

#endif

// Parallel processing variables to be known globally 
extern int mpi_rank, mpi_np; 

extern int nsubdomains;      // number of subdomains
extern int proc_east;        // mpi_rank of processor to the east
extern int proc_west;        // mpi_rank of processor to the west

#endif
