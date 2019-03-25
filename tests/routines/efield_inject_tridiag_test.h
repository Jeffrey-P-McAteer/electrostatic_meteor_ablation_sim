
/// all globals and includes for efield_inject_tridiag_test program

using namespace std;
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "eppic-system.h"
#include "eppic-efield.h" /* efield routines, includes eppic types */
#include "eppic-math.h" /* for Sqr function */
#include "eppic-mpi.h" /* for domain_decomp routines */
#include "eppic-io.h"
/* globals pertaining to above header */
ArrayNd<MPI_Comm,1> comm_across_subdomains;
ArrayNd<MPI_Comm,1> comm_within_subdomain;
int proc_west;
int proc_east;
int boundary_type;
int mpi_rank=0,mpi_np=1;
int nsubdomains=1;
subdomain_traits subdomain;

FTYPE Ex0_external=0.; 
FTYPE Ey0_external=0.; 
#include "eppic-times.h" /* for efield timing routines */
tms times_buf;
clock_t efield_time=0;
clock_t first_time=0;

#ifndef PI
#define PI  3.1415926535897932385
#endif
