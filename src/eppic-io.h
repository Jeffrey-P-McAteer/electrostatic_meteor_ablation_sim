
// This file holds definitions to key I/O based routines used in EPPIC

#ifndef EPPIC_IO_H
#define EPPIC_IO_H

#include <hdf5.h>
#include "eppic-types.h"
#include "eppic-system.h"
#ifdef USE_MPI
#include "eppic-mpi.h"

void read_domains(FArrayND_ranged &array,
		  MPI_Comm across_comm,
		  char *path_name,
		  int nghost_pts[2]);

#endif

void terminate(int n, const char *message);

int output_FTarray(FILE *fname, FArrayND &A, eppic_system &sys,
		   FTYPE Afrac, FTYPE kmax=0.);

void output_array(FILE*, FArrayND &, int ngrid[NDIM],int nout_avg);

void output_array(FILE*, FArrayND_ranged &);

void output_array_h5(hid_t, char *, FArrayND &, int ngrid[NDIM],int nout_avg);
void output_array_h5(hid_t, char *, FArrayND_ranged &);

void output_collective_array_h5(hid_t, char *, FArrayND &, int ngrid[NDIM],int nout_avg);

FILE* fopensafe(char* filename, char* mode);

FILE* fopensafe(char* filename, char* mode, unsigned long skip);

// Not yet implemented:
//void infile_dictionary(char *infile, dictionary &sysInParams, 
//		       dictVec &distInParams);

void trim_strip_quote(std::string &instr);

#endif
