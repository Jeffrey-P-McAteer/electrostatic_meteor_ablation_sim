

#ifndef USE_MPI
#include "srcconfig.h"
#endif
#ifdef EPPIC_FFTW_USE_D_PREFIX 
#include "drfftw_mpi.h"
#ifdef STAMPEDE2_FFTW
#include "fftw3-mpi.h"
#endif

#else
#include "rfftw_mpi.h"

#endif
