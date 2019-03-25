
include make.options

#-----------------------------------------------------------------------------
# Contains definitions customized for this system
#-----------------------------------------------------------------------------

# library locations, compiler/program definitions
CPPFLAGS = -D F_PRECISION=$(F_PRECISION) \
	-D P_PRECISION=$(P_PRECISION) \
	-D NDIM=$(NDIM) 

export MPI_COMPILER=gnu
CXX = g++
CXXCPP = 
CXXFLAGS_OPT = -O3
CXXFLAGS_DEBUG = -O0
FFLAGS = -O4
FFLAGS_DEBUG = -O0
FFTWLIBDIR = /projectnb/eregion/yannpaul/fftw2.1.5/install/lib
FFTWINCDIR = /projectnb/eregion/yannpaul/fftw2.1.5/install/include
MPICXX = mpic+++
TIMER = time
HDF5LIBDIR = /usr/lib64
HDF5INCDIR = /usr/include
HDF5LIBS = -lhdf5_hl -lhdf5

# I don't know what this is for?
CLASSLIB = 
BIN = $(BINDIR)/$(BINNAME)



FFTWLIBS = 
ifeq ($(USE_FFTW3),1) 
	ifeq ($(F_PRECISION),1)
		FFTLIBS = -lfftw3f
	else		
		FFTLIBS = -lfftw3
	endif
else
	ifeq ($(F_PRECISION),1)
		FFTLIBS = -lsrfftw libsfftw
	else		
		FFTLIBS = -lrfftw -lfftw
	endif
endif

ifeq ($(USE_DOMAINS),1)
	FFTWLIBS = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
	CPPFLAGS += -D USE_DOMAINS=1
endif

ifeq ($(DEBUG),1) 
	CXXFLAGS = -g $(CXXFLAGS_DEBUG)
	CPPFLAGS += -D DEBUG=1 
		-D CHECK=1
else
	CXXFLAGS = $(CXXFLAGS_OPT)
endif

ifeq ($(PROFILING),1)
	CXXFLAGS += -pg
endif

ifeq ($(EPPIC_FFTW_USE_D_PREFIX),1)
	CPPFLAGS += -D EPPIC_FFTW_USE_D_PREFIX=1
endif

ifeq ($(USE_MPI),1)
	CXX = $(MPICXX)
	CPPFLAGS += -D USE_MPI=1
endif



LIBS = -L $(FFTWLIBDIR) $(FFTWLIBS) -lm $(HDF5LIBS)

