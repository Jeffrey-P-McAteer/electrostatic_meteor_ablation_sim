USE_MPI = ON
NDIM = 3
#Define variables passed into program	
include ../make.inc

BIN = testAfft.x
OUT = log

LIBALL = -L /usr/local/IT/fftw-2.1.3/lib/aix5 -lrfftw_mpi -lfftw_mpi $(LIBS) 

$(OUT) : $(BIN) makefile ../make.inc
	poe $(BIN) -procs 4 > $(OUT) 

$(BIN) : $(CXXOBJS) $(CLASSLIB) makefile 
	$(CXX) $(CXXFLAGS) $(CXXOBJS) $(CLASSLIB) $(LIBALL) -o $(BIN)

$(CXXOBJS) : $(INCS) 
