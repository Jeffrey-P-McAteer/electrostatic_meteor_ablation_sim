#This is A GNU makefile - using GNU MAKE extensions and expecting the system
#variable HOSTTYPE to be defined

#Define variables passed into program	(first is default)

ifeq (scc, $(findstring scc,$(HOSTNAME)))
	include make.inc.buscc
else ifeq ($(TACC_SYSTEM),stampede)
        include make.inc.stampede
else ifeq ($(TACC_SYSTEM),stampede2)
	include make.inc.knl
else ifeq ($(HOSTNAME),docker)
	include make.inc.docker
else
	include make.inc
endif

#the list of source files
CXXSOURCES = $(wildcard *.cc)
CXXOBJS = $(CXXSOURCES:.cc=.o)
CSOURCES = $(wildcard *.c)
COBJS = $(CSOURCES:.c=.o)
FSOURCES = $(wildcard *.f)
FOBJS = $(FSOURCES:.f=.o)
INCS = $(wildcard *.h) $(wildcard classes/*.h)
CLASSDIR = classes
INCDIRS = -I../include -I$(CLASSDIR) -I. -I$(FFTWINCDIR) -I$(HDF5INCDIR)
ifeq ($(HAVE_PETSC),1)
   INCDIRS += -I$(TACC_PETSC_DIR)/$(PETSC_ARCH)/include -I$(TACC_PETSC_DIR)/include 
endif
ifeq ($(USE_P3DFFT),1)
   INCDIRS += -I$(P3DFFT_INC)
endif

LINUX_DIR = /usr/lib/gcc/x86_64-redhat-linux/4.4.7
# MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_COMPILE_FLAGS = $(shell $(CXX) --showme:compile)
#MPI_LINK_FLAGS = -L$(MPI_DIR) -L$(LINUX_DIR) -lmpi_f90 -lmpi_f77 -lgfortran -lm -lm -Wl,-rpath,$(MPI_DIR)/lib -lm -lm -lmpi_cxx -lstdc++ -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl
MPI_LINK_FLAGS = -L$(MPI_DIR) -L$(LINUX_DIR) -lmpi_f77 -lgfortran -lm -lm -Wl,-rpath,$(MPI_DIR)/lib -lm -lm -lmpi_cxx -lstdc++ -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl
make_files = makefile make.inc make.inc.* make.options

CXXFLAGS += $(INCDIRS)

#all: $(BIN)

#$(OUTDIR)/eppic.log: $(BIN) $(OUTDIR)/eppic.i $(OUTDIR)/eppic.tgz
#	cd $(RUNDIR)
#	$(TIMER) $(RUNPROGRAM) $(RUNFLAGS) $(BIN) $(OUTDIR)/eppic.i > $(OUTDIR)/eppic.log 
#	tail -40 $(OUTDIR)/eppic.log
#$(OUTDIR)/eppic.log: $(BIN) $(OUTDIR)/eppic.i $(OUTDIR)/eppic.tgz
#	cd $(RUNDIR)
#	$(TIMER) $(RUNPROGRAM) $(RUNFLAGS) $(BIN) $(OUTDIR)/eppic_test.i > $(OUTDIR)/eppic.log 
#	tail -40 $(OUTDIR)/eppic.log


#$(OUTDIR)/eppic.tgz: $(CXXOBJS) $(COBJS) $(FOBJS) $(CLASSLIB) $(make_files)
#	-/bin/rm -rf ../*/core* ../core* 
#	cd .. ;\
#	tar chf - src/*.h* src/*.c* src/makefile src/make.inc src/classes */*.com data/eppic.i bin |gzip - > src/$(OUTDIR)/eppic.tgz

$(BIN) : $(CXXOBJS) $(COBJS) $(FOBJS) $(CLASSLIB) $(make_files)
	-/bin/rm -f $(BINNAME)
	$(CXX) $(CXXFLAGS) $(CXXOBJS) $(COBJS) $(FOBJS) $(CLASSLIB) $(LIBS) -o $(BINNAME)
	cp $(BINNAME) $(BINDIR)/$(BIN)

tacc: $(BIN)
	cp $(BIN) $(SCRATCH)/eppic


#$(CLASSLIB): classes/makefile classes/*.cc classes/*.h
#	cd $(CLASSDIR); $(MAKE)

$(FFTLIB): 
	cd ../fftw-2.1/; csh ../bin/install-fftw.com

$(CXXOBJS) : $(INCS) $(make_files)

$(COBJS) : $(INCS) $(make_files)

$(FOBJS) : $(make_files)

clean:
	/bin/rm -rf *.o eppic.x* *.tmp restart/* core* lib* *.out
#	cd $(CLASSDIR); $(MAKE) clean

echo:
	@echo CPPFLAGS: $(CPPFLAGS)
	@echo CXXFLAGS: $(CXXFLAGS)
	@echo RTMFLAGS: $(RTMFLAGS)
	@echo FFTWLIBS: $(FFTWLIBS)
	@echo LIBS: $(LIBS)
#	@echo make_files: $(make_files)
