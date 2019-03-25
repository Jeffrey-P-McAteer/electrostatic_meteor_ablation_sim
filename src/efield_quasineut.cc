 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solves for phi, the electric potential, using the quasi-neutral 
     form of the inertialess electron momentum equation, in which the ion 
     density and flux can replace the analogous electron terms and all 
     other parameters are from the electrons' definition. 

     This equation resembles a modified Poisson's equation. 

     We use the PETSc library, which allows for distributed solutions to
     linear (used in this case, ie Ax=b), non-linear, and matrix free 
     problems. 

     A detailed description of how this code works is given in:
     JGR, VOL. 101, NO. A8, PAGES 17,273-17,286, AUGUST 1, 1996

     The following paper is another example using this algorithm:
     GRL, VOL. 22, NO. 4, PAGES 353-356, FEBRUARY 15, 1995
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#if 0 // I think this code is obsolete - may (14Mar2018)
#if HAVE_PETSC
#include "petscksp.h" // includes other header files, eg petscmat.h ...

/*
  ONLY SETUP FOR 2D, FOR NOW.... 
  And Bz (not Bx)

  Consider putting all the BC stuff here. Make n[x,y,z]_guards large enough for
  the necessary derivatives (see getPetsc[L,R]HS) and write appropriate BC data
  to those cells before passing to efield_quasineut_solver. Then, in getPetsc[L,R]HS,
  just access the guard cells explicitly.
*/

#include "efield_quasineut.h"

#undef __FUNCT__
#define __FUNCT__ "efield_quasineut"
void efield_quasineut(field &Efield,particle_dist *pic,fluid *fspecie)
{
  
  /* Error handling */
  PetscErrorCode perr=0;
  int merr=0;
  int status=0;
  MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);

  /* Local function declarations */
  void mpi_check(int error_code);
  void qden_divG(FArrayND &qden,FArrayND &divG,int id,
		 particle_dist &pic,fluid &fspecie,FTYPE scaler=1.0);
  void periodic_filter(FArrayND_ranged &array,const int nx_guard[]);
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  void pass_sum_guard(FArrayND &,int nx);
  void write_local_bin(FArrayND_ranged array, const char *name, int *wxguards);
  void write_local_bin(FArrayND array, const char *name);
  int efield_quasineut_solver(field &,FArrayND_ranged,FArrayND_ranged,fluid);

  /* Local arrays */
  const int nx_guard[] = {1,1}; // -->Also ny_guard[] and nz_guard[] for non-periodic BC?
  int wxguards[] = {0,0}; // For write_local_bin()
  int tmp_size[]={INDICIES(nx+xguard_size,ny,nz)};
  int ghosted_size[]={INDICIES(nx+nx_guard[0]+nx_guard[1],ny,nz)};
  int ghosted_start[]={INDICIES(0-nx_guard[0],0,0)};
  FArrayND divG_tmp,qden_tmp;
  FArrayND divG,qden,work;
  FArrayND_ranged divG_ranged,qden_ranged;
  divG_tmp = FArrayND(tmp_size);
  qden_tmp = FArrayND(tmp_size);
  divG = FArrayND(tmp_size);
  qden = FArrayND(tmp_size);
  work = FArrayND(tmp_size);
  divG_ranged = FArrayND_ranged(ghosted_start,ghosted_size);
  qden_ranged = FArrayND_ranged(ghosted_start,ghosted_size);
  
  /* Loop through each species: gather PIC */
  qden = 0.0;
  divG = 0.0;
  for (int id=0;id<ndist;id++) {
    if (id != electron_dist) {
      if (method[id] >= 0) { // Need to add subcycling (cf. charges.cc)
	qden_divG(qden_tmp,divG_tmp,id,pic[id],fspecie[id],1.0);
	qden += qden_tmp;
	divG += divG_tmp;
      }
    }
  }

  /* Average arrays within subdomains */
  // quasineutral density
  work = 0.0;
  merr=MPI_Allreduce((void*)&(qden(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     qden.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  mpi_check(merr);
  work /= subdomain.np;
  qden = work;
  // flux divergence
  work = 0.0;
  merr=MPI_Allreduce((void*)&(divG(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     divG.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  mpi_check(merr);
  work /= subdomain.np;
  divG = work;

  /* Loop through each species: copy fluids */
  for (int id=0;id<ndist;id++) {
    if (id != electron_dist) {
      if (method[id] < 0) { // Need to add subcycling (cf. charges.cc)
	qden_divG(qden_tmp,divG_tmp,id,pic[id],fspecie[id],1.0);
	for (int ix=0;ix<nx+xguard_size;ix++) {
	  for (int iy=0;iy<ny;iy++) {
	    for (int iz=0;iz<nz;iz++) {
	      qden(INDICIES(ix,iy,iz)) += qden_tmp(INDICIES(ix,iy,iz));
	      divG(INDICIES(ix,iy,iz)) += divG_tmp(INDICIES(ix,iy,iz));
	    }
	  }
	}
      }
    }
  }

  /* Pass pad cells around and sum */
  pass_sum_guard(qden,nx);
  pass_sum_guard(divG,nx);

  /* Write non-ghosted data to ghosted arrays */
  qden_ranged = 0.0;
  divG_ranged = 0.0;
  for (int ix=0;ix<nx;ix++) {
    for (int iy=0;iy<ny;iy++) {
      for (int iz=0;iz<nz;iz++) {
	qden_ranged(INDICIES(ix,iy,iz)) = qden(INDICIES(ix,iy,iz));
	divG_ranged(INDICIES(ix,iy,iz)) = divG(INDICIES(ix,iy,iz));
      }
    }
  }

  /* Take care of global BC in ghosted arrays */
  if (boundary_type[0] == PERIODIC) { // Calls pass_guards internally
    periodic_filter(divG_ranged,nx_guard);
    periodic_filter(qden_ranged,nx_guard);
  } else { // Add equivalent BC for divG??
    pass_guards(divG_ranged,nx_guard);
    pass_guards(qden_ranged,nx_guard);
    if (subdomain.id_number == 0) {
      for (int id=0;id<ndist;id++) {
  	if (id != electron_dist) {
  	  FTYPE n0dist_id = 0;
  	  if (method[id]>=0) n0dist_id = pic[id].n0avg;
  	  else n0dist_id = fspecie[id].n0;
  	  for (int iy=0;iy<ny;iy++) {
  	    for (int iz=0;iz<nz;iz++) {
  	      qden(INDICIES(-1,iy,iz)) += n0dist_id;
  	    }
  	  }
  	}
      }
    }
    if (subdomain.id_number == nsubdomains-1) {
      for (int id=0;id<ndist;id++) {
  	if (id != electron_dist) {
  	  FTYPE n0dist_id = 0;
  	  if (method[id]>=0) n0dist_id = pic[id].n0avg;
  	  else n0dist_id = fspecie[id].n0;
  	  for (int iy=0;iy<ny;iy++) {
  	    for (int iz=0;iz<nz;iz++) {
  	      qden(INDICIES(nx,iy,iz)) += n0dist_id;
  	    }
  	  }
  	}
      }
    }
  }

  /* Pass to PETSc solver */
  int efield_quasineut_solver(field &Efield,FArrayND_ranged qden,FArrayND_ranged divG);
  status = efield_quasineut_solver(Efield,qden_ranged,divG_ranged);
  if (status != 0) {
    if (mpi_rank == 0) printf("WARNING: efield_quasineutral_solver returned status = %d at step = %d\n",
			      status,it);
  }

  /* Reset MPI error handler */
  MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL);

}

#else
/* If not using PETSc, here is a dummy routine */
void efield_quasineut(field &Efield,particle_dist *pic,fluid *fspecie) 
{ return; }


#endif

#endif // Don't compile in routine (see comment at top)
