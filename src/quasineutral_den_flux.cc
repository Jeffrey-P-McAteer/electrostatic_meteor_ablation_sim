#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#if USE_QN
void quasineutral_den_flux(FArrayND_ranged &qden,
			   INDICIES(FArrayND_ranged &Gx,
				    FArrayND_ranged &Gy,
				    FArrayND_ranged &Gz),
			   particle_dist *pic,fluid *fspecie,
			   int it)
{

  // !!! SUBCYCLING (see charges.cc)

  // printf("[%d] %s:%d it = %d\n",mpi_rank,__func__,__LINE__,it);

  /* Let's time the amount of time in this routine */
  clock_t start_time=times(&times_buf);

  /* Error handling */
  int merr=0;
  int status=0;

  /* Local function declarations */
  void mpi_check(int error_code);
  void qden_flux(FArrayND &qden,INDICIES(FArrayND &Gx, \
					 FArrayND &Gy, \
					 FArrayND &Gz),
		 int id,
		 particle_dist *pic,fluid *fspecie,FTYPE scaler);
  void periodic_filter(FArrayND_ranged &array,const int nx_guard[]);
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  void background_dust(FArrayND &, FTYPE);
  void write_local_bin(FArrayND_ranged array, const char *name, int *wxguards);
  void write_local_bin(FArrayND array, const char *name);

  /* Local arrays */
  const int nx_guard[] = {1,1};
  int wxguards[] = {0,0}; // For write_local_bin()
  int tmp_size[]={INDICIES(nx+xguard_size,ny,nz)};

  FArrayND Gx_tmp,qden_tmp;
  FArrayND Gx_loc,qden_loc,work;
  Gx_tmp = FArrayND(tmp_size);
  Gx_loc = FArrayND(tmp_size);
  qden_tmp = FArrayND(tmp_size);
  qden_loc = FArrayND(tmp_size);
  work = FArrayND(tmp_size);
#if NDIM > 1
  FArrayND Gy_tmp;
  FArrayND Gy_loc;
  Gy_tmp = FArrayND(tmp_size);
  Gy_loc = FArrayND(tmp_size);
#if NDIM > 2
  FArrayND Gz_tmp;
  FArrayND Gz_loc;
  Gz_tmp = FArrayND(tmp_size);
  Gz_loc = FArrayND(tmp_size);
#endif
#endif

  /* Loop through each species: gather PIC */
  qden_loc = 0.0;
  Gx_loc = 0.0;
#if NDIM > 1
  Gy_loc = 0.0;
#if NDIM > 2
  Gz_loc = 0.0;
#endif
#endif
  for (int id=0;id<ndist;id++) {
    if (id != electron_dist) {
      if (method[id] >= 0) { // Need to add subcycling (cf. charges.cc)
	qden_flux(qden_tmp,INDICIES(Gx_tmp,Gy_tmp,Gz_tmp),
		  id,pic,fspecie,1.0);
	qden_loc += qden_tmp;
	Gx_loc += Gx_tmp;
#if NDIM > 1
	Gy_loc += Gy_tmp;
#if NDIM > 2
	Gz_loc += Gz_tmp;
#endif
#endif
      }
    }
  }

  /* Average arrays within subdomains */
  // quasineutral density
  work = 0.0;
  merr=MPI_Allreduce((void*)&(qden_loc(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     qden_loc.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  work /= subdomain.np;
  qden_loc = work;
  // x flux divergence
  work = 0.0;
  merr=MPI_Allreduce((void*)&(Gx_loc(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     Gx_loc.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  work /= subdomain.np;
  Gx_loc = work;
#if NDIM > 1
  // y flux divergence
  work = 0.0;
  merr=MPI_Allreduce((void*)&(Gy_loc(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     Gy_loc.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  work /= subdomain.np;
  Gy_loc = work;
#if NDIM > 2
  // z flux divergence
  work = 0.0;
  merr=MPI_Allreduce((void*)&(Gz_loc(INDICIES(0,0,0))),
		     (void*)&(work(INDICIES(0,0,0))),
		     Gz_loc.size(),
		     MPI_FTYPE,MPI_SUM,
		     subdomain.internal_comm);
  work /= subdomain.np;
  Gz_loc = work;
#endif
#endif

  /* Loop through each species: copy fluids */
  for (int id=0;id<ndist;id++) {
    if (id != electron_dist) {
      if (method[id] < 0) { // Need to add subcycling (cf. charges.cc)
	qden_flux(qden_tmp,INDICIES(Gx_tmp,Gy_tmp,Gz_tmp),
		  id,pic,fspecie,1.0);
	qden_loc += qden_tmp;
	Gx_loc += Gx_tmp;
#if NDIM > 1
	Gy_loc += Gy_tmp;
#if NDIM > 2
	Gz_loc += Gz_tmp;
#endif
#endif
      }
    }
  }

  /* Add background dust layer */
  if (dust_shape != 0) background_dust(qden_loc,(FTYPE)pic[1].n0avg);

  /* Pass pad cells around and sum */
#ifdef USE_DOMAINS
  void pass_sum_guard(FArrayND &, int nx);
  pass_sum_guard(qden_loc,nx);
  pass_sum_guard(Gx_loc,nx);
#if NDIM > 1
  pass_sum_guard(Gy_loc,nx);
#if NDIM > 2
  pass_sum_guard(Gz_loc,nx);
#endif
#endif
#endif

  /* all processors need to wait for root to finish its work */
  MPI_Barrier(MPI_COMM_WORLD);

  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"qden_loc-%06d-",it);
    write_local_bin(qden_loc,wlb_name);
    sprintf(wlb_name,"Gx_loc-%06d-",it);
    write_local_bin(Gx_loc,wlb_name);
    sprintf(wlb_name,"Gy_loc-%06d-",it);
    write_local_bin(Gy_loc,wlb_name);
  }

  /* Write non-ghosted data to ghosted arrays */
  qden = 0.0;
  Gx = 0.0;
#if NDIM > 1
  Gy = 0.0;
#if NDIM > 2
  Gz = 0.0;
#endif
#endif
  for (int ix=0;ix<nx;ix++) {
    for (int iy=0;iy<ny;iy++) {
      for (int iz=0;iz<nz;iz++) {
	qden(INDICIES(ix,iy,iz)) = qden_loc(INDICIES(ix,iy,iz));
	Gx(INDICIES(ix,iy,iz)) = Gx_loc(INDICIES(ix,iy,iz));
#if NDIM > 1
	Gy(INDICIES(ix,iy,iz)) = Gy_loc(INDICIES(ix,iy,iz));
#if NDIM > 2
	Gz(INDICIES(ix,iy,iz)) = Gz_loc(INDICIES(ix,iy,iz));
#endif
#endif
      }
    }
  }

  /* Take care of global BC in ghosted arrays */
  if (boundary_type[0] == PERIODIC) { // Calls pass_guards internally
    periodic_filter(qden,qnx_guard_size);
    periodic_filter(Gx,qnx_guard_size);
#if NDIM > 1
    periodic_filter(Gy,qnx_guard_size);
#if NDIM > 2
    periodic_filter(Gz,qnx_guard_size);
#endif
#endif
  } else if (boundary_type[0] == OPEN) {
    if (subdomain.id_number == 0) {
      for (int id=0;id<ndist;id++) {
  	if (id != electron_dist) {
	  for (int iy=0;iy<ny;iy++) {
	    for (int iz=0;iz<nz;iz++) {
	      qden(INDICIES(-1,iy,iz)) = qden(INDICIES(0,iy,iz));
	      Gx(INDICIES(-1,iy,iz)) = Gx(INDICIES(0,iy,iz));
#if NDIM > 1
	      Gy(INDICIES(-1,iy,iz)) = Gy(INDICIES(0,iy,iz));
#if NDIM > 2
	      Gz(INDICIES(-1,iy,iz)) = Gz(INDICIES(0,iy,iz));
#endif
#endif
	    }
	  }
	}
      }
    }
    if (subdomain.id_number == nsubdomains-1) {
      for (int id=0;id<ndist;id++) {
  	if (id != electron_dist) {
	  for (int iy=0;iy<ny;iy++) {
	    for (int iz=0;iz<nz;iz++) {
	      qden(INDICIES(nx,iy,iz)) = qden(INDICIES(nx-1,iy,iz));
	      Gx(INDICIES(nx,iy,iz)) = Gx(INDICIES(nx-1,iy,iz));
#if NDIM > 1
	      Gy(INDICIES(nx,iy,iz)) = Gy(INDICIES(nx-1,iy,iz));
#if NDIM > 2
	      Gz(INDICIES(nx,iy,iz)) = Gz(INDICIES(nx-1,iy,iz));
#endif
#endif
	    }
	  }
	}
      }
    }
    pass_guards(qden,qnx_guard_size);
    pass_guards(Gx,qnx_guard_size);
#if NDIM > 1
    pass_guards(Gy,qnx_guard_size);
#if NDIM > 2
    pass_guards(Gz,qnx_guard_size);
#endif
#endif
  } else if (boundary_type[0] == INJECT) {
    pass_guards(qden,qnx_guard_size);
    pass_guards(Gx,qnx_guard_size);
#if NDIM > 1
    pass_guards(Gy,qnx_guard_size);
#if NDIM > 2
    pass_guards(Gz,qnx_guard_size);
#endif
#endif
    if (subdomain.id_number == 0) {
      for (int id=0;id<ndist;id++) {
  	if (id != electron_dist) {
  	  FTYPE n0dist_id = 0;
  	  if (method[id]>=0) n0dist_id = pic[id].n0avg;
  	  else n0dist_id = fspecie[id].n0;
  	  for (int iy=0;iy<ny;iy++) {
  	    for (int iz=0;iz<nz;iz++) {
  	      // // qden_loc(INDICIES(-1,iy,iz)) += n0dist_id;
  	      // qden(INDICIES(-1,iy,iz)) += n0dist_id;
	      qden(INDICIES(-1,iy,iz)) = qden(INDICIES(0,iy,iz));
	      Gx(INDICIES(-1,iy,iz)) = Gx(INDICIES(0,iy,iz));
#if NDIM > 1
	      Gy(INDICIES(-1,iy,iz)) = Gy(INDICIES(0,iy,iz));
#if NDIM > 2
	      Gz(INDICIES(-1,iy,iz)) = Gz(INDICIES(0,iy,iz));
#endif
#endif
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
  	      // // qden_loc(INDICIES(nx,iy,iz)) += n0dist_id;
  	      // qden(INDICIES(nx,iy,iz)) += n0dist_id;
	      qden(INDICIES(nx,iy,iz)) = qden(INDICIES(nx-1,iy,iz));
	      Gx(INDICIES(nx,iy,iz)) = Gx(INDICIES(nx-1,iy,iz));
#if NDIM > 1
	      Gy(INDICIES(nx,iy,iz)) = Gy(INDICIES(nx-1,iy,iz));
#if NDIM > 2
	      Gz(INDICIES(nx,iy,iz)) = Gz(INDICIES(nx-1,iy,iz));
#endif
#endif
  	    }
  	  }
  	}
      }
    }
  } else {
    if (mpi_rank == 0) printf("%s:%d Invalid boundary_type (or boundary_type_x) =%d\n",
			      __func__,__LINE__,boundary_type[0]);
  }

  // Let's time the amount of time in this routine:
  collect_time += times(&times_buf)-start_time;

}
#else
void quasineutral_den_flux(FArrayND_ranged &qden,
			   INDICIES(FArrayND_ranged &Gx, \
				    FArrayND_ranged &Gy, \
				    FArrayND_ranged &Gz),
			   particle_dist *pic,fluid *fspecie,
			   int it)
{ /* Dummy routine */ }
#endif // USE_QN
