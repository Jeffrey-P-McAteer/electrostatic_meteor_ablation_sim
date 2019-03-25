// All I/O should be preformed here 
// Note: charge array will be used as workspace and will return modified

using namespace std;
#include <cmath>
#include <cstdio> 
#include <string.h> 

#include <sys/times.h> //Clock related functions
#include <sys/stat.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "eppic.h"
#include "eppic-io.h"
#include "eppic-mpi.h"
#include "eppic-times.h"
#include "eppic-system.h"


#ifdef HAVE_T3PIO
#include "t3pio.h"
#endif

extern hid_t H5_FID;
extern hid_t h5_gid;
FILE* fopensafe(char* filename, char* mode);
FILE* fopensafe(char* filename, char* mode, unsigned long skip);

extern int output_FTarray(FILE *fname, FArrayND &A, eppic_system &sys,
			  FTYPE Afrac, FTYPE kmax);
extern int output_FTarray_h5(hid_t h5_id, char *varname, FArrayND &A, 
			     eppic_system &sys, FTYPE Afrac, FTYPE kmax);
extern int output_FTarray_collective_h5(hid_t h5_id, char *varname, FArrayND &A, 
			     eppic_system &sys, FTYPE Afrac, FTYPE kmax);

void output_vector(FILE *file, PTYPEAVec &x, int imax, PTYPE scaler=1.0, FTYPE shifter=0.0)
{
  
  // Output a vector in single precision binary, scaling by scaler
  // Only output the first n/npout points.
  for (int i=0; i< imax; ++i) {
    float f=x[i]*scaler+shifter;
    fwrite(&f,sizeof(f),1,file);
  }   
}


void output_vector_h5(hid_t h5_id, char *varname, PTYPEAVec &x, int imax, 
                      PTYPE scaler=1.0)
{
  
  // Output a vector in single precision binary, scaling by scaler
  // Only output the first n/npout points.
  float f[imax];
  for (int i=0; i< imax; ++i) {
    f[i] = x[i]*scaler;
  }
  hsize_t dims[1]={(hsize_t)imax};
  H5LTmake_dataset(h5_id, varname, 1, dims, H5T_NATIVE_FLOAT, f);
  
}

// void output(char* indir, particle_dist *pic, fluid *fspecie, field &Efield, 
// 	    FArrayND &rho, FArrayND_ranged &qden, FArrayND_ranged &divG, 
// 	    particle_misc *misc, int it, int it0)
void output(char* indir, particle_dist *pic, fluid *fspecie, field &Efield, 
	    FArrayND &rho, FArrayND_ranged &qden, 
	    INDICIES(FArrayND_ranged &Gx,
		     FArrayND_ranged &Gy,
		     FArrayND_ranged &Gz),
	    particle_misc *misc, int it, int it0)
{
  
  static int first_entry=TRUE;
  static FILE *fphi=0;
  static FILE *fEx=0, *fEy=0, *fEz=0;
  //static FILE *fcharge = 0;
  //static FILE *fdivj=0;
  static FILE *fcons=0;
  static FILE *fvx[MAXDIST]={0}, *fvy[MAXDIST]={0}, *fvz[MAXDIST]={0};
  static FILE *fx[MAXDIST]={0}, *fy[MAXDIST]={0}, *fz[MAXDIST]={0}; 
  static FILE *fden[MAXDIST]={0};
  //int ngrid[] = {INDICIES(nx,ny,nz)}; // Changed by Glenn on 5/2/2018
  int ngrid[] = {INDICIES(nx,ny+yguard_size,nz+zguard_size)};
  int gridG[] = {INDICIES(nx*nsubdomains,ny,nz)};
  char name[256];
  
  static char h5_groupname[256] = "\0";
  static char h5_fphi[256] = "\0";
  static char h5_fphift[256] = "\0";
  static char h5_fExft[256] = "\0",
              h5_fEyft[256] = "\0",
              h5_fEzft[256] = "\0";
  static char h5_fEx[256] = "\0", 
              h5_fEy[256] = "\0", 
              h5_fEz[256] = "\0";
  static char h5_fcons[256] = "\0";
  static char h5_fvx[MAXDIST][256] = {"\0"}, 
              h5_fvy[MAXDIST][256] = {"\0"}, 
              h5_fvz[MAXDIST][256] = {"\0"};
  static char h5_fx[MAXDIST][256] = {"\0"}, 
              h5_fy[MAXDIST][256] = {"\0"}, 
              h5_fz[MAXDIST][256] = {"\0"}; 
  static char h5_fden[MAXDIST][256] = {"\0"};
  static char h5_fdenft[MAXDIST][256] = {"\0"};

  //definition of eppic_system sys
  sys.nx = nx;
  sys.ny = ny;
  sys.nz = nz;
  sys.dx = dx;
  sys.dy = dy;
  sys.dz = dz;
  sys.ndim = NDIM;

  // Create the MPI info object and initialize it.
  MPI_Info info = MPI_INFO_NULL;
  MPI_Info_create(&info);

  // The directory for output files will be stored here
  extern void density(FArrayND &den, int id, particle_dist *pic, 
		      fluid *, FTYPE);
  
  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);

  // Create parallel HDF files - these files must be created outside of first_entry block,
  // since new files must be created every timestep
  if (hdf_output_arrays == 2 && subdomain.rank == 0) {
      // create H5 parallel files
      // If HAVE_T3PIO is defined then use t3pio_set_info to the appropriate 
      // number of stripes for the current number of domains.
	
#ifdef HAVE_T3PIO
      t3pio_set_info(subdomain.neighbor_comm, info, outdir);
#endif
	
      //Give file access to all processors
      hid_t plist_id;
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, subdomain.neighbor_comm, info);
	
      // Create file collectively
      sprintf(name,"%sparallel/parallel%06d.h5", outdir, it);
      H5_FID = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
      if (H5_FID < 0) terminate(-1,"Parallel file creation unsuccessful\n");
      H5Fclose(H5_FID);
      H5Pclose(plist_id);
      MPI_Info_free(&info);
    }

  
  if (first_entry) {
    first_entry=FALSE;
    unsigned long skip2;

    if (subdomain.rank == 0) {
      
      // Open the necessary files 
      char *opentype;
      if (it == 0) opentype="w\0 ";
      else opentype="at\0";
      char *openbtype;
      if (it == 0) openbtype="wb\0 ";
      else openbtype="ab\0";

      unsigned long asize=sizeof(OTYPE)*max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1);
      //      if (charge_out_subcycle > 0)
      //	fcharge=fopensafe(strcat(strcpy(whole_name,outdir),"charge.bin"),openbtype,
      //			  asize*((it-1)/(nout*charge_out_subcycle)+1));

      if ( hdf_output_arrays == 1 ) {
	// create H5 domain files, or open files upon restart
	char temp[128];
	sprintf(temp,"%sdomain%03d.h5",outdir,subdomain.id_number);
	if (it0 == 1) H5_FID = H5Fcreate(temp, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (it0 != 1) H5_FID = H5Fopen(temp, H5F_ACC_RDWR, H5P_DEFAULT);
	if (H5_FID < 0) terminate(-1,"File creation unsucessful\n");
	H5Fclose(H5_FID);
      } // Parallel HDF file creation moved outside of first_entry block - LKT 8/12/16
      
      if (phi_out_subcycle > 0) {
	// for INJECT, phi has extra cells left and right (3 total)
	unsigned long phiasize=asize;
        // 3 lines below commented out by Glenn Sugar 20180412.  We don't need ghost cell output
	//if (boundary_type[0] == INJECT) {
        // phiasize = asize/max(nx/nout_avg,1)*max((nx+3)/nout_avg,1);
	//}
        phiasize = sizeof(OTYPE)*max(nx/nout_avg,1)*max((ny+yguard_size)/nout_avg,1)*max((nz+zguard_size)/nout_avg,1);
	if (hdf_output_arrays == 0) {
	  sprintf(name,"%sdomain%03d/phi.bin",outdir,subdomain.id_number);
	  fphi=fopensafe(name,
			 openbtype, phiasize*((it-1)/(nout*phi_out_subcycle)+1));
        }
	/* Make directories for phift */
	char temp[128];
	if (hdf_output_arrays == 0 && ft_output_arrays == 1) {
	  int len=sprintf(temp,"%sdomain%03d/phift",outdir,subdomain.id_number);
	  int mkerr=mkdir(temp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	strcpy(h5_fphi, "phi");
	strcpy(h5_fphift, "phift");
      }
      //	for (int id=0; id<ndist; ++id) {
      //	  if (phistore_out_subcycle[id] > 0 && subcycle[id] > 1) {
      //	    sprintf(name,"%sphistore%d.bin",outdir,id);
      //	    fphistore[id]=fopensafe(name,openbtype,
      //				    asize*((it-1)/(nout*phistore_out_subcycle[id])+1));
      //	  } }

      if (E_out_subcycle > 0) {
	if (hdf_output_arrays == 0) {
	  sprintf(name,"%sdomain%03d/Ex.bin",outdir,subdomain.id_number);
	  fEx=fopensafe(name,
			openbtype, asize*((it-1)/(nout*E_out_subcycle)+1));
	}
        strcpy(h5_fEx, "Ex");
        strcpy(h5_fExft, "Exft");
	
	if (ndim >= 2) {
	  if (hdf_output_arrays == 0) {
	    sprintf(name,"%sdomain%03d/Ey.bin",outdir,subdomain.id_number);
	    fEy=fopensafe(name,
			  openbtype, asize*((it-1)/(nout*E_out_subcycle)+1));
	  }
	  strcpy(h5_fEy, "Ey");
	  strcpy(h5_fEyft, "Eyft");
	}
	
	if (ndim == 3) {
	  if (hdf_output_arrays == 0) {
	    sprintf(name,"%sdomain%03d/Ez.bin",outdir,subdomain.id_number);
	    fEz=fopensafe(name,
			  openbtype, asize*((it-1)/(nout*E_out_subcycle)+1));
	  }
          strcpy(h5_fEz, "Ez");
	  strcpy(h5_fEzft, "Ezft");
	}
      }
      
      if (divj_out_subcycle > 0) {
	sprintf(name,"%sdomain%03d/divj.bin",outdir,subdomain.id_number);
	// The output routine that uses this file is commented out below:
	//if (hdf_output_arrays == 0) fdivj=fopensafe(name,openbtype, asize*((it-1)/(nout*divj_out_subcycle)+1));
      }
      
      for (int id2=0; id2<ndist; id2++) {
	if (method[id2] != -4) {
	  sprintf(name,"%sdomain%03d/den%d.bin",outdir,subdomain.id_number,id2);
	  if (den_out_subcycle[id2] > 0) {
	    if (hdf_output_arrays == 0) 
	      fden[id2]=fopensafe(name, openbtype, 
				  asize*((it-1)/(nout*den_out_subcycle[id2])+1));
	    sprintf(h5_fden[id2], "den%d", id2);
	    sprintf(h5_fdenft[id2], "denft%d", id2);
	    if (ft_output_arrays == 1) {
	      char temp[128];
	      if (hdf_output_arrays == 0) {
		int len=sprintf(temp,"%sdomain%03d/denft%d",outdir,subdomain.id_number,id2);
		int mkerr=mkdir(temp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	      }
	    }
	  }
	} // if (method[id2] != -4)
      } // for (int id2=0; id2<ndist; id2++)

      sprintf(name,"%sdomain%03d/conserved.out",outdir,subdomain.id_number);
      fcons  = fopensafe(name,opentype,(it-1)/nout+1);
      if (mpi_rank==0) {
	printf("; Time Step: %10s %10s %10s","W particle","W field","Nparts\n");
	fprintf(fcons,"; Time Step: %12s %12s","W particle","W field\n");
      }

      // Open files for particle data
      for (int id=0; id<ndist; ++id) {
	if (method[id] >= 0) {
	  skip2=((static_cast<unsigned long>(it)-1)/(nout*part_out_subcycle[id])+1)*
	    pic[id].x.length()/npout*sizeof(float);

	  if (hdf_output_arrays == 0) {
	    sprintf(name,"%sdomain%03d/x%d.bin",outdir,subdomain.id_number,id);
	    fx[id]  = fopensafe(name,openbtype, skip2);
	    sprintf(name,"%sdomain%03d/vx%d.bin",outdir,subdomain.id_number,id);
	    fvx[id]  = fopensafe(name, openbtype, skip2);
	  }

	  sprintf(h5_fvx[id], "vx%d", id);
	  sprintf(h5_fx[id], "x%d", id);

	  if (ndim >= 2) {
	    if (hdf_output_arrays == 0) {
	      sprintf(name,"%sdomain%03d/vy%d.bin",outdir,subdomain.id_number,id);
	      fvy[id]  = fopensafe(name, openbtype, skip2);
	      sprintf(name,"%sdomain%03d/y%d.bin",outdir,subdomain.id_number,id);
	      fy[id]  = fopensafe(name,openbtype, skip2);
	    }
	    sprintf(h5_fy[id], "y%d", id);
	    sprintf(h5_fvy[id], "vy%d", id);
	  }
	  
	  if (vel_dim[id]==3) {
	    if (hdf_output_arrays == 0) {
	      sprintf(name,"%sdomain%03d/vz%d.bin",outdir,subdomain.id_number,id);
	      fvz[id]  = fopensafe(name, openbtype, skip2);
	      sprintf(name,"%sdomain%03d/z%d.bin",outdir,subdomain.id_number,id);
	      fz[id]  = fopensafe(name,openbtype, skip2);
	    }
	    sprintf(h5_fz[id], "z%d", id);
	    sprintf(h5_fvz[id], "vz%d", id);
	  }
	}
      }
      
      void output_idl_file(char *);
      sprintf(name,"%sparams.pro",outdir);
      output_idl_file(name);
      
    } //end if (subdomain_rank == 0) 

    // Flux output subcycling should be >= den_out_subcycling)
    for (int id2=0; id2<ndist; id2++) 
      if (flux_out_subcycle[id2] < den_out_subcycle[id2]) 
	flux_out_subcycle[id2] = den_out_subcycle[id2];

    // Nvsqr output subcycling should be >= den_out_subcycling)
    for (int id3=0; id3<ndist; id3++) 
      if (nvsqr_out_subcycle[id3] > 0 && nvsqr_out_subcycle[id3] < den_out_subcycle[id3]) 
	nvsqr_out_subcycle[id3] = den_out_subcycle[id3];

  } // end if (first_entry) 

  if (mpi_rank == 0) printf("; %9d:",it);

  if (subdomain.rank == 0) {

    if (hdf_output_arrays > 0) {
      // give HDF files necessary properties
      if (hdf_output_arrays == 1) {
	char temp[128];
	sprintf(temp,"%sdomain%03d.h5",outdir,subdomain.id_number);
	H5_FID = H5Fopen(temp, H5F_ACC_RDWR, H5P_DEFAULT);
	sprintf(h5_groupname, "/time_%06d\0", it);
	h5_gid = H5Gcreate(H5_FID, h5_groupname, H5P_DEFAULT, H5P_DEFAULT, 
			   H5P_DEFAULT);
      } else if (hdf_output_arrays == 2) {
	hid_t plist_id;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, subdomain.neighbor_comm, info);
	sprintf(name,"%sparallel/parallel%06d.h5", outdir, it);
	H5_FID = H5Fopen(name, H5F_ACC_RDWR, plist_id);
	if (H5_FID < 0) terminate(-1,"Parallel file opened unsuccessful\n");
	H5Pclose(plist_id);
      }
    }
    
    if (Efield.Ex.size() == Efield.phi_rho.size() && it%(nout*E_out_subcycle)==0) {
      if (ft_output_arrays == 0 || it%(full_array_nout*E_out_subcycle)==0) {
	if ( hdf_output_arrays == 0 ) 
	  output_array(fEx, Efield.Ex,ngrid,nout_avg);
	else if ( hdf_output_arrays == 1 ) 
	  output_array_h5(h5_gid, h5_fEx, Efield.Ex,ngrid,nout_avg);
	else if ( hdf_output_arrays == 2 ) output_collective_array_h5(H5_FID, 
								      h5_fEx, Efield.Ex, ngrid, nout_avg);      
      } if (ft_output_arrays == 1) {
	if (hdf_output_arrays == 0) {
	  sprintf(name,"%sdomain%03d/Exft/d%09d.bin",outdir,subdomain.id_number,it);
	  FILE *fExFT=fopensafe(name, "wb\0", 0);
	  output_FTarray(fExFT, Efield.Ex, sys, Eft_out_min, Eft_out_kmax);
	  fclose(fExFT);
	}
	else if (hdf_output_arrays == 1) output_FTarray_h5(h5_gid, h5_fExft, Efield.Ex, 
							   sys, Eft_out_min, Eft_out_kmax);
	else if (hdf_output_arrays == 2) output_FTarray_collective_h5(H5_FID, h5_fExft, Efield.Ex, 
								      sys, Eft_out_min, Eft_out_kmax);
      }
    }

    int nghost[2]={1,2};
    // 4 lines below commented out by Glenn 20180412
    //if (boundary_type[0] == INJECT) {
    //  nghost[0]=0;
    //  nghost[1]=0;
    //}

    if (it%(nout*phi_out_subcycle)==0) {
      if (ft_output_arrays == 0 || it%(full_array_nout*phi_out_subcycle)==0) { //output full arrays
	if ( hdf_output_arrays == 0 ) {
          // Added by Glenn for debugging memory bug
          //printf("Outputting Efield.phi_rho\n");
          Efield.phi_rho.bin_output_ghost(fphi,nghost,nout_avg);
        }
	else if ( hdf_output_arrays == 1 ) 
	  Efield.phi_rho.bin_output_ghost_h5(h5_gid, h5_fphi,nghost,nout_avg);
	else if ( hdf_output_arrays == 2 ) 
	  Efield.phi_rho.bin_output_ghost_collective_h5(H5_FID, h5_fphi,nsubdomains, subdomain.id_number,
							nghost ,nout_avg);
      } if (ft_output_arrays == 1) { //output fourier transformed arrays
	if (hdf_output_arrays == 0) {
	  sprintf(name,"%sdomain%03d/phift/d%09d.bin",outdir,subdomain.id_number,it);
	  FILE *fphiFT=fopensafe(name, "wb\0", 0);
	  output_FTarray(fphiFT, Efield.phi_rho, sys, phift_out_min, phift_out_kmax);
	  fclose(fphiFT);
	}
	else if (hdf_output_arrays == 1) output_FTarray_h5(h5_gid, h5_fphift, Efield.phi_rho, 
							   sys, phift_out_min, phift_out_kmax);
	else if (hdf_output_arrays == 2) output_FTarray_collective_h5(H5_FID, h5_fphift, Efield.phi_rho, 
								      sys, phift_out_min, phift_out_kmax);
      }    
    }
        
    // Output vector quantities such as particles:
    for (int id=0; id<ndist; ++id) {
      if (method[id]>=0) 
	if (it%(nout*part_out_subcycle[id]) == 0) {
	  if ( hdf_output_arrays == 0 )
	    output_vector(fvx[id], pic[id].vx, pic[id].vx.length()/npout, dx/dt);
	  else if ( hdf_output_arrays == 1 ) 
            output_vector_h5(h5_gid, h5_fvx[id], pic[id].vx, 
                             pic[id].vx.length()/npout, dx/dt);
	  if (ndim >= 2) {
	    if ( hdf_output_arrays == 0 ) 
	      output_vector(fvy[id], pic[id].vy, pic[id].vy.length()/npout, dy/dt);
	    else  if ( hdf_output_arrays == 1 ) 
              output_vector_h5(h5_gid, h5_fvy[id], pic[id].vy, 
                               pic[id].vy.length()/npout, dy/dt);
	  }
	  if (vel_dim[id]==3) {
	    if ( hdf_output_arrays == 0 ) 
	      output_vector(fvz[id], pic[id].vz, pic[id].vz.length()/npout, dz/dt);
	    else  if ( hdf_output_arrays == 1 ) 
              output_vector_h5(h5_gid, h5_fvz[id], pic[id].vz, 
                               pic[id].vz.length()/npout, dz/dt);
	  }	  
	  if ( hdf_output_arrays == 0 )
	    output_vector(fx[id], pic[id].x, pic[id].x.length()/npout, dx, nx*subdomain.id_number*dx);
	  else  if ( hdf_output_arrays == 1 ) 
            output_vector_h5(h5_gid, h5_fx[id], pic[id].x, 
                             pic[id].x.length()/npout, dx);
	  if (ndim >= 2) {
	    if ( hdf_output_arrays == 0 ) 
	      output_vector(fy[id], pic[id].y, pic[id].y.length()/npout, dy);
	    else  if ( hdf_output_arrays == 1 ) 
              output_vector_h5(h5_gid, h5_fy[id], pic[id].y, 
                               pic[id].y.length()/npout, dy);
	  }
	  if (ndim==3) {
	    if ( hdf_output_arrays == 0 ) 
	      output_vector(fz[id], pic[id].z, pic[id].z.length()/npout, dz);
	    else  if ( hdf_output_arrays == 1 ) 
              output_vector_h5(h5_gid, h5_fz[id], pic[id].z, 
                               pic[id].z.length()/npout, dz);
	  }
	}
    }
    if ( hdf_output_arrays == 1 ) H5Gclose(h5_gid);
  }// End if (subdomain.rank == 0)

  // Calculate and output div(J) 
  if (divj_out_subcycle>0 && it%(nout*divj_out_subcycle) == 0) {
    // Use rho as array
    void divj(FArrayND &rho, particle_dist *pic, fluid *fspecie);
    // divj is not working
    // a) it's still periodic
    // b) it's broken for (domain) fluids 
    //    divj(rho, pic, fspecie);
    
    //    if (subdomain.rank == 0) output_array(fdivj, rho);
  } // End of div(J) section

  // Output perturbed densities 
  for (int id=0; id<ndist; ++id) {
    if ( (den_out_subcycle[id] >0 && it%(nout*den_out_subcycle[id]) == 0) ) {
      int tmp_size[]={INDICIES(nx+xguard_size,ny+yguard_size,nz+zguard_size)};
      FArrayND den_tmp=FArrayND(tmp_size);
      // if added by Glenn 5/2/2018
      if (unscale_density){
        // We want the absolute density
        density(den_tmp, id, pic, fspecie, 1.);
      } else {
        // We want the perturbed density
        density(den_tmp, id, pic, fspecie, 1./pic[id].n0avg);
      }
#ifdef USE_MPI
      // Pass the arrays around the networked processes 
      if (method[id] >= 0) {
	FArrayND den_all=den_tmp;
	den_all=0;
	
	int mpi_err=MPI_Reduce((void*) &(den_tmp(INDICIES(0,0,0))),
			       (void*) &(den_all(INDICIES(0,0,0))),
			       den_all.size(), MPI_FTYPE, MPI_SUM, 0, 
			       subdomain.internal_comm);
	if ( mpi_err != MPI_SUCCESS) 
	  mpi_error(mpi_rank, mpi_err,"Reduce call in output");
	
        // den_all /= subdomain.np was commented.  Comment removed by Glenn 4/19/2018
	den_all /= subdomain.np;
	den_tmp=den_all;
      }
#endif
      if (subdomain.rank == 0) {
        if ( hdf_output_arrays == 1 ) {
	  h5_gid = H5Gopen(H5_FID, h5_groupname, H5P_DEFAULT);
	}
#ifdef USE_DOMAINS
	void pass_sum_guard(FArrayND &, int);
	pass_sum_guard(den_tmp,nx);
#endif

        if (!unscale_density) den_tmp -= 1.0;

	//Produce an artificial density for testing purposes
	// int nmodes=3;

	/*
	if (mpi_rank == 0) 
	  printf("!!!!\n\t\tWARNING: %s using artificial density for testing \n!!!!\n",__func__);
	// int nmodes=1;
	// double kx[]={2*2.*M_PI/(nx*nsubdomains),20*2.*M_PI/(nx*nsubdomains),4*2.*M_PI/(nx*nsubdomains)};
	// double ky[]={2*(2.*M_PI/ny),8*(2.*M_PI/ny),4*(2.*M_PI/ny)};
	// double kz[]={2*(2.*M_PI/nz),8*(2.*M_PI/nz),4*(2.*M_PI/nz)};
	// for (int ix = 0; ix < nx; ++ix)
	//   for (int iy = 0; iy < ny; ++iy)
	//     for (int iz = 0; iz < nz; ++iz) {
	//       den_tmp(INDICIES(ix, iy, iz))= 0;
	//       for (int im=0; im<nmodes; im++)
	// 	den_tmp(INDICIES(ix, iy, iz)) += sin(kx[im]*(ix+subdomain.id_number*nx))
	// 	  *cos(ky[im]*iy)*cos(kz[im]*iz);
	//     }
	FTYPE Lx=nx*nsubdomains-1,Ly=ny-1,Lz=nz-1;
	FTYPE Px=2,Py=2,Pz=2;
	FTYPE Hx=Px/Lx,Hy=Py/Ly,Hz=Pz/Lz;
	for (int ix = 0; ix < nx; ++ix)
	  for (int iy = 0; iy < ny; ++iy)
	    for (int iz = 0; iz < nz; ++iz)
	      den_tmp(INDICIES(ix,iy,iz)) = cos(Hx*M_PI*ix)*cos(Hy*M_PI*iy)*cos(Hz*M_PI*iz);
	*/

	if (it%(nout*den_out_subcycle[id]) == 0) {
	  if (ft_output_arrays == 0 || it%(full_array_nout*den_out_subcycle[id])==0) {
	    if ( hdf_output_arrays == 0 ) 
	      output_array(fden[id], den_tmp, ngrid, nout_avg);
	    else if ( hdf_output_arrays == 1 ) 
	      output_array_h5(h5_gid, h5_fden[id], den_tmp, ngrid, nout_avg);
	    else if ( hdf_output_arrays == 2 ) 
	      output_collective_array_h5(H5_FID, h5_fden[id], den_tmp, ngrid, nout_avg);
	  } 
	  if (ft_output_arrays == 1) {
	    if (hdf_output_arrays == 0) {
	      sprintf(name,"%sdomain%03d/denft%d/d%09d.bin",outdir,subdomain.id_number,id,it);
	      FILE *fdenFT=fopensafe(name, "wb\0", 0);
	      output_FTarray(fdenFT, den_tmp, sys, denft_out_min[id], denft_out_kmax[id]);
	      fclose(fdenFT);
	    }
	    else if (hdf_output_arrays == 1) output_FTarray_h5(h5_gid, h5_fdenft[id], den_tmp, 
			       sys, denft_out_min[id], denft_out_kmax[id]);
	    else if (hdf_output_arrays == 2) output_FTarray_collective_h5(H5_FID, h5_fdenft[id], 
			  den_tmp, sys, denft_out_min[id], denft_out_kmax[id]);
	  }
	}
	if ( hdf_output_arrays == 1 ) H5Gclose(h5_gid);
      }
    }
  } // End of density output section

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif	  

  // Fix rho (non-quasineutral) or qden & divG (quasineutral)
  // in case later routines need it
  if (efield_algorithm == 2)
    // quasineutral_den_flux(qden, divG, pic, fspecie, it);
    quasineutral_den_flux(qden, INDICIES(Gx,Gy,Gz), pic, fspecie, it);
  else
    charges(rho, pic, fspecie, it);

  // Output velocity distribution, ignoring spacial localization
  void output_vdist(particle_dist *pic, int it, char *outdir, int it0);
  sprintf(name,"%sdomain%03d/",outdir,subdomain.id_number);
  output_vdist(pic, it, name, it0);

  // Output full distribution functions
  extern void output_fdist(particle_dist *pic, int, char *);
  output_fdist(pic, it, name);

  // Output flux arrays
  extern void output_fluxes(particle_dist *pic, fluid *, int, char *, int it0);
  output_fluxes(pic, fspecie, it, name, it0);

  // Output nvsqr arrays
  extern void output_nvsqr(particle_dist *pic, int, char *, int it0);
  output_nvsqr(pic, it, name, it0);

  // Conserved properties:
  void  output_momentum(char*, particle_dist *, fluid *, int);
  output_momentum(name, pic, fspecie, it);
  
  void output_energy(FILE* fcons, particle_dist *, field &);
  output_energy(fcons, pic, Efield);

  // Output particle statistics
  void output_particle_stats(int, char *, particle_dist *);
  if (method.max() >= 0) output_particle_stats(it, name, pic);
  // Broadcast l_debye to all processors. Probably breaks if FTYPE = float
  MPI_Bcast(&l_debye, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Diagnostic check for coulomb collisions to output l_debye into eppic.o:
  // if (mpi_rank == 0) printf("; d0= %6g", l_debye);
 

  // Add a CR to the output file
  if (mpi_rank == 0) printf("\n");

  fflush((FILE *) 0);

  // DUMP to enable a restart capability 
  static int dcount=0;
  dcount += nout;
  FTYPE time_lapsed = difftime(time(NULL),sim_start_time)/60;
  // every thing is counted in minutes, explaining '/60' above
  // noutput_before_limit explained:
  //  noutput_before_limit = (time_limit-time_lapsed)/time_per_output;
  //  time_per_output is undefined for it=0, so default to 1
  FTYPE noutput_before_limit = 1;
  if (it>0) noutput_before_limit = (time_limit-time_lapsed)/(time_lapsed/it*nout);
  void dump(particle_dist *pic, fluid *fspecie, field &Efield, int it);
  if (it != 0 && iwrite !=0 && (dcount > iwrite || it == nt || 
       (noutput_before_limit < 1 && time_limit>0))) {

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif	  
    
    // This loop insures that only nfiles_parallel processors 
    // outputs at a time. Something neccesary for the BU machines.
    int ndumps=static_cast<int>(ceil(FTYPE(mpi_np/nfiles_parallel)));
    if (ndumps == 0) ndumps = 1;
    for (int ipro=0; ipro<ndumps; ++ipro) {
      if (ipro == mpi_rank%ndumps) dump(pic,fspecie,Efield,it);
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif	  
    }
    dcount=0;
    if (noutput_before_limit < 1 && time_limit>0) {
      nt = it;
      if (mpi_rank==0) 
	cout << endl
	     << "Ending program because of time limit ("
	     << time_limit
	     << " min.), with time lapsed = "
	     << time_lapsed
	     << " min."
	     << endl;
    }
  }
  fflush((FILE *) 0);
  if (subdomain.rank == 0 && hdf_output_arrays > 0 ) H5Fclose(H5_FID);


  /* all processors need to wait for first to finish it's work */
  MPI_Barrier(MPI_COMM_WORLD);
  // Time information:
  // Let's time the amount of time in this routine:
  output_time += times(&times_buf)-start_time;
  
  void output_times(int it, const char *outdir);
  if (mpi_rank == 0) output_times(it, name);


} //end output

void output_times(int it, const char *dir)
{
  
  static FILE* ftime;
  static int first_entry=TRUE;
  if (first_entry) {
    first_entry=FALSE;
    
    // Open the files 
    char name[256];
    char *opentype;
    if (it == 0) opentype="w\0 ";
    else opentype="at\0";
    sprintf(name,"%stimers.dat",dir);
    ftime=fopensafe(name, opentype);
    fprintf(ftime, "%8s %11s %11s %11s %11s %11s %11s %11s %11s %11s\n","iter #",
	    "Wall Clock","Sys Clock","vadv time", "xadv time",
	    "charge","collect", "efield","output","fluid");
  }

  static clock_t prev_time=times(&times_buf);
  static tms prev_times_buf=times_buf;
  clock_t from_time=times(&times_buf);
    
  fprintf(ftime, "%8d %11ld %11ld %11ld %11lu %11lu %11lu %11lu %11lu  %11lu\n", it,
	  from_time-prev_time, times_buf.tms_utime-prev_times_buf.tms_utime,
	  vadvance_time, xadvance_time, charge_time, collect_time, efield_time,
	  output_time, fluid_time);
  // Reset all the times:
  vadvance_time=0;       // Time spent scattering and advancing velocities
  xadvance_time=0;       // Time spent advancing positions
  charge_time=0;         // Time spent gathering densities
  collect_time=0;        // Time spent collecting densities across processors
  efield_time=0;         // Time spent calculating electric field
  output_time=0;         // Time spent in output routine
  fluid_time=0;          // Time spent in fluid routine

  prev_time=from_time;
  prev_times_buf.tms_utime =  times_buf.tms_utime;
  
}

void output_idl_file(char *name)
{
  
  // We need to store some parameters where IDL can read them 
  FILE *fidl;
  fidl = fopensafe(name,"w\0");
  // Create a file with the velocity distributions readable by IDL 
  for (int id=0; id < ndist; ++id) {
    fprintf(fidl,";Distribution %d velocities\n",id);
    fprintf(fidl,"pnvx%d = %d \npnvy%d = %d\n",id,pnvx[id],id,pnvy[id]);
    fprintf(fidl,"pvxmin%d = %g\npvxmax%d = %g\n",id,pvxmin[id],id,pvxmax[id]);
    fprintf(fidl,"pvymin%d = %g\npvymax%d = %g\n",id,pvymin[id],id,pvymax[id]);
  }
  fclose(fidl);
  
}

void output_particle_stats(int it, char *dir, particle_dist *pic)
{

  // Let's print out some statistical tests of the particle velocity data 
  const int nmom=4;
  int moment(PTYPEAVec &x, PTYPEAVec &moments, PTYPEAVec &flag_array, int np=0);

  static FILE* pstats[MAXDIST];
  static int first_entry=TRUE;
  if (first_entry){
    first_entry=FALSE;

    if (mpi_rank == 0) {
      char name[256];
      char *opentype;
      if (it == 0) opentype="w\0 ";
      else opentype="at\0";
      for (int id=0; id<ndist; ++id) {
	if (method[id]>=0) {
	  sprintf(name,"%smoments%d.out", dir, id);
	  pstats[id]  = fopensafe(name, opentype, (it-1)/nout+2);
	  
	  //Label Line
	  if (iread==0) {
	    const char* mom_names[]= {"Mean","Variance","Moment 3","Moment 4","Moment 5","Moment ..."};
	    const char* label[]={"Vx","Vy","Vz"};
	    fprintf(pstats[id], "it    ");
	    for (int ivel=0; ivel<vel_dim[id]; ivel++) {
	      fprintf(pstats[id], "  %2s:%5s",label[ivel],mom_names[0]);
	      for (int imom=1;imom<nmom;imom++) fprintf(pstats[id], " %10s",mom_names[imom]);
	    }
	    fprintf(pstats[id], "\n");
	  }
	} // end not fluid if
      }
    }
  } // first entry

  PTYPEAVec moments(nmom);
  FTYPE v2_therm = 0;  // Thermal velocity, squared
  //  FTYPE v_bulk = 0; //RM
  l_debye = 0;  // reset global variable

  for (int id=0; id<ndist; ++id) {
    if (method[id]>=0) {
      if (mpi_rank == 0) fprintf(pstats[id],"%6d",it);
      // Vx stats:
      moment(pic[id].vx, moments, pic[id].x, pic[id].np);
      //  v_bulk = moments(0)*dx/dt;  //RM
      v2_therm += moments(1)*dx*dx/dt/dt;
      if (mpi_rank == 0) {
	PTYPE norm = 1;
	for (int imom=0;imom<nmom;imom++) {
	  norm *= dx/dt;
	  moments(imom) *= norm;
	  fprintf(pstats[id]," %10.4g",moments(imom));
	}
      }

    // Vy stats:
    if (pic[id].vy.size() > 1) {
	moment(pic[id].vy, moments, pic[id].x, pic[id].np);
	v2_therm += moments(1)*dy*dy/dt/dt;
      if (mpi_rank == 0) {
	PTYPE norm = 1;
	for (int imom=0;imom<nmom;imom++) {
	  norm *= dy/dt;
	  moments(imom) *= norm;
	  fprintf(pstats[id]," %10.4g",moments(imom));
	}
      }
    }

    // Vz stats:
    if (pic[id].vz.size() > 1) {
      moment(pic[id].vz, moments, pic[id].x, pic[id].np);
      v2_therm += moments(1)*dz*dz/dt/dt;
      if (mpi_rank == 0) {
	PTYPE norm = 1;
	for (int imom=0;imom<nmom;imom++) {
	  norm *= dz/dt;
	  moments(imom) *= norm;
	  fprintf(pstats[id]," %10.4g",moments(imom));
	}
      }
    }
    else{v2_therm *= 1.5;}
    
    if (mpi_rank == 0) fprintf(pstats[id],"\n");

    //    vth_gb[id] = pow(v2_therm, 0.5);  // Record the thermal speed as a global variable
    if (id == 0) { // FIX ME, assumes 0 distribution is always electrons
      l_debye = (3*pic[id].n0avg*qd[id]*qd[id])/(eps*md[id]*v2_therm);
    }
    
    } // end not fluid if
  }

  // Finish calculating Debye length
  l_debye = pow(l_debye, -0.5);

}

void output_momentum(char *dir, particle_dist *pic, fluid *fspecie, int it)
{

  // Calculate and output momentum
  static FILE *fmoment[MAXDIST];
  static int first_entry=TRUE;

  if (mpi_rank == 0)
    if (first_entry) {
      first_entry = FALSE;

      char name[256];
      char *opentype;
      if (it == 0) opentype="w\0 ";
      else opentype="at\0";
      for (int id=0; id<ndist; ++id) {
	if (method[id] != -4) {
	  sprintf(name,"%smisc%d.out", dir, id);
	  fmoment[id]  = fopensafe(name, opentype, (it-1)/nout+1);
	}
      }
    }

  FTYPE volume=nx*dx*ny*dy*nz*dz;
  FTYPEAVec pxpart(ndist), pypart(ndist), pzpart(ndist);

  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0) { // Particle methods:
      for (int i=0;i<pic[id].np;i++) {
	pxpart[id] += pic[id].vx(i)*(dx/dt)*pic[id].m*(pic[id].n0avg/pic[id].np*volume);
	if (vel_dim[id] >= 2) pypart[id] += pic[id].vy(i)*(dy/dt)*
				pic[id].m*(pic[id].n0avg/pic[id].np*volume);
	if (vel_dim[id] >= 3) pzpart[id] += pic[id].vz(i)*(dz/dt)*pic[id].m*
				(pic[id].n0avg/pic[id].np*volume);
      }
    // } else { // Fluid method:
    } else if (method[id] < 0 && method[id] != -4) {

#ifdef USE_DOMAINS // don't sum over guard cells
      extern FTYPE sum_noguard(FArrayND_ranged &in);
      pxpart[id]=sum_noguard(fspecie[id].vx)*fspecie[id].m*dx*dy*dz;
      if (vel_dim[id] >= 2) 
	pypart[id]=sum_noguard(fspecie[id].vy)*fspecie[id].m*dx*dy*dz;
      if (vel_dim[id] == 3)
	pzpart[id]=sum_noguard(fspecie[id].vz)*fspecie[id].m*dx*dy*dz;
#else
      pxpart[id]=fspecie[id].vx.sum()*fspecie[id].m*dx*dy*dz;
      if (vel_dim[id] >= 2) pypart[id]=fspecie[id].vy.sum()*fspecie[id].m*dx*dy*dz;
      if (vel_dim[id] >= 3) pzpart[id]=fspecie[id].vz.sum()*fspecie[id].m*dx*dy*dz;
#endif
    }
  }
  
#ifdef USE_MPI
  {
    FTYPEAVec pall(ndist);
    int mpi_err=MPI_Reduce(pxpart.address(), pall.address(),ndist,
			   MPI_FTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce p call in output_momentum");
    
    pxpart=pall;
    
    mpi_err=MPI_Reduce((void*) &(pypart[0]),(void*) &(pall[0]),ndist,
		       MPI_FTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce p call in output_momentum");
    pypart=pall;
    
    mpi_err=MPI_Reduce((void*) &(pzpart[0]),(void*) &(pall[0]),ndist,
		       MPI_FTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce p call in output_momentum");
    pzpart=pall;
  } 
#endif    

  if (mpi_rank == 0) {
    FTYPE ptotal=pxpart.sum()+pypart.sum()+pzpart.sum();
    //printf(" p=%10.5g", ptotal);
    //    fprintf(fcons," %10.5g", ptotal);

    // Print the momentums in fmisc.
    for (int id=0; id<ndist; ++id) {
      if (method[id] != -4) {
	fprintf(fmoment[id]," %14.9g %14.9g %14.9g %14.9g",pxpart[id],
		pypart[id],pzpart[id],pxpart[id]+pypart[id]+pzpart[id]);
      }
    }
  }
  
}
