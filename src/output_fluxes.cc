// Output of fluxes

#include <stdio.h> 
#include <string.h> 
#include <hdf5.h>
#include <hdf5_hl.h>

#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-io.h"

extern hid_t H5_FID;
extern void output_collective_array_gen_h5(hid_t file_id, char *varname, ArrayNd<OTYPE,NDIM> &a, int ngrid[NDIM], int nout_avg);

void output_fluxes(particle_dist *pic, fluid *fspecie, int it, char *dir, int it0)
{
  
  extern void gather_weight(ArrayNd<OTYPE, NDIM> &den, particle_dist &pic, 
			    PTYPEAVec &b, FTYPE *xmin, 
			    FTYPE *xmax, int *nmesh, int *wrap, FTYPE scale);
  /*
  extern void gather_weight_inject(ArrayNd<OTYPE, NDIM> &den, particle_dist &pic, 
			    PTYPEAVec &b, FTYPE *xmin, 
				   FTYPE *xmax, int *nmesh, int *wrap, FTYPE scale, FTYPE nxmin);
  */

  //  extern void output_array(FILE* fname, FArrayND &a);

  void pass_sum_guard(ArrayNd<OTYPE,NDIM> &, int);

  static FILE *ffluxx[MAXDIST], *ffluxy[MAXDIST], *ffluxz[MAXDIST];
  
  extern hid_t h5_gid;
  char h5_fname[256], h5_groupname[256];
  static char h5_ffluxx[MAXDIST][256] = {"\0"}, 
              h5_ffluxy[MAXDIST][256] = {"\0"}, 
              h5_ffluxz[MAXDIST][256] = {"\0"};
              
  if (subdomain.rank == 0 && hdf_output_arrays == 1 ) {
    sprintf(h5_groupname, "/time_%06d\0", it);
    h5_gid = H5Gopen(H5_FID, h5_groupname, H5P_DEFAULT);
  }

  static int first_entry=TRUE;
  if (first_entry) {
    first_entry=FALSE;

    // Open the files
    if (subdomain.rank == 0) {
    
      for (int id2=0; id2<ndist; id2++) {
	if (method[id2] != -4) {
	  char name[256];
	  char *openbtype;
	  FILE* fopensafe(char* filename, char* mode, unsigned long skip);
	  if (it == 0) openbtype="wb";
	  else openbtype="ab";
	  
	  if (flux_out_subcycle[id2] > 0) {
	    sprintf(name,"%sfluxx%d.bin",dir,id2);
	    unsigned long skip=sizeof(OTYPE)*
	      max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
	      ((static_cast<unsigned long>(it)-1)/(nout*flux_out_subcycle[id2])+1);
	    if (hdf_output_arrays == 0) ffluxx[id2]  = fopensafe(name,openbtype,skip);
	    if ( hdf_output_arrays > 0 )sprintf(h5_ffluxx[id2], "fluxx%d", id2);
	  }

	  if (vel_dim[id2] >= 2) {
	    if (flux_out_subcycle[id2] > 0) {
	      sprintf(name,"%sfluxy%d.bin",dir,id2);
	      int skip=sizeof(OTYPE)*
		max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
		((static_cast<unsigned long>(it)-1)/
		 (nout*flux_out_subcycle[id2])+1);
	      if (hdf_output_arrays == 0) ffluxy[id2]  = fopensafe(name,openbtype,skip);
	      if ( hdf_output_arrays > 0 ) sprintf(h5_ffluxy[id2], "fluxy%d", id2);
	    }
	    
	    if (vel_dim[id2] >= 3) {
	      if (flux_out_subcycle[id2] > 0) {
		sprintf(name,"%sfluxz%d.bin",dir,id2);
		int skip=sizeof(OTYPE)*
		  max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
		  ((static_cast<unsigned long>(it)-1)/
		   (nout*flux_out_subcycle[id2])+1);
		if (hdf_output_arrays == 0) ffluxz[id2]  = fopensafe(name,openbtype,skip);
		if ( hdf_output_arrays > 0 ) sprintf(h5_ffluxz[id2], "fluxz%d", id2);
	      }
	    }
	  }
	} // if (method[id2] != -4)
      } // for (int id2=0; id2<ndist; id2++)
    }

  } // end if (first_entry) 

  for (int id=0; id<ndist; ++id) {
    if (it%(nout*flux_out_subcycle[id]) == 0) {

      // Define array limits
      int gn[]={INDICIES(nx/nout_avg,ny/nout_avg,nz/nout_avg)};
      FTYPE gmin[]={INDICIES(0,0,0)};
      FTYPE gmax[]={INDICIES((FTYPE) nx,(FTYPE) ny,(FTYPE) nz)};
      // Generate a set of pointers to vectors of particles
      int wrap[]={INDICIES(nx/nout_avg,ny/nout_avg,nz/nout_avg)};
      
      // These need to be adjusted if there are subdomains:
      if ((nsubdomains > 1) || (boundary_type[0] == INJECT)) {
	gn[0]+=1;
	wrap[0]+=1; // This eliminates the wrap
	gmax[0]+=1*nout_avg; 
      }
      

      // Calculate total output gridsize
      FTYPE gridsize=nx/nout_avg;
      for (int idim=1; idim<ndim; ++idim) gridsize *= (gmax[idim]-gmin[idim])/nout_avg;
      
      // Define the array to sum onto
      ArrayNd<OTYPE,NDIM> flux(gn);
#ifdef USE_MPI
      // Define an mpi work array
      ArrayNd<OTYPE,NDIM> flux_sum(gn);
#endif

      
      // Output Flux_x:
      
      if (method[id] >=0) { // PIC gather:
	
	//FTYPE scale=gridsize*pic[id].n0avg/pic[id].np*dx/dt;
	FTYPE scale=gridsize*pic[id].n0avg/npd[id]*dx/dt;
	flux=0.;
	/*
	FTYPE nxmin = -1/nout_avg;
	if ((boundary_type[0] == INJECT) && (subdomain.id_number == 0)) 
	  gather_weight_inject(flux, pic[id], pic[id].vx, gmin, gmax, gn, wrap, scale,nxmin);
	else
	*/
	  gather_weight(flux, pic[id], pic[id].vx, gmin, gmax, gn, wrap, scale);

#ifdef USE_MPI
	// Sum the flux arrays accross all processors 
	int mpi_err=MPI_Reduce((void*)&flux(0),(void*)&flux_sum(0), flux.length(), 
			       MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	if ( mpi_err != MPI_SUCCESS) 
	  mpi_error(mpi_rank, mpi_err,"Reduce flux call in output_fluxes failed");
	flux = flux_sum;
	flux /= subdomain.np;
#ifdef USE_DOMAINS
	  pass_sum_guard(flux, nx/nout_avg);
#endif
#endif    

	//Print flux_x:
	  if (subdomain.rank == 0) {
	    if (hdf_output_arrays == 0) fwrite(&flux(0), sizeof(OTYPE), static_cast<int>(gridsize), ffluxx[id]);
	    hsize_t dims[1] = {(hsize_t)gridsize};
	    if ( hdf_output_arrays == 1 )
	      H5LTmake_dataset(h5_gid, h5_ffluxx[id], 1, dims, H5T_NATIVE_FLOAT, &flux(0));
	    int ngrid[] = {INDICIES(nx,ny,nz)}; 
	    if (hdf_output_arrays == 2) 
	      output_collective_array_gen_h5(H5_FID, h5_ffluxx[id], flux, ngrid, 1);
	  }

	
      // } else { //Fluid method
      } else if (method[id] < 0 && method[id] != -4) {
	if (subdomain.rank == 0) {
#ifdef USE_DOMAINS
	  void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
			      INDICIES(int startx, int starty, int startz),
			      INDICIES(int endbeforex, 
				       int endbeforey, 
				       int endbeforez));
	  FArrayND no_guards(FArrayND_ranged &in);
	  FArrayND flux_tmp=no_guards(fspecie[id].vx);
	  sub_array_mult(fspecie[id].den,flux_tmp,
			 INDICIES(0,0,0),INDICIES(nx,ny,nz));
	  int ngrid[] = {INDICIES(nx,ny,nz)}; 
	  if ( hdf_output_arrays == 0 ) 
	    output_array(ffluxx[id], flux_tmp, ngrid, nout_avg);
	  else if ( hdf_output_arrays == 1 ) 
            output_array_h5(h5_gid, h5_ffluxx[id], flux_tmp, ngrid, nout_avg);
	  else if ( hdf_output_arrays == 2 ) 
            output_collective_array_h5(H5_FID, h5_ffluxx[id], flux_tmp, ngrid, nout_avg);

#else
	  FArrayND flux_tmp=fspecie[id].vx;
	  flux_tmp *= fspecie[id].den;
	  if ( hdf_output_arrays == 0 )
	    output_array(ffluxx[id], flux_tmp);
	  else  if ( hdf_output_arrays == 1 ) 
            output_array_h5(h5_gid, h5_ffluxx[id], flux_tmp);
#endif
	}
      }
      // End of flux_x output

      // Output flux_y:
      if (vel_dim[id] > 1) {
	if (method[id] >=0) { //PIC gather:
	  FTYPE scale=gridsize*pic[id].n0avg/npd[id]*dy/dt;
	  flux=0.;
	  gather_weight(flux, pic[id], pic[id].vy, gmin, gmax, gn, wrap, scale);
	  
#ifdef USE_MPI
	  // Sum the flux arrays accross all processors 
	  int mpi_err=MPI_Reduce((void*)&flux(0),(void*)&flux_sum(0), flux.length(), 
			     MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	  if ( mpi_err != MPI_SUCCESS) 
	    mpi_error(mpi_rank, mpi_err,"Reduce flux call in output_fluxes failed");
	  flux = flux_sum;
	  flux /= subdomain.np;
#ifdef USE_DOMAINS
	  pass_sum_guard(flux, nx/nout_avg);
#endif    
#endif    
	  
	  //Print y:
	  if (subdomain.rank == 0) {
	    if (hdf_output_arrays == 0) fwrite(&flux(0), sizeof(OTYPE), static_cast<int>(gridsize), ffluxy[id]);
	    hsize_t dims[1] = {(hsize_t)gridsize};
	    if ( hdf_output_arrays == 1 ) 
	      H5LTmake_dataset(h5_gid, h5_ffluxy[id], 1, dims, H5T_NATIVE_FLOAT, &flux(0));
	    int ngrid[] = {INDICIES(nx,ny,nz)}; 
	    if (hdf_output_arrays == 2)
	      output_collective_array_gen_h5(H5_FID, h5_ffluxy[id], flux, ngrid, nout_avg);
          }
	    
	// } else { // Fluid method
	} else if (method[id] < 0 && method[id] != -4) {
	  if (subdomain.rank == 0) {
	    int ngrid[] = {INDICIES(nx,ny,nz)}; 
	    
#ifdef USE_DOMAINS
	    void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
				INDICIES(int startx, int starty, int startz),
			      INDICIES(int endbeforex, 
				       int endbeforey, 
				       int endbeforez));
	    FArrayND no_guards(FArrayND_ranged &in);
	    FArrayND flux_tmp=no_guards(fspecie[id].vy);
	    sub_array_mult(fspecie[id].den,flux_tmp,
			   INDICIES(0,0,0),INDICIES(nx,ny,nz));
	    if ( hdf_output_arrays == 0 ) 
              output_array(ffluxy[id], flux_tmp,ngrid,nout_avg);
	    else if ( hdf_output_arrays == 1 ) 
              output_array_h5(h5_gid, h5_ffluxy[id], flux_tmp,ngrid,nout_avg);
	    else if ( hdf_output_arrays == 2 ) 
              output_collective_array_h5(H5_FID, h5_ffluxy[id], flux_tmp,ngrid,nout_avg);

#else
	      FArrayND flux_tmp=fspecie[id].vy;
	      flux_tmp *= fspecie[id].den;
	      if ( hdf_output_arrays == 0 ) 
                output_array(ffluxy[id], flux_tmp,ngrid,nout_avg);
	      else if ( hdf_output_arrays == 1 ) 
                output_array(h5_gid, h5_ffluxy[id], flux_tmp,ngrid,nout_avg);
              
#endif
	    }
	}
      } // End of flux_y output
      
      // Output flux_z:
      if (vel_dim[id] > 2) {
	if (method[id] >=0) { //PIC gather:
	  FTYPE scale=gridsize*pic[id].n0avg/npd[id]*dz/dt;
	  flux=0.;
	  gather_weight(flux, pic[id], pic[id].vz, gmin, gmax, gn, wrap, scale);
	  
#ifdef USE_MPI
	  // Sum the flux arrays accross all processors 
	  int mpi_err=MPI_Reduce((void*)&flux(0),(void*)&flux_sum(0), flux.length(), 
			     MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	  if (mpi_err != MPI_SUCCESS) 
	    mpi_error(mpi_rank, mpi_err,"Reduce flux call in output_fluxes failed");
	  flux = flux_sum;
	  flux /= subdomain.np;
#ifdef USE_DOMAINS
	  pass_sum_guard(flux, nx/nout_avg);
#endif    
#endif    
	  //Print particle flux:
	  if (subdomain.rank == 0) {
	    if (hdf_output_arrays == 0) fwrite(&flux(0), sizeof(OTYPE), static_cast<int>(gridsize), ffluxz[id]);
	    hsize_t dims[1] = {(hsize_t)gridsize};
	    if ( hdf_output_arrays == 1 ) 
	      H5LTmake_dataset(h5_gid, h5_ffluxz[id], 1, dims, H5T_NATIVE_FLOAT, &flux(0));
	    int ngrid[] = {INDICIES(nx,ny,nz)}; 
	    if ( hdf_output_arrays == 2 ) 
              output_collective_array_gen_h5(H5_FID, h5_ffluxz[id], flux, ngrid, nout_avg);
          }
	// } else { // Fluid method
	} else if (method[id] < 0 && method[id] != -4) {
	  if (subdomain.rank == 0) {
	    int ngrid[] = {INDICIES(nx,ny,nz)}; 
	    
#ifdef USE_DOMAINS
	    void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
				INDICIES(int startx, int starty, int startz),
				INDICIES(int endbeforex, 
					 int endbeforey, 
					 int endbeforez));
	    FArrayND no_guards(FArrayND_ranged &in);
	    FArrayND flux_tmp=no_guards(fspecie[id].vz);
	    sub_array_mult(fspecie[id].den,flux_tmp,
			   INDICIES(0,0,0),INDICIES(nx,ny,nz));
	    if ( hdf_output_arrays == 0 ) 
              output_array(ffluxz[id], flux_tmp,ngrid,nout_avg);
	    else if ( hdf_output_arrays == 1 ) 
              output_array_h5(h5_gid, h5_ffluxz[id], flux_tmp,ngrid,nout_avg);
	    else if ( hdf_output_arrays == 2 ) 
              output_collective_array_h5(H5_FID, h5_ffluxz[id], flux_tmp,ngrid,nout_avg);

#else
	    FArrayND flux_tmp=fspecie[id].vz;
	    flux_tmp *= fspecie[id].den;
	    if ( hdf_output_arrays == 0 ) 
              output_array(ffluxz[id], flux_tmp,ngrid,nout_avg);
	    else  if ( hdf_output_arrays == 1 ) 
              output_array_h5(h5_gid, h5_ffluxz[id], flux_tmp,ngrid,nout_avg);
#endif
	  }
	}
	
      } // End of flux_z output

    }
      
  }
  
  if (subdomain.rank == 0 && hdf_output_arrays == 1 ) H5Gclose(h5_gid);

  return;
  
}
