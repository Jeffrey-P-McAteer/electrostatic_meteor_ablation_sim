// Output of nvsqr

#include <stdio.h> 
#include <string.h> 

#include "eppic.h"
#include "eppic-mpi.h"

extern hid_t H5_FID;
extern void output_collective_array_gen_h5(hid_t file_id, char *varname, ArrayNd<OTYPE,NDIM> &a, int ngrid[NDIM], int nout_avg);

#undef __FUNCT__
#define __FUNCT__ "output_nvsqr"
void output_nvsqr(particle_dist *pic, int it, char *dir, int it0)
{
  
  extern void gather_weight_sqr(ArrayNd<OTYPE, NDIM> &den, particle_dist &pic,
			    PTYPEAVec &b, FTYPE *xmin, 
			    FTYPE *xmax, int *nmesh, int *wrap, FTYPE scale);

  static FILE *fnvsqrx[MAXDIST], *fnvsqry[MAXDIST], *fnvsqrz[MAXDIST];
  
  extern hid_t h5_gid;
  char h5_fname[256], h5_groupname[256];
  static char h5_fnvsqrx[MAXDIST][256] = {"\0"},
              h5_fnvsqry[MAXDIST][256] = {"\0"},
              h5_fnvsqrz[MAXDIST][256] = {"\0"};
  
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
      	if (method[id2] >= 0) {

	  char name[256];
	  char *openbtype;
	  FILE* fopensafe(char* filename, char* mode, unsigned long skip);
	  if (it == 0) openbtype="wb\0 ";
	  else openbtype="ab\0";
	  if (nvsqr_out_subcycle[id2] > 0) {
	    sprintf(name,"%snvsqrx%d.bin",dir,id2);
	    unsigned long skip=sizeof(OTYPE)*
	      max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
	      ((static_cast<unsigned long>(it)-1)/
	       (nout*nvsqr_out_subcycle[id2])+1);
	    if (hdf_output_arrays == 0) fnvsqrx[id2]  = fopensafe(name,openbtype,skip);
	    sprintf(h5_fnvsqrx[id2], "nvsqrx%d", id2); 
	  }

	  if (vel_dim[id2] >= 2) {
	    if (nvsqr_out_subcycle[id2] > 0) {
	      sprintf(name,"%snvsqry%d.bin",dir,id2);
	      unsigned long skip=sizeof(OTYPE)*
		max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
		((static_cast<unsigned long>(it)-1)/
		 (nout*nvsqr_out_subcycle[id2])+1);
	      if (hdf_output_arrays == 0) fnvsqry[id2]  = fopensafe(name,openbtype,skip);
	      sprintf(h5_fnvsqry[id2], "nvsqry%d", id2);
	    }

	    if (vel_dim[id2] >= 3) {
	      if (nvsqr_out_subcycle[id2] > 0) {
		sprintf(name,"%snvsqrz%d.bin",dir,id2);
		unsigned long skip=sizeof(OTYPE)*
		  max(nx/nout_avg,1)*max(ny/nout_avg,1)*max(nz/nout_avg,1)*
		  ((static_cast<unsigned long>(it)-1)/
		   (nout*nvsqr_out_subcycle[id2])+1);
		if (hdf_output_arrays == 0) fnvsqrz[id2]  = fopensafe(name,openbtype,skip);
		sprintf(h5_fnvsqrz[id2], "nvsqrz%d", id2);
	      }
	    }
	  }
	} // end if (method[id2] >= 0)
      }
    }

  } // end if (first_entry) 

  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0) {
      if (it%(nout*nvsqr_out_subcycle[id]) == 0) {
	  // Define array limits
	  int gn[]={INDICIES(nx/nout_avg,ny/nout_avg,nz/nout_avg)};
	  FTYPE gmin[]={INDICIES(0,0,0)};
	  FTYPE gmax[]={INDICIES((FTYPE)nx,(FTYPE)ny,(FTYPE)nz)};
	  int wrap[]={INDICIES(nx/nout_avg,ny/nout_avg,nz/nout_avg)};
	  
	  // These need to be adjusted if there are subdomains:
	  if (nsubdomains > 1) {
	      gn[0]+=1;
	      wrap[0]+=1; // This eliminates the wrap
	      gmax[0]+=1*nout_avg;
	  }

	  // Calculate total output gridsize
	  FTYPE gridsize=nx/nout_avg;
	  for (int idim=1; idim<ndim; ++idim) gridsize *= (gmax[idim]-gmin[idim])/nout_avg;

	  // Define the array to sum onto
	  ArrayNd<OTYPE,NDIM> nvsqr(gn);

	  // Output Nvsqr_x:
	  FTYPE scale=gridsize*pic[id].n0avg/npd[id]*Sqr(dx/dt);
	  nvsqr=0.;
	  gather_weight_sqr(nvsqr, pic[id], pic[id].vx, gmin, gmax, gn, wrap, scale);

#ifdef USE_MPI
	  // Sum the nvsqr arrays accross all processors 
	  ArrayNd<OTYPE,NDIM> nvsqr_sum(gn);
	  int mpi_err=MPI_Reduce((void*)&nvsqr(0),(void*)&nvsqr_sum(0), nvsqr.length(), 
				 MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	  if ( mpi_err != MPI_SUCCESS) mpi_error(mpi_rank, mpi_err,"Reduce nvsqr call in output_nvsqr failed");
	  nvsqr = nvsqr_sum;
	  nvsqr /= subdomain.np;

#ifdef USE_DOMAINS
	  void pass_sum_guard(ArrayNd<OTYPE,NDIM> &, int);
	  pass_sum_guard(nvsqr, nx/nout_avg);
#endif
#endif    

	  //Print nvsqr_x:
	  if (subdomain.rank == 0) {
            if ( hdf_output_arrays == 0 ) 
	      fwrite(&nvsqr(0), sizeof(OTYPE),  static_cast<int>(gridsize), fnvsqrx[id]);
            else  if ( hdf_output_arrays == 1 ) { 
	      hsize_t dims[1] = {(hsize_t)gridsize};
	      H5LTmake_dataset(h5_gid, h5_fnvsqrx[id], 1, dims, H5T_NATIVE_FLOAT, &nvsqr(0));
	    }
	    else if ( hdf_output_arrays == 2) {
	      int ngrid[] = {INDICIES(nx,ny,nz)}; 
	      output_collective_array_gen_h5(H5_FID, h5_fnvsqrx[id], nvsqr, ngrid, 1);
	    }
          }
	  // End of nvsqr_x output
	  
	  // Output nvsqr_y:
	  if (vel_dim[id] > 1) {
	      FTYPE scale=gridsize*pic[id].n0avg/npd[id]*Sqr(dy/dt);
	      nvsqr=0.;
	      gather_weight_sqr(nvsqr, pic[id], pic[id].vy, gmin, gmax, gn, wrap, scale);
#ifdef USE_MPI
	      // Sum the nvsqr arrays accross all processors 
	      mpi_err=MPI_Reduce((void*)&nvsqr(0),(void*)&nvsqr_sum(0), nvsqr.length(), 
				 MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	      if ( mpi_err != MPI_SUCCESS) mpi_error(mpi_rank, mpi_err,"Reduce nvsqr call in output_nvsqr failed");
	      nvsqr = nvsqr_sum;
	      nvsqr /= subdomain.np;	      
#ifdef USE_DOMAINS
	      pass_sum_guard(nvsqr, nx/nout_avg);
#endif
#endif    

	      /*cout << mpi_rank<<":"<< pic[id].n0avg <<" "<< pic[id].np <<" "<< npd[id]
		   <<" "<< scale <<" "<<nvsqr.sum()<< endl;
	      */
	      //Print nvsqr_y
	      if (subdomain.rank == 0) {
		if ( hdf_output_arrays == 0 ) 
		  fwrite(&nvsqr(0), sizeof(OTYPE),  static_cast<int>(gridsize), fnvsqry[id]);
		else if ( hdf_output_arrays == 1 ) {
		  hsize_t dims[1] = {(hsize_t)gridsize}; 
                  H5LTmake_dataset(h5_gid, h5_fnvsqry[id], 1, dims, H5T_NATIVE_FLOAT, &nvsqr(0));
		}
		else if ( hdf_output_arrays == 2) {
		  int ngrid[] = {INDICIES(nx,ny,nz)}; 
		  output_collective_array_gen_h5(H5_FID, h5_fnvsqry[id], nvsqr, ngrid, 1);
	    }

              }
	  } // End of nvsqr_y output
      
	  // Output nvsqr_z:
	  if (vel_dim[id] > 2) {
	      FTYPE scale=gridsize*pic[id].n0avg/npd[id]*Sqr(dz/dt);
	      nvsqr=0.;
	      gather_weight_sqr(nvsqr, pic[id], pic[id].vz, gmin, gmax, gn, wrap, scale);
	  
#ifdef USE_MPI
	      // Sum the nvsqr arrays accross all processors 
	      mpi_err=MPI_Reduce((void*)&nvsqr(0),(void*)&nvsqr_sum(0), nvsqr.length(), 
				 MPI_OTYPE, MPI_SUM, 0, subdomain.internal_comm);
	      if ( mpi_err != MPI_SUCCESS) mpi_error(mpi_rank, mpi_err,"Reduce nvsqr call in output_nvsqr failed");
	      nvsqr = nvsqr_sum;
	      nvsqr /= subdomain.np;	      
#ifdef USE_DOMAINS
	      pass_sum_guard(nvsqr, nx/nout_avg);
#endif
#endif    

	  //Print nvsqr_z:
	      if (subdomain.rank == 0) {
		if ( hdf_output_arrays == 0 ) 
		  fwrite(&nvsqr(0), sizeof(OTYPE),  static_cast<int>(gridsize), fnvsqrz[id]);
		else if ( hdf_output_arrays == 1 ) {
                  hsize_t dims[1] = {(hsize_t)gridsize};
		  H5LTmake_dataset(h5_gid, h5_fnvsqrz[id], 1, dims, H5T_NATIVE_FLOAT, &nvsqr(0));
		}
		else if ( hdf_output_arrays == 2) {
		  int ngrid[] = {INDICIES(nx,ny,nz)}; 
		  output_collective_array_gen_h5(H5_FID, h5_fnvsqrz[id], nvsqr, ngrid, 1);
		}
	      }
	  } // End of nvsqr_z output
      }
    } // end if (method[id] >= 0)
  }

  if (subdomain.rank == 0) {
    if ( hdf_output_arrays == 1 ) H5Gclose(h5_gid);
  }
  return;
  
}

