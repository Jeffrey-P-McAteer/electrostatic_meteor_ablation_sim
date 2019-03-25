// Output the distribution in velocity space of (marker) particles

#include <stdio.h> 
#include <string.h> 

#include "eppic.h"
#include "eppic-mpi.h"

extern hid_t H5_FID;
extern void output_collective_array_gen_h5(hid_t file_id, char *varname, ArrayNd<OTYPE,NDIM> &a, int ngrid[NDIM], int nout_avg);

FILE* fopensafe(char* filename, char* mode);
FILE* fopensafe(char* filename, char* mode, unsigned long skip);

void gather_scale(ArrayNd<OTYPE, NDIM> &den, PTYPEAVec *a[], FTYPE *xmin, 
		  FTYPE *xmax, int *nmesh, int *iwrap, FTYPE scale, int np, PTYPEAVec xmissing);
void gather_scale(ArrayNd<OTYPE, 3> &den, PTYPEAVec *a[], FTYPE *xmin, 
		  FTYPE *xmax, int *nmesh, int *iwrap, FTYPE scale, int np, PTYPEAVec xmissing);

void output_vdist(particle_dist *pic, int it, char *dir, int it0)
{
  static FILE *fvdist[MAXDIST];
  static int first_entry=TRUE;
  
  extern hid_t h5_gid;
  char h5_fname[256], h5_groupname[256];
  static char h5_fvdist[MAXDIST][256] = {"\0"};
  
  if (first_entry) {
    first_entry=FALSE;

    if (subdomain.rank==0) {
      // Open the files 
      char name[256];
      char *openbtype;

      if (it == 0) openbtype="wb\0 ";
      else openbtype="ab\0";

      for (int id=0; id<ndist; ++id) {
	if (method[id] != -4) {
	  if (vdist_out_subcycle[id] > 0) {
	    unsigned long asize=sizeof(float) * pnvx[id] * pnvy[id] * pnvz[id];
	    if (hdf_output_arrays == 0) {
	      sprintf(name,"%svdist%d.bin",dir,id);
	      fvdist[id]=fopensafe(name, openbtype, 
				   asize*((it-1)/(nout*vdist_out_subcycle[id])+1));
	    }
	    sprintf(h5_fvdist[id], "vdist%d", id);
	  }
	}
      }
    }
  } // End of if (first_entry)
  
  if (subdomain.rank == 0 && hdf_output_arrays == 1 ) {
    sprintf(h5_groupname, "/time_%06d\0", it);
    h5_gid = H5Gopen(H5_FID, h5_groupname, H5P_DEFAULT);
  }
	  
  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0) {
      if (it%(nout*vdist_out_subcycle[id]) == 0) {
	// Create arrays to control gather_scale:
	if (vel_dim[id] == NDIM) {
	  ArrayNd<OTYPE, NDIM> vdist(INDICIES(pnvx[id],pnvy[id],pnvz[id]));
	  PTYPEAVec *a[NDIM]={INDICIES(&pic[id].vx,&pic[id].vy,&pic[id].vz)};
	  FTYPE vxmin[NDIM]=
	    {INDICIES(pvxmin[id]*dt/dx, pvymin[id]*dt/dy, pvzmin[id]*dt/dz)};
	  FTYPE vxmax[NDIM]=
	    {INDICIES(pvxmax[id]*dt/dx, pvymax[id]*dt/dy, pvzmax[id]*dt/dz)};
	  int nmesh[NDIM]={INDICIES(pnvx[id], pnvy[id], pnvz[id])};
	  int iwrap[NDIM]={INDICIES(nmesh[0],nmesh[1],nmesh[2])};
	  FTYPE fscale=pic[id].n0avg/pic[id].np;
	  for (int idim=0; idim<vel_dim[id]; idim++) fscale *= nmesh[idim];
	
          // Changed by Glenn 5/8/2018 to deal with changing np from injecting/open boundaries
          if (unscale_density){
            gather_scale(vdist, a, vxmin, vxmax, nmesh, iwrap, 1.0, pic[id].np, pic[id].x );
          }
          else{
            gather_scale(vdist, a, vxmin, vxmax, nmesh, iwrap, fscale, pic[id].np, pic[id].x );
          }
	
#ifdef USE_MPI
	  // Gather from all processors in the domain
	  ArrayNd<OTYPE, NDIM> vdist_all(INDICIES(pnvx[id],pnvy[id],pnvz[id]));
	  int mpi_err=MPI_Reduce(vdist.address(),
				 vdist_all.address(),
				 vdist.length(),MPI_OTYPE,MPI_SUM,0,
				 subdomain.internal_comm);
	  if ( mpi_err != MPI_SUCCESS) 
	    mpi_error(mpi_rank, mpi_err,"Reduce vdist call in output");
	  vdist_all /= subdomain.np;
	  vdist=vdist_all;
#endif

	  if (subdomain.rank == 0) {
	    if ( hdf_output_arrays == 0 ) 
	      fwrite(vdist.address(),sizeof(OTYPE),vdist.length(),fvdist[id]);
	    else if ( hdf_output_arrays == 1 ) {
	      hsize_t dims[1] = {(hsize_t)vdist.length()};
	      H5LTmake_dataset(h5_gid, h5_fvdist[id], 1, dims, H5T_NATIVE_FLOAT, vdist.address());
	    }
	    else if ( hdf_output_arrays ==2) {
	      int ngrid[] = {INDICIES(pnvx[id],pnvy[id],pnvz[id])};
	      output_collective_array_gen_h5(H5_FID, h5_fvdist[id], vdist, ngrid, 1);
	    }
	  }
      } else if (vel_dim[id] == 3) {
	  ArrayNd<OTYPE, 3> vdist(pnvx[id],pnvy[id],pnvz[id]);
	  PTYPEAVec *a[3]={&pic[id].vx,&pic[id].vy,&pic[id].vz};
	  FTYPE vxmin[3]={pvxmin[id]*dt/dx, pvymin[id]*dt/dy, pvzmin[id]*dt/dz};
	  FTYPE vxmax[3]={pvxmax[id]*dt/dx, pvymax[id]*dt/dy, pvzmax[id]*dt/dz};
	  int nmesh[3]={pnvx[id], pnvy[id], pnvz[id]};
	  int iwrap[3]={nmesh[0],nmesh[1],nmesh[2]};
	  FTYPE fscale=pic[id].n0avg/pic[id].np;
	  for (int idim=0; idim<vel_dim[id]; idim++) fscale *= nmesh[idim];
          // Changed by Glenn 5/8/2018 to deal with changing np from injecting/open boundaries
          if (unscale_density){
            gather_scale(vdist, a, vxmin, vxmax, nmesh, iwrap, 1.0,  pic[id].np, pic[id].x);
          }
          else{
            gather_scale(vdist, a, vxmin, vxmax, nmesh, iwrap, fscale,  pic[id].np, pic[id].x);
          }

#ifdef USE_MPI
	  // Gather from all processors in the domain
	  ArrayNd<OTYPE, 3> vdist_all(pnvx[id],pnvy[id],pnvz[id]);
	  int mpi_err=MPI_Reduce(vdist.address(),
				 vdist_all.address(),
				 vdist.length(),MPI_OTYPE,MPI_SUM,0,
				 subdomain.internal_comm);
	  if ( mpi_err != MPI_SUCCESS) 
	    mpi_error(mpi_rank, mpi_err,"Reduce vdist call in output");
	  vdist_all /= subdomain.np;
	  vdist=vdist_all;
#endif
	  if (subdomain.rank == 0) {
	    if ( hdf_output_arrays == 0 ) 
	      fwrite(vdist.address(),sizeof(OTYPE),vdist.length(),fvdist[id]);
	    else if ( hdf_output_arrays == 1 ) {
	      hsize_t dims[1] = {(hsize_t)vdist.length()};
	      H5LTmake_dataset(h5_gid, h5_fvdist[id], 1, dims, H5T_NATIVE_FLOAT, vdist.address());
	    }
	    else if ( hdf_output_arrays ==2) {
	      int ngrid[] = {INDICIES(pnvx[id],pnvy[id],pnvz[id])}; 
	      output_collective_array_gen_h5(H5_FID, h5_fvdist[id], vdist, ngrid, 1);
	    }
	  }

       } else {
 	terminate(-1,"output_vdist does not currently support this dimensional output");
       }
        
      }
    }  //if (method[id] >= 0)
  } //  for (id...
 
  if (subdomain.rank == 0) {
    if ( hdf_output_arrays == 1 ) H5Gclose(h5_gid);
  } 
  
}
