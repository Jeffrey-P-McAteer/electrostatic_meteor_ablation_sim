
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "eppic.h"
#include "eppic-types.h"
#include "eppic-mpi.h"
#include <math.h>

void output_array(FILE* fname, FArrayND &a, int ngrid[NDIM], int nout_avg){
  
  // Open a file, fname, if needed, and output an array in single precision 
  // binary, averaging over nout points.

  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM == 2 )
    ny = ngrid[1];
  if ( NDIM == 3 ) {
    ny = ngrid[1];
    nz = ngrid[2];
  }

  if ( fname != 0 ) {
    for (int ix=0; ix<nx; ix += nout_avg)
      for (int iy=0; iy<ny; iy += nout_avg) 
	for (int iz=0; iz<nz; iz += nout_avg) {
	  float f=0.;
	  int navg=0;
	  for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++) 
	    for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++) 
	      for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
		navg++;
		f += a(INDICIES(ix2,iy2,iz2));
	      } 
	  f /= navg;
	  fwrite(&f,sizeof(f),1,fname);
	} 
  }
  
}

void output_array(FILE* fname, ArrayNd_ranged<FTYPE,NDIM> &a,int ngrid[NDIM], int nout_avg){
  // Open a file, fname, if needed, and output an array in single precision 
  // binary, averaging over nout points.
  
  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM == 2 )
    ny = ngrid[1];
  if ( NDIM == 3 ) {
    ny = ngrid[1];
    nz = ngrid[2];
  }
  
  if ( fname != 0 ) {
    for (int ix=0; ix<nx; ix += nout_avg)
      for (int iy=0; iy<ny; iy += nout_avg) 
	for (int iz=0; iz<nz; iz += nout_avg) {
	  float f=0.;
	  int navg=0;
	  for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++) 
	    for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++) 
	      for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
		navg++;
		f += a(INDICIES(ix2,iy2,iz2));
	      } 
	  f /= navg;
	  fwrite(&f,sizeof(f),1,fname);
	} 
  }
  
}

void output_array_h5(hid_t h5_id, char *varname, FArrayND &a, int ngrid[NDIM], int nout_avg){
  // Output an array in single precision binary, averaging over nout points to 
  // an HDF5 file h5_id.

  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM == 2 )
    ny = ngrid[1];
  if ( NDIM == 3 ) {
    ny = ngrid[1];
    nz = ngrid[2];
  }
  
  if ( strlen(varname) > 0 ) {
    int nx2=nx/nout_avg, ny2=ny/nout_avg, nz2=nz/nout_avg;             //output grid size
    typedef ArrayNd<OTYPE,NDIM> OArrayND;
    OArrayND f(INDICIES(nx2,ny2,nz2));                                 //ArrayND array f
    float navg=pow(nout_avg,OTYPE(NDIM));                                    //number of cells to average over

    for (int ix=0; ix<nx; ix += nout_avg)
      for (int iy=0; iy<ny; iy += nout_avg)
        for (int iz=0; iz<nz; iz += nout_avg) {
	  int ix3=ix/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;       //output cell location
          for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++)
            for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++)
              for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
		f(INDICIES(ix3,iy3,iz3)) += a(INDICIES(ix2,iy2,iz2));  //write f
              }

	  f(INDICIES(ix3,iy3,iz3)) /= navg;                            //averaging
        }

    hsize_t dims[NDIM]={INDICIES(nx2,ny2,nz2)};
    hid_t space_id = H5Screate_simple(NDIM, dims, NULL);
    htri_t dset_exists;
    dset_exists = H5Lexists(h5_id, varname, H5P_DEFAULT);
    hid_t data_id;
    if (dset_exists == 0) data_id = H5Dcreate(h5_id, varname, H5T_NATIVE_FLOAT, space_id, 
					      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    else data_id = H5Dopen(h5_id, varname, H5P_DEFAULT);
    H5Dwrite(data_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f.address());
  }
  
}

void output_array_h5(hid_t h5_id, char *varname, ArrayNd_ranged<FTYPE,NDIM> &a,int ngrid[NDIM], int nout_avg){
  // Output an array in single precision binary, averaging over nout points to 
  // an HDF5 file h5_id.
  
  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM == 2 )
    ny = ngrid[1];
  if ( NDIM == 3 ) {
    ny = ngrid[1];
    nz = ngrid[2];
  }
  
  if ( strlen(varname) == 0 ) return;
      
  int nx2=nx/nout_avg, ny2=ny/nout_avg, nz2=nz/nout_avg;             //output grid size
  typedef ArrayNd<OTYPE,NDIM> OArrayND;
  OArrayND f(INDICIES(nx2,ny2,nz2));
  float navg=pow(nout_avg,OTYPE(NDIM));                                    //number of cells to average over
  
  for (int ix=0; ix<nx; ix += nout_avg)
    for (int iy=0; iy<ny; iy += nout_avg)
      for (int iz=0; iz<nz; iz += nout_avg) {
	int ix3=ix/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;       //output cell location
	for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++)
	  for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++)
	    for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
	      f(INDICIES(ix3,iy3,iz3)) += a(INDICIES(ix2,iy2,iz2));  //write f
	    }
	
	f(INDICIES(ix3,iy3,iz3)) /= navg;                            //averaging
      }
  
  hsize_t dims[NDIM]={INDICIES(nx2,ny2,nz2)};
  hid_t space_id = H5Screate_simple(NDIM, dims, NULL);
  htri_t dset_exists;
  dset_exists = H5Lexists(h5_id, varname, H5P_DEFAULT);
  cout << "dset_exists = " << dset_exists << "\n";
  hid_t data_id;
  if (dset_exists == 0) data_id = H5Dcreate(h5_id, varname, H5T_NATIVE_FLOAT, space_id, 
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  else data_id = H5Dopen(h5_id, varname, H5P_DEFAULT);
  H5Dwrite(data_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f.address());
  
}

void output_collective_array_h5(hid_t file_id, char *varname, FArrayND &a, int ngrid[NDIM],
                                int nout_avg) {

  hid_t   dset_id;       // Dataset identifier
  hid_t   filespace;     // Dataspace id in file
  hid_t   memspace;      // Dataspace id in memory.
  hid_t   plist_id;      // Property List id

  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM > 1 ) ny = ngrid[1];
  if ( NDIM > 2 ) nz = ngrid[2];
  
  hsize_t sx = maxInt(nx/nout_avg, 1);
  hsize_t sy = maxInt(ny/nout_avg, 1);
  hsize_t sz = maxInt(nz/nout_avg, 1);
    
  hsize_t tx        = nsubdomains*sx;       // the total number of cells in x-direction.
  hsize_t szL[NDIM] = {INDICIES(sx,sy,sz)}; // the local size of the subdomain
  hsize_t szG[NDIM] = {INDICIES(tx,sy,sz)}; // the global domain size.
  hsize_t starts[3];

  starts[0] = sx*subdomain.id_number;
  starts[1] = 0;
  starts[2] = 0;

  if ( strlen(varname) == 0 ) return;
  
  typedef ArrayNd<OTYPE,NDIM> OArrayNd;
  OArrayNd f(INDICIES(sx,sy,sz));
  float navg=pow(nout_avg,OTYPE(NDIM));                                    //number of cells to average over

  for (int ix=0; ix<nx; ix += nout_avg)
    for (int iy=0; iy<ny; iy += nout_avg) 
      for (int iz=0; iz<nz; iz += nout_avg) {
	int ix3=ix/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;       //output cell location
	for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++) 
	  for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++) 
	    for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
	      f(INDICIES(ix3,iy3,iz3)) += a(INDICIES(ix2,iy2,iz2));  //write f
	    } 
	f(INDICIES(ix3,iy3,iz3)) /= navg;                            //averaging
      }

  filespace = H5Screate_simple(NDIM, szG, NULL);
  memspace  = H5Screate_simple(NDIM, szL, NULL);
  dset_id   = H5Dcreate(file_id, varname, H5T_NATIVE_FLOAT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_id < 0) terminate(-1,"Parallel dataset creation unsuccessful\n");
  H5Sclose(filespace);
  
  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, szL , NULL);
  
  // Tell HDF5 (and MPI I/O) to use collective writes
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // Write the dataset
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
	   plist_id, f.address());

  // Close HDF5 identifiers
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);

}

void output_collective_array_gen_h5(hid_t file_id, char *varname, ArrayNd<OTYPE,NDIM> &a, int ngrid[NDIM], int nout_avg) {
  // this routine is exactly the same as output_collective_array_h5, but is passed an ArrayNd 
  // instead of a FArrayND. Used for output_fluxes, output_vdist, output_nvsqr.
  //  There was probably a better way to do this, but this worked.
  // -LKT 3/6/17

  hid_t   dset_id;       // Dataset identifier
  hid_t   filespace;     // Dataspace id in file
  hid_t   memspace;      // Dataspace id in memory.
  hid_t   plist_id;      // Property List id

  int nx = ngrid[0];
  int ny = 1;
  int nz = 1;
  if ( NDIM > 1 ) ny = ngrid[1];
  if ( NDIM > 2 ) nz = ngrid[2];
  
  hsize_t sx = maxInt(nx/nout_avg, 1);
  hsize_t sy = maxInt(ny/nout_avg, 1);
  hsize_t sz = maxInt(nz/nout_avg, 1);
    
  hsize_t tx        = nsubdomains*sx;       // the total number of cells in x-direction.
  hsize_t szL[NDIM] = {INDICIES(sx,sy,sz)}; // the local size of the subdomain
  hsize_t szG[NDIM] = {INDICIES(tx,sy,sz)}; // the global domain size.
  hsize_t starts[3];

  starts[0] = sx*subdomain.id_number;
  starts[1] = 0;
  starts[2] = 0;

  if ( strlen(varname) == 0 ) return;

  typedef ArrayNd<OTYPE,NDIM> OArrayNd;
  OArrayNd f(INDICIES(sx/nout_avg,sy/nout_avg,sz/nout_avg));
  float navg=pow(nout_avg,OTYPE(NDIM));                          //number of cells to average over

  for (int ix=0; ix<nx; ix += nout_avg)
    for (int iy=0; iy<ny; iy += nout_avg) 
      for (int iz=0; iz<nz; iz += nout_avg) {
	int ix3=ix/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;       //output cell location
	for (int ix2=ix; ix2<ix+nout_avg && ix2 < nx; ix2++) 
	  for (int iy2=iy; iy2<iy+nout_avg && iy2 < ny; iy2++) 
	    for (int iz2=iz; iz2<iz+nout_avg && iz2 < nz; iz2++) {
	      f(INDICIES(ix3,iy3,iz3)) += a(INDICIES(ix2,iy2,iz2));  //write f
	    } 
	f(INDICIES(ix3,iy3,iz3)) /= navg;                            //averaging
      }
  
  filespace = H5Screate_simple(NDIM, szG, NULL);
  memspace  = H5Screate_simple(NDIM, szL, NULL);
  dset_id   = H5Dcreate(file_id, varname, H5T_NATIVE_FLOAT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_id < 0) terminate(-1,"Parallel dataset creation unsuccessful\n");
  H5Sclose(filespace);
  
  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, szL , NULL);
  
  // Tell HDF5 (and MPI I/O) to use collective writes
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // Write the dataset
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
	   plist_id, a.address());

  // Close HDF5 identifiers
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);

}
