#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-efield.h"
#include <stdio.h> 
#include <string.h> 
#include <cmath> 

void load_phi(char* outdir, field &Efield, int it0)
{

  /* Read phi from existing data */
  /*
  hid_t f_id,d_id,m_id;
  herr_t h5_err;
  hsize_t dims[NDIM] = {INDICIES(Efield.phi_rho.xsize(),
				 Efield.phi_rho.ysize(),
				 Efield.phi_rho.zsize())};
  hsize_t start[NDIM] = {INDICIES(phix_guard_size[0],0,0)};
  hsize_t count[NDIM] = {INDICIES(nx*nsubdomains,ny,nz)};
  int rank = NDIM;
  char path[256];
  
  sprintf(path,"%sparallel/parallel%06d.h5",outdir,it0-1);
  f_id = H5Fopen(path,H5F_ACC_RDONLY,H5P_DEFAULT);
  d_id = H5Dopen(f_id,"/phi",H5P_DEFAULT);
  m_id = H5Screate_simple(rank,dims,NULL);
  h5_err = H5Sselect_hyperslab(m_id,H5S_SELECT_SET,start,NULL,count,NULL);
  
  if (mpi_rank == 0 && h5_err < 0)
    printf("%s:%d WARNING: H5Sselect_hyperslab error %d\n",
	   __func__,__LINE__,h5_err);
  h5_err = H5Dread(d_id,H5T_NATIVE_FLOAT,m_id,H5S_ALL,H5P_DEFAULT,
		   (void*)&(Efield.phi_rho(INDICIES(0,0,0))));
  if (mpi_rank == 0 && h5_err < 0)
    printf("%s:%d WARNING: H5Dread error %d\n",
	   __func__,__LINE__,h5_err);
  */

  char path[256];
  sprintf(path,"%sparallel/parallel%06d.h5",outdir,it0-1);

  herr_t h5_err;
  int rank = NDIM;
  hid_t file_id = H5Fopen(path,H5F_ACC_RDONLY,H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id,"/phi",H5P_DEFAULT);
  hsize_t dims_mem[NDIM] = {INDICIES(Efield.phi_rho.xsize(),
				     Efield.phi_rho.ysize(),
				     Efield.phi_rho.zsize())};
  hsize_t start_mem[NDIM] = {INDICIES(phix_guard_size[0],0,0)};
  hsize_t stride_mem[NDIM] = {INDICIES(nout_avg,nout_avg,nout_avg)};
  hsize_t count[NDIM] = {INDICIES(nx/nout_avg,ny/nout_avg,nz/nout_avg)};

  hsize_t dims_dat[NDIM] = {INDICIES((nx*nsubdomains)/nout_avg,ny/nout_avg,nz/nout_avg)};
  hsize_t start_dat[NDIM] = {INDICIES(subdomain.id_number*nx/nout_avg,0,0)};

  hid_t memspace_id = H5Screate_simple(rank,dims_mem,NULL);
  hid_t dataspace_id = H5Screate_simple(rank,dims_dat,NULL);

  h5_err = H5Sselect_hyperslab(memspace_id,H5S_SELECT_SET,start_mem,stride_mem,count,NULL);
  h5_err = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET,start_dat,NULL,count,NULL);
  h5_err = H5Dread(data_id,H5T_NATIVE_FLOAT,memspace_id,dataspace_id,H5P_DEFAULT,
		   (void*)&(Efield.phi_rho(INDICIES(0,0,0))));

}
