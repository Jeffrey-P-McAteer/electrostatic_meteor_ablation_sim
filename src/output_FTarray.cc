#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions

#include "eppic-system.h"
#include "eppic-types.h"
#include "eppic-io.h" // for eppics version of terminate
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#ifdef USE_DOMAINS
#include "classes/realfftwmpi_transpose.h"
typedef realfftwmpi<FTYPE,NDIM> FArrayND_fftw;
#else 
#include <complex>
#include "classes/realfftw2.h"
typedef realfftw2<FTYPE,NDIM> FArrayND_fftw;
#endif

int output_FTarray(FILE *fname, FArrayND &A, eppic_system &sys,
		   FTYPE Afrac, FTYPE kmax)
{
  /* This subroutine takes an array A, copies it to Awork (to add
     necessary memory buffers) Fast Fourier Transforms it, finds the
     highest amplitude mode (excluding k=0), Amax, outputs to fname
     all components of A and their indicies which exceed Afrac*Amax.
  */
  
  /* Local variables */



  static FTYPE size_fact=1.0/(sys.nx*nsubdomains)/sys.ny/sys.nz;
  
  const int nsize[]={INDICIES(sys.nx*nsubdomains,sys.ny,sys.nz)};
#ifdef USE_DOMAINS
  static FArrayND_fftw A_trans(subdomain.neighbor_comm, nsize);
  if (A_trans.local_nx != sys.nx) 
    terminate(-1,"Error: fftw_mpi is not assigning A_trans to have nx=local_nx as required\n");
  if (A_trans.workspace == 0) A_trans.make_workspace(nsize,0.);
#else
  FArrayND A_trans(nsize);
#endif

  // A must have some padding in the last dimension to 
  // allow fftw to perform an inplace transform so
  // Copy A into A_trans, setting all non-overlapping values to zero
  A_trans.cdata=std::complex<FTYPE>(0.,0.);
  for (int ix = 0;  ix < sys.nx; ++ix)
    for (int iy = 0; iy < sys.ny; ++iy)
      for (int iz = 0; iz < sys.nz; ++iz) {
	A_trans.data(INDICIES(ix,iy,iz))=A(INDICIES(ix,iy,iz));
      }
  
  // Transform into k space
  A_trans.transform();
  A_trans.data_transposed *= size_fact;

  // Remove the DC component:
  complex<FTYPE> A_dc;
  int ic=0; // A counter of the number of output

  if (mpi_rank == 0) {
    A_dc=A_trans.cdata(INDICIES(0,0,0));
    A_trans.cdata(INDICIES(0,0,0))= complex<FTYPE>(0.,0.);
  }

  //Set the max extent of the arrays - divide the last dim by 2 (except nky
  int nx_global=sys.nx*nsubdomains;
  int nk[]={1,1,1};
  for (int id=0; id < NDIM; id++) nk[id]=A_trans.cdata.size(id);

  //Obtain the maximum of A^2 values
  FTYPE A2_max=0, A2_cutoff, tmp;
  if (A_trans.cdata.size() > 0) { // Skip if this processor contains no data:
    for (int y = 0; y < nk[0]; ++y)
      for (int x = 0; x < nk[1]; ++x)
	for (int z = 0; z < nk[2]; ++z)
	  if ( (tmp=norm(A_trans.cdata(INDICIES(y,x,z)))) > A2_max) A2_max=tmp;
  }
  A2_cutoff=A2_max*Afrac*Afrac;

  // Pass around A2_cutoff to all domains to obtain the global maximum.
#ifdef USE_MPI
  {
    FTYPE A2_cutoff_all=0;  
    int mpi_err=MPI_Allreduce(&A2_cutoff, &A2_cutoff_all, 1, 
			   MPI_FTYPE, MPI_MAX, subdomain.neighbor_comm);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce call in compress_array");
    A2_cutoff = A2_cutoff_all;
  }
#endif
 
  // Restore the DC component.
  if (mpi_rank == 0) A_trans.cdata(INDICIES(0,0,0))= A_dc; 
  
  FTYPE ksqr, kx, ky, kz; 
  FTYPE kxbase=2*M_PI/(nx_global*sys.dx/2.);
  FTYPE kybase=2*M_PI/(sys.ny*sys.dy/2.);
  FTYPE kzbase=2*M_PI/(sys.nz*sys.dz/2.);
  FTYPE kmax_sqr = kmax*kmax;

  for (int iky=0; iky<nk[0]; iky++) {
    int iky_global=iky+A_trans.y_start_transpose;
    if (iky_global <= sys.ny/2) ky=kybase*iky_global;
    else ky=kybase*(-(sys.ny-iky_global));

    for (int ikx=0; ikx<nk[1]; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));

      for (int ikz=0; ikz<nk[2]; ikz++) {
	if (ikz <= sys.nz/2) kz=kzbase*ikz;
	else kz=kzbase*(-(sys.nz-ikz));
	ksqr=kx*kx+ky*ky+kz*kz;

	if (norm(A_trans.cdata(INDICIES(iky,ikx,ikz))) >= A2_cutoff) 
	  if (kmax == 0.0 || ksqr <= kmax_sqr) {
	    unsigned short int ik[NDIM]={INDICIES(ikx_global,iky_global,ikz)};
	    fwrite(&ik, sizeof(ik[0]), NDIM, fname);
	    // Copy data to single precision
	    float data[2]={real(A_trans.cdata(INDICIES(iky,ikx,ikz))),
			   imag(A_trans.cdata(INDICIES(iky,ikx,ikz)))};
	    fwrite(&data,sizeof(data[0]),2,fname);
	    ic++;
	  }
      }
    }
  }

  // Clear and restore the A_trans array
  A_trans.del_workspace();
  A_trans.ordering=A_trans.NOT_TRANSPOSED;

  // Return the count
  return ic;    

}


int output_FTarray_h5(hid_t h5_id, char *varname, FArrayND &A, 
		      eppic_system &sys, FTYPE Afrac, FTYPE kmax)
{

  /* This subroutine takes an array A, copies it to Awork (to add
     necessary memory buffers) Fast Fourier Transforms it, finds the
     highest amplitude mode (excluding k=0), Amax, outputs to fname
     all components of A and their indicies which exceed Afrac*Amax.

     This routine uses HDF writes.  For each component, the HDF 
     dataset is extended to contain the new data.
  */
  
  static FTYPE size_fact=1.0/(sys.nx*nsubdomains)/sys.ny/sys.nz;
  const int nsize[]={INDICIES(sys.nx*nsubdomains,sys.ny,sys.nz)};

#ifdef USE_DOMAINS
  static FArrayND_fftw A_trans(subdomain.neighbor_comm, nsize);
  if (A_trans.local_nx != sys.nx) 
    terminate(-1,"Error: fftw_mpi is not assigning A_trans to have nx=local_nx as required\n");
  if (A_trans.workspace == 0) A_trans.make_workspace(nsize,0.);
#else
  FArrayND A_trans(nsize);
#endif

  // A must have some padding in the last dimension to 
  // allow fftw to perform an inplace transform so
  // Copy A into A_trans, setting all non-overlapping values to zero
  A_trans.cdata=std::complex<FTYPE>(0.,0.);
  for (int ix = 0;  ix < sys.nx; ++ix)
    for (int iy = 0; iy < sys.ny; ++iy)
      for (int iz = 0; iz < sys.nz; ++iz) {
	A_trans.data(INDICIES(ix,iy,iz))=A(INDICIES(ix,iy,iz));
      }
  
  // Transform into k space
  A_trans.transform();
  A_trans.data_transposed *= size_fact;

  // Remove the DC component:
  complex<FTYPE> A_dc;

  if (mpi_rank == 0) {
    A_dc=A_trans.cdata(INDICIES(0,0,0));
    A_trans.cdata(INDICIES(0,0,0))= complex<FTYPE>(0.,0.);
  }

  //Set the max extent of the arrays - divide the last dim by 2 (except nky
  int nx_global=sys.nx*nsubdomains;
  int nk[]={1,1,1};
  for (int id=0; id < NDIM; id++) nk[id]=A_trans.cdata.size(id);

  //Obtain the maximum of A^2 values
  FTYPE A2_max=0, A2_cutoff, tmp;
  if (A_trans.cdata.size() > 0) { // Skip if this processor contains no data:
    for (int y = 0; y < nk[0]; ++y)
      for (int x = 0; x < nk[1]; ++x)
	for (int z = 0; z < nk[2]; ++z)
	  if ( (tmp=norm(A_trans.cdata(INDICIES(y,x,z)))) > A2_max) A2_max=tmp;
  }

  A2_cutoff=A2_max*Afrac*Afrac;

  // Pass around A2_cutoff to all domains to obtain the global maximum.
#ifdef USE_MPI
  {
    FTYPE A2_cutoff_all;  
    A2_cutoff_all=0;
    int mpi_err=MPI_Allreduce(&A2_cutoff, &A2_cutoff_all, 1, 
			   MPI_FTYPE, MPI_MAX, subdomain.neighbor_comm);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce call in compress_array");
    A2_cutoff = A2_cutoff_all;
  }
#endif
 
  // Restore the DC component.
  if (mpi_rank == 0) A_trans.cdata(INDICIES(0,0,0))= A_dc; 
  
  FTYPE ksqr, kx, ky, kz; 
  FTYPE kxbase=2*M_PI/(nx_global*sys.dx/2.);
  FTYPE kybase=2*M_PI/(sys.ny*sys.dy/2.);
  FTYPE kzbase=2*M_PI/(sys.nz*sys.dz/2.);
  FTYPE kmax_sqr = kmax*kmax;

  // Declare HDF variables
  // writes 2 sets of data: indices in a list of [ikx, iky (,ikz)]
  // and data in a list of [real part, imaginary part]
  hid_t   ind_dataset,     dat_dataset;          //dataset identifiers
  hid_t   ind_dataspace,   dat_dataspace;        //dataspace indentifiers - refers to total size of dataset
  hid_t   ind_filespace,   dat_filespace;        //filespace identifiers - dataspace id in file
  hid_t   ind_memspace,    dat_memspace;         //memoryspace identifiers - dataspace id in memory
  hid_t   ind_plist,       dat_plist;            //dataset property list identifiers
  hsize_t ind_size_tmp[2], dat_size_tmp[2];      //dimensions of extended dataspace
  hsize_t ind_offset[2],   dat_offset[2];        //hyperslab selection parameters    
  hsize_t ind_block[2]={1,NDIM};
  hsize_t dat_block[2]={1,2};                    //...   
  hsize_t count[2]={1,1};
  hsize_t ind_dims[2]={1,NDIM};                  //dataset dimensions
  hsize_t ind_maxdims[2]={H5S_UNLIMITED,NDIM};
  hsize_t dat_dims[2]={1,2};                     //...
  hsize_t dat_maxdims[2]={H5S_UNLIMITED,2};      
  herr_t  status;                                //error check variable - if an HDF function (file open, 
                                                 //  write, etc.) returns a value<0 it failed

  // Create unique index name
  char index_name[80];
  strcpy(index_name,varname);
  strcat(index_name,"_index");

  // Create index dataspace: initialize to 1xNDIM, maximum size UNLIMITEDxNDIM
  // Create index memoryspace: 1xNDIM
  ind_dataspace = H5Screate_simple(2,ind_dims,ind_maxdims);
  ind_memspace  = H5Screate_simple(2,ind_dims,NULL);

  // Create data dataspace: initialize to 1x2, maximum size UNLIMITEDx2
  // Create data memoryspace: 1x2
  dat_dataspace = H5Screate_simple(2,dat_dims,dat_maxdims);
  dat_memspace  = H5Screate_simple(2,dat_dims,NULL);

  // Enable chunking - tells HDF how big each write wil be
  // Necesary to extend datasets
  ind_plist = H5Pcreate(H5P_DATASET_CREATE);
  dat_plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(ind_plist,2,ind_dims); 
  H5Pset_chunk(dat_plist,2,dat_dims);

  ind_offset[1]=0;                                
  dat_offset[1]=0;

  // Create datasets
  ind_dataset = H5Dcreate(h5_id, index_name, H5T_NATIVE_USHORT, ind_dataspace, 
			  H5P_DEFAULT, ind_plist, H5P_DEFAULT);
  dat_dataset = H5Dcreate(h5_id, varname, H5T_NATIVE_FLOAT, dat_dataspace, 
			  H5P_DEFAULT, dat_plist, H5P_DEFAULT);

  // Release dataset property lists & dataspaces
  H5Pclose(ind_plist);
  H5Pclose(dat_plist);
  H5Sclose(ind_dataspace);
  H5Sclose(dat_dataspace);

  int ic=0; // counter
  
  for (int iky=0; iky<nk[0]; iky++) {
    int iky_global=iky+A_trans.y_start_transpose;
    if (iky_global <= sys.ny/2) ky=kybase*iky_global;
    else ky=kybase*(-(sys.ny-iky_global));

    for (int ikx=0; ikx<nk[1]; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));

      for (int ikz=0; ikz<nk[2]; ikz++) {
	if (ikz <= sys.nz/2) kz=kzbase*ikz;
	else kz=kzbase*(-(sys.nz-ikz));
	ksqr=kx*kx+ky*ky+kz*kz;

	if (norm(A_trans.cdata(INDICIES(iky,ikx,ikz))) >= A2_cutoff) 
	  if (kmax == 0.0 || ksqr <= kmax_sqr) {
	    
	    // Index part:
	    // Extend dataset
	    ind_size_tmp[0] = ic+1;
	    ind_size_tmp[1] = NDIM;
	    H5Dset_extent(ind_dataset, ind_size_tmp);

	    // Select filespace
	    ind_filespace = H5Dget_space(ind_dataset);
	    ind_offset[0] = ic;
	    H5Sselect_hyperslab(ind_filespace,H5S_SELECT_SET,ind_offset,NULL,count,ind_block);

	    unsigned short int ik[NDIM] = {INDICIES(ikx_global, iky_global, ikz)};

	    // Write indices
	    H5Dwrite(ind_dataset,H5T_NATIVE_USHORT,ind_memspace,ind_filespace,H5P_DEFAULT,ik);

	    // Release filespace
	    H5Sclose(ind_filespace);

	    // Data part:
	    // Extend dataset
	    dat_size_tmp[0] = ic+1;
	    dat_size_tmp[1] = 2;
	    H5Dset_extent(dat_dataset,dat_size_tmp);

	    // Select filespace
	    dat_filespace = H5Dget_space(dat_dataset);
	    dat_offset[0] = ic;
	    H5Sselect_hyperslab(dat_filespace,H5S_SELECT_SET,dat_offset,NULL,count,dat_block);
	    
	    float data[2] = {real(A_trans.cdata(INDICIES(iky,ikx,ikz))),
			     imag(A_trans.cdata(INDICIES(iky,ikx,ikz)))};
	    
	    // Write data
	    H5Dwrite(dat_dataset,H5T_NATIVE_FLOAT,dat_memspace,dat_filespace,H5P_DEFAULT,data);

	    // Releaes filespace
	    H5Sclose(dat_filespace);
	    ic++;
	  }
      }
    }
  }
  
  //Release datasets & memoryspaces
  H5Dclose(ind_dataset);
  H5Dclose(dat_dataset);  
  H5Sclose(ind_memspace);
  H5Sclose(dat_memspace);

  // Clear and restore the A_trans array
  A_trans.del_workspace();
  A_trans.ordering=A_trans.NOT_TRANSPOSED;

  // Return the count
  return ic;    

}

int output_FTarray_collective_h5(hid_t h5_id, char *varname, FArrayND &A, 
		      eppic_system &sys, FTYPE Afrac, FTYPE kmax)
{

  /* This subroutine takes an array A, copies it to Awork (to add
     necessary memory buffers) Fast Fourier Transforms it, finds the
     highest amplitude mode (excluding k=0), Amax, outputs to fname
     all components of A and their indicies which exceed Afrac*Amax.

     This routine uses parallel HDF writes.  It obtains the components
     of A to be output by each processor, communicates with neignboring
     processors to determine the total number of outputs, and writes to 
     a global file containing data from the entire domain.
  */
  
  static FTYPE size_fact=1.0/(sys.nx*nsubdomains)/sys.ny/sys.nz;
  const int nsize[]={INDICIES(sys.nx*nsubdomains,sys.ny,sys.nz)};

#ifdef USE_DOMAINS
  static FArrayND_fftw A_trans(subdomain.neighbor_comm, nsize);
  if (A_trans.local_nx != sys.nx) 
    terminate(-1,"Error: fftw_mpi is not assigning A_trans to have nx=local_nx as required\n");
  if (A_trans.workspace == 0) A_trans.make_workspace(nsize,0.);
#else
  FArrayND A_trans(nsize);
#endif

  // A must have some padding in the last dimension to 
  // allow fftw to perform an inplace transform so
  // Copy A into A_trans, setting all non-overlapping values to zero
  A_trans.cdata=std::complex<FTYPE>(0.,0.);
  for (int ix = 0;  ix < sys.nx; ++ix)
    for (int iy = 0; iy < sys.ny; ++iy)
      for (int iz = 0; iz < sys.nz; ++iz) {
	A_trans.data(INDICIES(ix,iy,iz))=A(INDICIES(ix,iy,iz));
      }
  
  // Transform into k space
  A_trans.transform();
  A_trans.data_transposed *= size_fact;

  // Remove the DC component:
  complex<FTYPE> A_dc;

  if (mpi_rank == 0) {
    A_dc=A_trans.cdata(INDICIES(0,0,0));
    A_trans.cdata(INDICIES(0,0,0))= complex<FTYPE>(0.,0.);
  }

  //Set the max extent of the arrays - divide the last dim by 2 (except nky
  int nx_global=sys.nx*nsubdomains;
  int nk[]={1,1,1};
  for (int id=0; id < NDIM; id++) nk[id]=A_trans.cdata.size(id);

  //Obtain the maximum of A^2 values
  FTYPE A2_max=0, A2_cutoff, tmp;
  if (A_trans.cdata.size() > 0) { // Skip if this processor contains no data:
    for (int y = 0; y < nk[0]; ++y)
      for (int x = 0; x < nk[1]; ++x)
	for (int z = 0; z < nk[2]; ++z)
	  if ( (tmp=norm(A_trans.cdata(INDICIES(y,x,z)))) > A2_max) A2_max=tmp;
  }

  A2_cutoff=A2_max*Afrac*Afrac;

  // Pass around A2_cutoff to all domains to obtain the global maximum.
#ifdef USE_MPI
  {
    FTYPE A2_cutoff_all;  
    A2_cutoff_all=0;
    int mpi_err=MPI_Allreduce(&A2_cutoff, &A2_cutoff_all, 1, 
			   MPI_FTYPE, MPI_MAX, subdomain.neighbor_comm);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Reduce call in compress_array");
    A2_cutoff = A2_cutoff_all;
  }
#endif
 
  // Restore the DC component.
  if (mpi_rank == 0) A_trans.cdata(INDICIES(0,0,0))= A_dc; 
  
  FTYPE ksqr, kx, ky, kz; 
  FTYPE kxbase=2*M_PI/(nx_global*sys.dx/2.);
  FTYPE kybase=2*M_PI/(sys.ny*sys.dy/2.);
  FTYPE kzbase=2*M_PI/(sys.nz*sys.dz/2.);
  FTYPE kmax_sqr = kmax*kmax;

  // Declare HDF variables
  // writes 2 sets of data: indices in a list of [ikx, iky (,ikz)]
  // and data in a list of [real part, imaginary part]
  hid_t   ind_dataset,     dat_dataset;          //dataset identifiers
  hid_t   ind_dataspace,   dat_dataspace;        //dataspace indentifiers - refers to total size of dataset
  hid_t   ind_memspace,    dat_memspace;         //memoryspace identifiers - dataspace identifier in memory
  hid_t   ind_filespace,   dat_filespace;        //filespace identifiers - dataspace identifier in the file
  hid_t   plist_id;                              //property list identifier
  herr_t  status;                                //error check variable - if an HDF function (file open, 
                                                 //  write, etc.) returns a value<0 it failed

  // Create unique index name
  char index_name[80];
  strcpy(index_name,varname);
  strcat(index_name,"_index");

  long ic=0;  //counter

  //run through to get local array size
  for (int iky=0; iky<nk[0]; iky++) {
    int iky_global=iky+A_trans.y_start_transpose;
    if (iky_global <= sys.ny/2) ky=kybase*iky_global;
    else ky=kybase*(-(sys.ny-iky_global));

    for (int ikx=0; ikx<nk[1]; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));

      for (int ikz=0; ikz<nk[2]; ikz++) {
	if (ikz <= sys.nz/2) kz=kzbase*ikz;
	else kz=kzbase*(-(sys.nz-ikz));
	ksqr=kx*kx+ky*ky+kz*kz;

	if (norm(A_trans.cdata(INDICIES(iky,ikx,ikz))) >= A2_cutoff) 
	  if (kmax == 0.0 || ksqr <= kmax_sqr) {
	    ic++;
	  }
      }
    }
  }
  
  //account for ic=0
  if (ic == 0) ic = 1;

  //declare local arrays
  ArrayNd<OTYPE,2> data_array(ic,2);
  ArrayNd<unsigned short int,2> index_array(ic,NDIM);
  ic = 0; //reset counter
  
  //run through to assign values  
  for (int iky=0; iky<nk[0]; iky++) {
    int iky_global=iky+A_trans.y_start_transpose;
    if (iky_global <= sys.ny/2) ky=kybase*iky_global;
    else ky=kybase*(-(sys.ny-iky_global));
      
    for (int ikx=0; ikx<nk[1]; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));
	
      for (int ikz=0; ikz<nk[2]; ikz++) {
	if (ikz <= sys.nz/2) kz=kzbase*ikz;
	else kz=kzbase*(-(sys.nz-ikz));
	ksqr=kx*kx+ky*ky+kz*kz;
	  
	if (norm(A_trans.cdata(INDICIES(iky,ikx,ikz))) >= A2_cutoff) 
	  if (kmax == 0.0 || ksqr <= kmax_sqr) {
	      
	    unsigned short int ik[NDIM] = {INDICIES(ikx_global, iky_global, ikz)};
	    if (ik[0] > 30000 || ik[1] > 30000 || ik[NDIM-1] > 30000) 
	      cout << "ik " << ik[0] << " close to or exceeds bounds of unsigned short int" << endl;
	    float data[2] = {real(A_trans.cdata(INDICIES(iky,ikx,ikz))),
			     imag(A_trans.cdata(INDICIES(iky,ikx,ikz)))};
	      
	    data_array(ic,0) = data[0];
	    data_array(ic,1) = data[1];
	    index_array(ic,0) = ik[0];
	    index_array(ic,1) = ik[1];
#if NDIM == 3
	    index_array(ic,2) = ik[2];
#endif
	      
	    ic++;
	  }
      }
    }
  }

  //create array nlocal containing ic from each processor
  //& calculate nglobal, the total size of the array
  //  if (ic == 0) ic = 1;
  long nglobal=0;
  long *nlocal;
  int count=1;
  nlocal = (long *)malloc(nsubdomains*count*sizeof(long));
  MPI_Allgather(&ic, count, MPI_LONG, nlocal, count, MPI_LONG, subdomain.neighbor_comm);
  for (int i=0; i<nsubdomains; i++) nglobal += nlocal[i]; 

  // set offset for data write.  each processor should be offset by the total
  // number of data points written by the processors before it.
  // (i.e. each processor's write starts where the last one left off)
  long off=0;
  //    if (subdomain.id_number == 0) off = 0;
  for (int i = 1; i<=subdomain.id_number; i++) off += nlocal[i-1];
  hsize_t offset[2]={(unsigned long long)off,0};
    
  // define local and global array sizes
  hsize_t szL_i[2]={ic,NDIM};
  hsize_t szL_d[2]={ic,2};
  hsize_t szG_i[2]={nglobal,NDIM};
  hsize_t szG_d[2]={nglobal,2};
    
  // Create datasets
  ind_dataspace = H5Screate_simple(2, szG_i, NULL);  //total size of the dataset
  ind_memspace  = H5Screate_simple(2, szL_i, NULL);  //size of each write
  ind_dataset = H5Dcreate(h5_id, index_name, H5T_NATIVE_USHORT, ind_dataspace, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
  dat_dataspace = H5Screate_simple(2, szG_d, NULL);
  dat_memspace  = H5Screate_simple(2, szL_d, NULL);
  dat_dataset = H5Dcreate(h5_id, varname, H5T_NATIVE_FLOAT, dat_dataspace, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
  // Release dataset property lists & dataspaces
  H5Sclose(ind_dataspace);
  H5Sclose(dat_dataspace);
    
  // Select filespace - portion of dataset that will be written to
  ind_filespace = H5Dget_space(ind_dataset);
  H5Sselect_hyperslab(ind_filespace, H5S_SELECT_SET, offset, NULL, szL_i, NULL); 
  dat_filespace = H5Dget_space(dat_dataset);
  H5Sselect_hyperslab(dat_filespace, H5S_SELECT_SET, offset, NULL, szL_d, NULL);
    
  // Tell HDF5 (and MPI I/O) to use collective writes
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
  // Write indices
  H5Dwrite(ind_dataset, H5T_NATIVE_USHORT, ind_memspace, ind_filespace, plist_id, index_array.address());
  H5Dwrite(dat_dataset, H5T_NATIVE_FLOAT, dat_memspace, dat_filespace, plist_id, data_array.address());
    
  // Release filespaces, datasets & memoryspaces
  H5Sclose(ind_filespace);
  H5Sclose(dat_filespace);
  H5Dclose(ind_dataset);
  H5Dclose(dat_dataset);  
  H5Sclose(ind_memspace);
  H5Sclose(dat_memspace);

  // Clear and restore the A_trans array
  A_trans.del_workspace();
  A_trans.ordering=A_trans.NOT_TRANSPOSED;

  // Return the count
  return ic;    
}
