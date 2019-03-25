/* ArrayNd_ranged<Type,DIM>:  A simple N-Dimensional array class - 
   designed to be entirely inline. 

   Unlike ArrayNd, this array class can start at non zero integer indicies.

   It is an inherited class from ArrayNd

   Written by Meers Oppenheim starting in 2006
*/

#ifndef ArrayNd_ranged_H 
#define ArrayNd_ranged_H
//
//#ifdef DEBUG
#include <iostream>
#include <sstream>
#include <string.h>
#include <csignal>
#include <hdf5.h>
#include <hdf5_hl.h>
//#endif

#include "ArrayNd.h"
//
#if NDIM == 1
#define INDICIES(ix,iy,iz) ix
#elif NDIM == 2
#define INDICIES(ix,iy,iz) ix, iy
#else
#define INDICIES(ix,iy,iz) ix, iy, iz
#endif

using namespace std;

template <class TYPE,int A_NDIM>
class ArrayNd_ranged : public ArrayNd<TYPE, A_NDIM>
{

 protected:
  //Inherits the data from ArrayND, plus needs a starting index array
  int  nstart[A_NDIM];
  // Data ptr adjusted so that data(nstart) points at the first element
  // data_adj = data - index, index = -1*nelem[1] = -1*nz*ny
  // make data(nstart) point at the beginning data address: data_adj(1) = data(0)
  TYPE *data_adj;
 public:

  using ArrayNd<TYPE,A_NDIM>::operator=;
  using ArrayNd<TYPE,A_NDIM>::nelem;
  using ArrayNd<TYPE,A_NDIM>::nsize;
  using ArrayNd<TYPE,A_NDIM>::data;
  using ArrayNd<TYPE,A_NDIM>::check_error;
  
  //  Constructors:
  ArrayNd_ranged() : ArrayNd<TYPE,A_NDIM> () {
    for (int i=0; i < A_NDIM; i++) nstart[i]=0;
    data_adj=0;
  }
    
    ArrayNd_ranged(const int in_start[], const index_type in_size[], TYPE x=0, TYPE *workspace=0) 
      : ArrayNd<TYPE,A_NDIM>() {
      /*      : ArrayNd<TYPE,A_NDIM>(in_size, x, workspace) {*/
      build_ArrayNd_ranged(in_start,in_size,x,workspace);
    }
      
      
      ArrayNd_ranged(const int n_in[A_NDIM][2], TYPE x=0, TYPE *workspace=0)
      // n_in should be a 2D array where the first element of each dim 
    // is the starting point and the second is the ending point.
    {
      index_type nsize_tmp[A_NDIM];
      for (index_type i=0; i<A_NDIM; i++){
	nsize_tmp[i]=n_in[i][1]-n_in[i][0]+1;
	nstart[i]=n_in[i][0];
      }
      //      *this=ArrayNd_ranged(nstart, nsize_tmp, x, workspace);
      //ArrayNd<TYPE,A_NDIM>::build_ArrayNd(nsize_tmp, x, workspace);
      build_ArrayNd_ranged(nstart,nsize_tmp,x,workspace);
    }

  void build_ArrayNd_ranged(const int n_begin[], const index_type n_in[], TYPE x=0, TYPE *workspace=0) 
    
    {
      //ArrayNd<TYPE,A_NDIM>::build_ArrayNd(n_in, x, workspace);
      ArrayNd<TYPE,A_NDIM>::build_ArrayNd(n_in, x, workspace);
      for (index_type i=0; i<A_NDIM; i++) nstart[i]=n_begin[i];
      // make data(nstart) point at the beginning data address:
      int index=nstart[A_NDIM-1];
      for (index_type i=0; i<A_NDIM-1; i++) index += nstart[i]*nelem[i+1];
      data_adj = data - index;
    }

  
  // Copy Constructor
  ArrayNd_ranged(const ArrayNd_ranged &in) : ArrayNd<TYPE,A_NDIM>(in){
    for (index_type i=0; i<A_NDIM; i++) nstart[i]=in.nstart[i];
    data_adj=in.data_adj;
  };

  
  //  Assignment
    /*
  inline const ArrayNd_ranged& operator= (TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i]=x;
    return *this;
  };
    */

    /* Copy from ArrayNd, still a work in progress */
    /*
      const ArrayNd_ranged& operator= (const ArrayNd<TYPE,A_NDIM> &in) {
    // If they're the same size copy the data
    //otherwise, reinitialize:
    index_type old_nelem=nelem[0];
    for (index_type i=0; i<A_NDIM; i++) {
      nsize[i]=in.size(i);
    }
    nelem[A_NDIM-1]=nsize[A_NDIM-1];
    for (int id = A_NDIM-2; id>=0; id--)  nelem[id]=nsize[id]*nelem[id+1];

    if (old_nelem != in.length()) {
      delete [] data;
      data_adj = data;
      if  (nelem[0] > 0) data=new TYPE[nelem[0]];
      else data=0;
    }
    TYPE *in_data = in.address();
    for (index_type id=0;id<nelem[0];id++) data[id]=in_data[id];
    return *this;
  };
    */


    
  const ArrayNd_ranged& operator= (const ArrayNd_ranged &in) {
    delete [] data;
    build_ArrayNd_ranged(in.nstart,in.nsize,0.,0);
    for (index_type id=0;id<nelem[0];id++) data[id]=in.data[id];
    return *this;
  }


  //  Destructor
  ~ArrayNd_ranged(){
    data_adj = 0;
    }


  //  Access and setting -- Speed is everything
  inline TYPE& operator() (int n0){   
    // One can always access any dimension ArrayNd with a single index
#ifdef DEBUG
    // Need to zero all dimensions past the first because check bounds needs the ND Array
    const int coord[]={n0,0,0,0,0,0};  //Works only to 6D.
    check_bounds(coord);
#endif
    return data_adj[n0];
  }


  inline TYPE& operator[] (int n0){   
    // One can always access any dimension ArrayNd with a single index
#ifdef DEBUG
    check_bounds(n0);
#endif
    return data_adj[n0];
  }

  inline TYPE& operator() (int n0, int n1){   
    // 2-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 2,
		"operator(int n0, int n1) cannot be applied for A_NDIM<2");
    const int coord[]={n0,n1};
    check_bounds(coord);
#endif
    return data_adj[n0*nelem[1]+n1];
  }

  inline TYPE& operator() (int n0, int n1, int n2){   
    // 3-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 3,
		"operator(int n0, int n1, int n2) cannot be applied for A_NDIM<3");
    const int coord[]={n0,n1,n2};
    check_bounds(coord);
#endif
    return data_adj[n0*nelem[1]+n1*nelem[2]+n2];
  }

  inline TYPE& operator() (int n0, int n1, int n2, int n3){   
    // 4-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 4,
		"operator(int n0, int n1, int n2, int n3) cannot be applied for A_NDIM<4");
    const int coord[]={n0,n1,n2,n3};
    check_bounds(coord);
#endif
    return data_adj[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3];
  }
  
  inline TYPE& operator() (int n0, int n1, int n2, int n3, int n4){   
    // 5-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 5,
		"operator(int n0, int n1, int n2, int n3, int n4) cannot be applied for A_NDIM<5");
    const int coord[]={n0,n1,n2,n3,n4};
    check_bounds(coord);
#endif
    return data_adj[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3*nelem[4]+n4];
  }
  
  inline TYPE& operator() (int n0, int n1, int n2, int n3, int n4, int n5){   
    // 6-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 6,
		"operator(int n0, int n1, int n2, int n3, int n4, int n5) cannot be applied for A_NDIM<6");
    const int coord[]={n0,n1,n2,n3,n4,n5};
    check_bounds(coord);
#endif
    return data_adj[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3*nelem[4]+n4*nelem[5]+n5];
  }

  // Generalize access operator with and index argument
  inline TYPE& operator() (int ind[]) {
    int index=ind[A_NDIM-1];  // The last term need not be multiplied.
    for (index_type i=0; i<A_NDIM-1; i++) index += ind[i]*nelem[i+1];
    return operator[](index);
  }

  // To pass matricies to canned routines one sometimes needs the address:
  TYPE* address(int const  ind[]) {
#ifdef DEBUG
    int coord[A_NDIM];
    for (int i=0;i<A_NDIM;i++) coord[i] = ind[i];
    check_bounds(coord);
#endif

    int index=ind[A_NDIM-1];  
    for (index_type i=0; i<A_NDIM-1; i++) index += ind[i]*nelem[i+1];
    return data_adj+index;
  }

  inline index_type xstart() const {return nstart[0];};
  inline index_type ystart() const {return nstart[1];};
  inline index_type zstart() const {return nstart[2];};


  /*
    void assign_values(ArrayNd<TYPE,NDIM> &in, TYPE val_default) {
    // Make the unassined values be val_default
    for (index_type i=0;i<nelem[0];i++) data[i]=val_default;

    // The following code loops through all elements of in 
    // for any dimension array and 
    // assigns the same element number of in to data
    int ind[A_NDIM], index_local;
    for (index_type j=0; j<A_NDIM-1; j++) ind[j]=0;
    ind[A_NDIM-1]=-1;
    for (index_type i=0;i<in.length();i++) {
      //Advance the local index counter based on the index of in
      int idim=A_NDIM-1;
      while ((++ind[idim]) == in.size(idim)) {
	ind[idim]=0;
	idim--;
	check_error(idim < 0,"ArrayNd assign_values input array too small");
      }
      index_local=ind[A_NDIM-1];
      for (index_type j=0; j<A_NDIM-1; j++) index_local += ind[j]*nelem[j+1];
      data[index_local]=in(i);
    }
    }
  */

 protected:

  void error_msg(string complaint) const{
    cout << "\nError: Problem in ArrayNd_ranged\n"
	 << "\t" << complaint << endl;
    cerr << "\nError: Problem in ArrayNd_ranged\n"
	 << "\t" << complaint << endl;

    raise(SIGINT);
    exit(1);
  }

  inline void check_bounds(const int n0) const{
    if (n0>=int(nelem[0])+nstart[0] || n0<nstart[0]) {
      ostringstream msg;
      msg <<"subscript # 1 (of 1D array) = "<<n0<<" out of range ("
	  <<nstart[0]<<","<<nelem[0]+nsize[0]-1<<"]\n";
      error_msg(msg.str());
    }
  }

  inline void check_bounds(const int n0[]) const{
    for (int ii = 0; ii<A_NDIM; ii++) {
      if (n0[ii]>=int(nsize[ii])+nstart[ii] || n0[ii]<nstart[ii]) {
	ostringstream msg;
	msg <<"subscript # "<<ii<<" of "<<A_NDIM<<" = "<<n0[ii]<<" out of range ("
	    <<nstart[ii]<<","<<nstart[ii]+nsize[ii]-1<<"]\n";
	error_msg(msg.str());
      }
    }
  }

 public:  

  inline index_type start(index_type i) const {return nstart[i];};

  
  friend ostream &operator<<(ostream& ostr, const ArrayNd_ranged& x){

    ostr << "Array( ";
    for (int nd=0; nd<A_NDIM; nd++) 
      ostr << "\t" << x.nstart[nd] << ":"<<x.nsize[nd]-x.nstart[nd]-1;     
    ostr << ")=" << endl;
    
    for (int i=0; i<x.nelem[0]; i++) {
      ostr << i << " \t" ;
      for (int nd=0; nd<A_NDIM-1; nd++){
	int ic=(i%x.nelem[nd])/x.nelem[nd+1];
	ostr << ic + x.nstart[nd] << " \t" ;
      }
      if (A_NDIM > 1) 	ostr <<  i%x.nelem[A_NDIM-1] + x.nstart[A_NDIM-1] << " \t" ;
      ostr << x.data[i] << endl;
    }
    
    return ostr;
    
  };

  // A simple output to a file:
  int output(const char *filename) {
    ofstream out;
    out.open( filename);
    out << *this;
    out.close();
    return 0; 
  }; 

  /* Corresponding input from a file made by output: */
  int input(const char *filename) {
    TYPE *workspace=0;
    ifstream in;
    in.open( filename);
    in.ignore(6);
    int startx; in >> startx;
    in.ignore(1); /* skip ':' */
    int nx; in >> nx;
    int starty; in >> starty;
    in.ignore(1); /* skip ':' */
    int ny; in >> ny;
    if (in.peek() == ')') {
      index_type n_in2[2]={nx,ny};
      index_type n_start2[2] = {startx,starty};
      //ArrayNd<TYPE,A_NDIM>::build_ArrayNd(n_in2, 0, workspace);
      build_ArrayNd_ranged(n_start2,n_in2,0,workspace);
    } else { /* assume 3D */
      int startz; in >> startz;
      in.ignore(1); /* skip ':' */
      int nz; in >> nz;
      index_type n_in3[3]={nx,ny,nz};
      index_type n_start3[3] = {startx,starty,startz};
      //ArrayNd<TYPE,A_NDIM>::build_ArrayNd(n_in3, 0, workspace);
      build_ArrayNd_ranged(n_start3,n_in3,0,workspace);
    }
    in.ignore(2); /* skip ')\n' */
    for (int idata=0; idata<nelem[0]; idata++) {
      int count,ix,iy, iz;
      in >> count;
      in >> ix; 
      in  >> iy; 
      if (A_NDIM==3) in  >> iz; 
      in >> data[idata];
    }
    in.close();
    return 0; 
  }; 



  void bin_output_ghost(FILE* fname, int nghost_cells[2]=0, int nout_avg=1){
    // Open a file, fname, if needed, and output an array in single precision 
    // binary, averaging over nout points. nghost_cells defines the number
    // of ghost cells to the left (0) and right (1); these cells are not output
    
    
    int array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix]=1;
    int array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix]=1;
    
    for (index_type i=0; i<A_NDIM; i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i+1];
    
    if ( fname != 0 ) {
      for (int ix=nghost_cells[0]; ix<array_size[0]-nghost_cells[1]; 
	   ix += nout_avg)
	for (int iy=0; iy<array_size[1]; iy += nout_avg) 
	  for (int iz=0; iz<array_size[2]; iz += nout_avg) 
	    for (int iu=0; iu<array_size[3]; iu += nout_avg) 
	      for (int iv=0; iv<array_size[4]; iv += nout_avg) 
		for (int iw=0; iw<array_size[5]; iw += nout_avg) {
		  float f=0.;
		  int navg=0;
		  for (int ix2=ix; ix2<ix+nout_avg && 
			 ix2 < array_size[0]-nghost_cells[1]; ix2++) 
		    for (int iy2=iy; iy2<iy+nout_avg && 
			   iy2 < array_size[1]; iy2++) 
		      for (int iz2=iz; iz2<iz+nout_avg && 
			     iz2 < array_size[2]; iz2++) 
			for (int iu2=iu; iu2<iu+nout_avg && 
			       iu2 < array_size[3]; iu2++) 
			  for (int iv2=iv; iv2<iv+nout_avg && 
				 iv2 < array_size[4]; iv2++) 
			    for (int iw2=iw; iw2<iw+nout_avg && 
				   iw2 < array_size[2]; iw2++) {
			      navg++;
			      f += data[ix2*array_step[0]+
					iy2*array_step[1]+
					iz2*array_step[2]+
					iu2*array_step[3]+
					iv2*array_step[4]+
					iw2*array_step[5]];
			    } 
		  f /= navg;
		  fwrite(&f,sizeof(f),1,fname);
		} 
    }
    
  }
  
  void bin_output_ghost_h5(hid_t h5_id, char *name, int nghost_cells[2]=0, int nout_avg=1){
    // Open a file, fname, if needed, and output an array in single precision 
    // binary, averaging over nout points. nghost_cells defines the number
    // of ghost cells to the left (0) and right (1); these cells are not output
    
    int array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix]=1;
    int array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix]=1;
    
    for (index_type i=0; i<A_NDIM; i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i+1];
    
    int sx = array_size[0]/nout_avg-nghost_cells[0];
    int sy = array_size[1]/nout_avg;
    int sz = array_size[2]/nout_avg;

    hsize_t dims[A_NDIM]={INDICIES(sx,sy,sz)};

    if ( strlen(name) > 0 ) {
      typedef ArrayNd<OTYPE,NDIM> OArrayNd;
      OArrayNd f(INDICIES(sx,sy,sz));
      float navg=pow(nout_avg,OTYPE(NDIM));                                    //number of cells to average over
      
      for (int ix=nghost_cells[0]; ix<array_size[0]-nghost_cells[1]; 
	   ix += nout_avg)
	for (int iy=0; iy<array_size[1]; iy += nout_avg) 
	  for (int iz=0; iz<array_size[2]; iz += nout_avg) 
	    for (int iu=0; iu<array_size[3]; iu += nout_avg) 
	      for (int iv=0; iv<array_size[4]; iv += nout_avg) 
		for (int iw=0; iw<array_size[5]; iw += nout_avg) {
		  //output cell location
		  int ix3=(ix-nghost_cells[0])/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;
		  for (int ix2=ix; ix2<ix+nout_avg && 
			 ix2 < array_size[0]-nghost_cells[1]; ix2++) 
		    for (int iy2=iy; iy2<iy+nout_avg && 
			   iy2 < array_size[1]; iy2++) 
		      for (int iz2=iz; iz2<iz+nout_avg && 
			     iz2 < array_size[2]; iz2++) 
			for (int iu2=iu; iu2<iu+nout_avg && 
			       iu2 < array_size[3]; iu2++) 
			  for (int iv2=iv; iv2<iv+nout_avg && 
				 iv2 < array_size[4]; iv2++) 
			    for (int iw2=iw; iw2<iw+nout_avg && 
				   iw2 < array_size[2]; iw2++) {
			      f(INDICIES(ix3,iy3,iz3)) += data[ix2*array_step[0]+
							       iy2*array_step[1]+
							       iz2*array_step[2]+
							       iu2*array_step[3]+
							       iv2*array_step[4]+
							       iw2*array_step[5]];
			    }
		  f (INDICIES(ix3,iy3,iz3)) /= navg;
		}
      hid_t space_id = H5Screate_simple(NDIM, dims, NULL);
      hid_t data_id = H5Dcreate(h5_id, name, H5T_NATIVE_FLOAT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(data_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f.address());
      
    }
    
  }
  void bin_output_ghost_collective_h5(hid_t file_id, char *name, int nsubdomains, int id_number,
				      int nghost_cells[2]=0, int nout_avg=1){
    // output an array in single precision, averaging over nout points. nghost_cells defines the number
    // of ghost cells to the left (0) and right (1); these cells are not output
    
    hid_t   dset_id;       //Dataset identifier
    hid_t   filespace;     //Dataspace id in file
    hid_t   memspace;      //Dataspace id in memory.
    hid_t   plist_id;      //Property List id
    
    hsize_t szL[6];        for(int ix=0;ix<6;ix++) szL[ix]        = 1;
    hsize_t array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix] = 1;
    hsize_t array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix] = 1;
    hsize_t szG[6], starts[6];

    for (index_type i=0; i<A_NDIM;   i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i+1];
    for (index_type i=1; i<A_NDIM;   i++) szL[i]        = maxInt(nsize[i]/nout_avg,1);
    for (index_type i=1; i<A_NDIM;   i++) szG[i]        = szL[i];
    for (index_type i=1; i<A_NDIM;   i++) starts[i]     = 0;
                                            
    szL[0]    = maxInt((nsize[0] - nghost_cells[0] - nghost_cells[1])/nout_avg,1);
    szG[0]    = szL[0]*nsubdomains;
    starts[0] = szL[0]*id_number;

    if ( strlen(name) > 0 ) {
      typedef ArrayNd<OTYPE,NDIM> OArrayNd;
      OArrayNd f(INDICIES(szG[0],szG[1],szG[2]));
      float navg=pow(nout_avg,OTYPE(NDIM));                                    //number of cells to average over

      /* for (int ix=0; ix<array_size[0]-nghost_cells[0]-nghost_cells[1]; ix+=nout_avg) */
      for (int ix=nghost_cells[0]; ix<array_size[0]-nghost_cells[1]; ix+=nout_avg)
	for (int iy=0; iy<array_size[1]; iy += nout_avg) 
	  for (int iz=0; iz<array_size[2]; iz += nout_avg) 
	    for (int iu=0; iu<array_size[3]; iu += nout_avg) 
	      for (int iv=0; iv<array_size[4]; iv += nout_avg) 
		for (int iw=0; iw<array_size[5]; iw += nout_avg) {
		  //output cell location
		  int ix3=(ix-nghost_cells[0])/nout_avg, iy3=iy/nout_avg, iz3=iz/nout_avg;
		  for (int ix2=ix; ix2<ix+nout_avg && 
			 ix2 < array_size[0]-nghost_cells[1]; ix2++) 
		    for (int iy2=iy; iy2<iy+nout_avg && 
			   iy2 < array_size[1]; iy2++) 
		      for (int iz2=iz; iz2<iz+nout_avg && 
			     iz2 < array_size[2]; iz2++) 
			for (int iu2=iu; iu2<iu+nout_avg && 
			       iu2 < array_size[3]; iu2++) 
			  for (int iv2=iv; iv2<iv+nout_avg && 
				 iv2 < array_size[4]; iv2++) 
			    for (int iw2=iw; iw2<iw+nout_avg && 
				   iw2 < array_size[2]; iw2++) {
			      f(INDICIES(ix3,iy3,iz3)) += data[ix2*array_step[0]+
			      				       iy2*array_step[1]+
			      				       iz2*array_step[2]+
			      				       iu2*array_step[3]+
			      				       iv2*array_step[4]+
			      				       iw2*array_step[5]];
			    } 
		  f(INDICIES(ix3,iy3,iz3)) /= navg;
		} 
      
      filespace = H5Screate_simple(A_NDIM, szG, NULL);
      memspace  = H5Screate_simple(A_NDIM, szL, NULL);
      dset_id   = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(filespace);

      // Select hyperslab in the file.
      filespace = H5Dget_space(dset_id);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, szL, NULL);
  
      // Tell HDF5 (and MPI I/O) to use collective writes
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      
      // Write the dataset
      H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, f.address());

      // Close HDF5 identifiers
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
      H5Pclose(plist_id);
    }

  }

  void bin_input_ghost(FILE* fname, int nghost_cells[2]=0){
    // Open a file, fname, if needed, and output an array in single precision 
    // binary, averaging over nout points.

    int array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix]=1;
    int array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix]=1;

    for (index_type i=0; i<A_NDIM; i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i+1];

    if ( fname != 0 ) {
      for (int ix=nghost_cells[0]; ix<array_size[0]-nghost_cells[1];ix ++)
	for (int iy=0; iy<array_size[1]; iy ++) 
	  for (int iz=0; iz<array_size[2]; iz ++) 
	    for (int iu=0; iu<array_size[3]; iu ++) 
	      for (int iv=0; iv<array_size[4]; iv ++) 
		for (int iw=0; iw<array_size[5]; iw ++) {
		  float f=0.;
		  fread(&f,sizeof(f),1,fname);
		  data[ix*array_step[0]+
		       iy*array_step[1]+
		       iz*array_step[2]+
		       iu*array_step[3]+
		       iv*array_step[4]+
		       iw*array_step[5]]=f;
		} 
    }
    
  }


};

#endif // ArrayNd_ranged_H 
