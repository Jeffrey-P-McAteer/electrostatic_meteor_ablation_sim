// Array3d<Type>:  A simple 3D array class - designed to be entirely inline.

#ifndef Array3d_H 
#define Array3d_H

#include <iostream>
#include <stdio>
#include <signal>

template <class TYPE>
class Array3d {

 protected:
  TYPE*	data;	           // Array elements
  int	nx, ny, nz;        // Number of elements
  int   nynz, nxnynz;      // products used for indexing

  /*
  void check_error(bool, const char *) const;
  void check_error(bool, const char *, int) const;
  void check_error(bool, const char *, int, int) const;
  void check_error(bool, const char *, int, int, int) const;
  void check_error(bool, const char *, int, int, int, int) const;
  */
 public:

  //  Constructor:

  Array3d(int nx_in=1, int ny_in=1, int nz_in=1, TYPE x=0.): 
    nx(nx_in), ny(ny_in), nz(nz_in), 
    nynz(ny_in*nz_in), nxnynz(nx_in*ny_in*nz_in)
    {
      data=new TYPE[nx*ny*nz];
      for (int i=0;i<nx*ny*nz;i++) data[i]=x;
    };

  Array3d(int n_in[3]): 
    nx(n_in[0]), ny(n_in[1]), nz(n_in[2]), 
    nynz(n_in[1]*n_in[2]), nxnynz(n_in[0]*n_in[1]*n_in[2])
    {
      data=new TYPE[nx*ny*nz];
      for (int i=0;i<nx*ny*nz;i++) data[i]=0;
    };
  
  // Copy Constructor
  Array3d(const Array3d &in) {
    nx=in.nx;
    ny=in.ny;
    nz=in.nz;
    nynz=ny*nz;
    nxnynz=nx*ny*nz;
    data=new TYPE[nxnynz];
    for (int i=0;i<nx*ny*nz;i++) data[i]=in.data[i];
  };

  
  //  Destructor
  ~Array3d(){
    delete [] data;
  }
  
  //  Access and setting
  TYPE elem(int ix, int iy, int iz) { 
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
    check_error(iz >= 0 && iz < nz,"iz = %d out of range [%d,%d]",iz,0,nz);
#endif
    return data[iz + iy*nz + ix*nynz];
  }

  void elem(int ix, int iy, int iz, TYPE elem){
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
    check_error(iz >= 0 && iz < nz,"iz = %d out of range [%d,%d]",iz,0,nz);
#endif
    data[iz + iy*nz + ix*nynz] = elem;
  }
  
  TYPE& operator ()(int ix, int iy, int iz) {
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
    check_error(iz >= 0 && iz < nz,"iz = %d out of range [%d,%d]",iz,0,nz);
#endif
    return data[iz + iy*nz + ix*nynz];
  }

TYPE& operator ()(int i) {
#ifdef DEBUG
    check_error(i >= 0 && i < nxnynz,"i = %d out of range [%d,%d]",i,0,nxnynz);
#endif
    return data[i];
  }

  // To pass matricies to canned routines one sometimes needs the address:
  TYPE *address() {return data;}
  
  int size() {return nxnynz;};
  int xsize() {return nx;};
  int ysize() {return ny;};
  int zsize() {return nz;};

  //  Assignment
  const Array3d& operator=(TYPE x) {
    for (int i=0;i<nxnynz;i++) data[i]=x;
    return *this;
  };
  
  const Array3d& operator=(const Array3d &in) {
    // If they're the same size copy the data
    //otherwise, reinitialize:
    if (nx*ny*nz != in.nxnynz) {
      delete [] data;
      nx=in.nx;
      ny=in.ny;
      nz=in.nz;
      nynz=ny*nz;
      nxnynz=nx*ny*nz;
      data=new TYPE[nxnynz];
    } 
    //    check_error (nx*ny*nz == in.nxnynz,
    //    	       "Array3d = assignment not given matching sizes");
    for (int i=0;i<nx*ny*nz;i++) data[i]=in.data[i];
    return *this;
  };



//  Arithmetic
  //    Array3d operator+ (const Array3d &);
  //    Array3d operator- (const Array3d &);
  //    Array3d operator* (const Array3d &);
  //    Array3d operator/ (const Array3d &);
  const Array3d& operator+= (const Array3d &x) {
    check_error (nxnynz == x.nxnynz,
	       "Array3d += operator not given same size input (%d vs. %d)",
		 nxnynz, x.nxnynz);
    for (int i=0;i<nx*ny*nz;i++) data[i] += x.data[i];
    return *this;
  };
  
  const Array3d& operator-= (const Array3d &x) {
    check_error (nxnynz == x.nxnynz,
	       "Array3d -= operator not given same size input (%d vs. %d)",
		 nxnynz, x.nxnynz);
    for (int i=0;i<nx*ny*nz;i++) data[i] -= x.data[i];
    return *this;
  };
  
  const Array3d& operator*= (const Array3d &x) {
    check_error (nxnynz == x.nxnynz,
	       "Array3d *= operator not given same size input(%d vs. %d)",
		 nxnynz, x.nxnynz);
    for (int i=0;i<nx*ny*nz;i++) data[i] *= x.data[i];
    return *this;
  };

  const Array3d& operator/= (const Array3d &x) {
    check_error (nxnynz == x.nxnynz,
		 "Array3d /= operator not given same size input (%d vs. %d)",
		 nxnynz, x.nxnynz);
    for (int i=0;i<nx*ny*nz;i++) data[i] /= x.data[i];
    return *this;
  };

  Array3d& operator*= (const TYPE x) {
    for (int i=0;i<nx*ny*nz;i++) data[i] *= x;
    return *this;
  };

  Array3d& operator/= (const TYPE x) {
    for (int i=0;i<nx*ny*nz;i++) data[i] /= x;
    return *this;
  };

  Array3d& operator+= (const TYPE x) {
    for (int i=0;i<nx*ny*nz;i++) data[i] += x;
    return *this;
  };

  Array3d& operator-= (const TYPE x) {
    for (int i=0;i<nx*ny*nz;i++) data[i] -= x;
    return *this;
  };

 protected:

  void check_error(bool cond, const char *complaint) const{
    if (!cond) {
      cout << "Error: Problem in Array3d\n"
	   << "\t" << complaint << endl;
      cerr << "Error: Problem in Array3d\n"
	   << "\t" << complaint << endl;

      raise(SIGINT);
      terminate(1,"Error in Array3d");
    }
  };
  void check_error(bool cond, const char *complaint, int i1) const{
    if (!cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1);
      check_error(1==0,errmess);
    }
  };
  void check_error(bool cond, const char *complaint, int i1, int i2) const{
    if (!cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2);
      check_error(1==0, errmess);
    }
  };
  void check_error(bool cond, const char *complaint, int i1, int i2, int i3) 
    const {
    if (!cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2, i3);
      check_error(1==0, errmess);
    }
  };
  void check_error(bool cond, const char *complaint, 
		   int i1, int i2, int i3, int i4) const {
    if (!cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2, i3, i4);
      check_error(1==0, errmess);
    }
  };

public:  

  friend ostream& operator<<(ostream& ostr, const Array3d& x){
    for (int ix=0; ix<x.nx; ix++) 
      for (int iy=0; iy<x.ny; iy++) 
	for (int iz=0; iz<x.nz; iz++) 
	  ostr << ix << " \t" << iy << " \t" << iz << " \t"
	       << x.data[iz + iy*x.nz + ix*x.nynz] << endl;
    return ostr;
    
  };
  
    
      
};

#endif // Array3d_H 
