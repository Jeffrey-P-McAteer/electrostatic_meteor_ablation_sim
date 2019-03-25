// Array2d<Type>:  A simple 3D array class - designed to be entirely inline.

#ifndef Array2d_H 
#define Array2d_H

#include <iostream>
#include <stdio>
#include <signal>

template <class TYPE>
class Array2d {

 protected:
  TYPE*	data;	           // Array elements
  int	nx, ny    ;        // Number of elements
  int   nxny;              // products used for indexing

  /*
  void check_error(bool, const char *) const;
  void check_error(bool, const char *, int) const;
  void check_error(bool, const char *, int, int) const;
  void check_error(bool, const char *, int, int, int) const;
  void check_error(bool, const char *, int, int, int, int) const;
  */
 public:

  //  Constructor:

  Array2d(int nx_in=1, int ny_in=1, TYPE x=0.): 
    nx(nx_in), ny(ny_in), nxny(nx_in*ny_in)
    {
      data=new TYPE[nx*ny];
      for (int i=0;i<nx*ny;i++) data[i]=x;
    };

  Array2d(int n_in[3]): 
    nx(n_in[0]), ny(n_in[1]),  nxny(n_in[0]*n_in[1])
    {
      data=new TYPE[nx*ny];
      for (int i=0;i<nx*ny;i++) data[i]=0;
    };
  
  // Copy Constructor
  Array2d(const Array2d &in) {
    nx=in.nx;
    ny=in.ny;
    nxny=nx*ny;
    data=new TYPE[nxny];
    for (int i=0;i<nx*ny;i++) data[i]=in.data[i];
  };

  
  //  Destructor
  ~Array2d(){
    delete [] data;
  }
  
  //  Access and setting
  TYPE elem(int ix, int iy) { 
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
#endif
    return data[iy + ix*ny];
  }

  void elem(int ix, int iy, TYPE elem){
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
#endif
    data[iy + ix*ny] = elem;
  }
  
  TYPE& operator ()(int ix, int iy) {
#ifdef DEBUG
    check_error(ix >= 0 && ix < nx,"ix = %d out of range [%d,%d]",ix,0,nx);
    check_error(iy >= 0 && iy < ny,"iy = %d out of range [%d,%d]",iy,0,ny);
#endif
    return data[iy + ix*ny];
  }

TYPE& operator ()(int i) {
#ifdef DEBUG
    check_error(i >= 0 && i < nxny,"i = %d out of range [%d,%d]",i,0,nxny);
#endif
    return data[i];
  }

  // To pass matricies to canned routines one sometimes needs the address:
  TYPE *address() {return data;}
  
  int size() {return nxny;};
  int xsize() {return nx;};
  int ysize() {return ny;};
  int zsize() {return 1;};

  //  Assignment
  const Array2d& operator=(TYPE x) {
    for (int i=0;i<nxny;i++) data[i]=x;
    return *this;
  };
  
  const Array2d& operator=(const Array2d &in) {
    // If they're the same size copy the data
    //otherwise, reinitialize:
    if (nx*ny != in.nxny) {
      delete [] data;
      nx=in.nx;
      ny=in.ny;
      nxny=nx*ny;
      data=new TYPE[nxny];
    } 
    for (int i=0;i<nx*ny;i++) data[i]=in.data[i];
    return *this;
  };



//  Arithmetic
  //    Array2d operator+ (const Array2d &);
  //    Array2d operator- (const Array2d &);
  //    Array2d operator* (const Array2d &);
  //    Array2d operator/ (const Array2d &);
  const Array2d& operator+= (const Array2d &x) {
    check_error (nxny == x.nxny,
	       "Array2d += operator not given same size input (%d vs. %d)",
		 nxny, x.nxny);
    for (int i=0;i<nx*ny;i++) data[i] += x.data[i];
    return *this;
  };
  
  const Array2d& operator-= (const Array2d &x) {
    check_error (nxny == x.nxny,
	       "Array2d -= operator not given same size input (%d vs. %d)",
		 nxny, x.nxny);
    for (int i=0;i<nxny;i++) data[i] -= x.data[i];
    return *this;
  };
  
  const Array2d& operator*= (const Array2d &x) {
    check_error (nxny == x.nxny,
	       "Array2d *= operator not given same size input(%d vs. %d)",
		 nxny, x.nxny);
    for (int i=0;i<nxny;i++) data[i] *= x.data[i];
    return *this;
  };

  const Array2d& operator/= (const Array2d &x) {
    check_error (nxny == x.nxny,
		 "Array2d /= operator not given same size input (%d vs. %d)",
		 nxny, x.nxny);
    for (int i=0;i<nxny;i++) data[i] /= x.data[i];
    return *this;
  };

  Array2d& operator*= (const TYPE x) {
    for (int i=0;i<nxny;i++) data[i] *= x;
    return *this;
  };

  Array2d& operator/= (const TYPE x) {
    for (int i=0;i<nxny;i++) data[i] /= x;
    return *this;
  };

  Array2d& operator+= (const TYPE x) {
    for (int i=0;i<nxny;i++) data[i] += x;
    return *this;
  };

  Array2d& operator-= (const TYPE x) {
    for (int i=0;i<nxny;i++) data[i] -= x;
    return *this;
  };

 protected:

  void check_error(bool cond, const char *complaint) const{
    if (!cond) {
      cout << "Error: Problem in Array2d\n"
	   << "\t" << complaint << endl;
      cerr << "Error: Problem in Array2d\n"
	   << "\t" << complaint << endl;

      raise(SIGINT);
      terminate(1,"Error in Array2d");
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

  friend ostream& operator<<(ostream& ostr, const Array2d& x){
    for (int ix=0; ix<x.nx; ix++) 
      for (int iy=0; iy<x.ny; iy++) 
	  ostr << ix << " \t" << iy << " \t" 
	       << x.data[iy + ix*x.ny] << endl;
    return ostr;
    
  };
  
    
      
};

#endif // Array2d_H 
