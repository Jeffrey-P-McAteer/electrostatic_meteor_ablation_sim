/* ArrayNd<Type,DIM>:  A simple N-Dimensional array class - 
   designed to be entirely inline. 

   Written by Meers Oppenheim starting in 1999
*/

#ifndef ArrayNd_H 
#define ArrayNd_H

#include <cstdlib>


#include <iostream>
#include <fstream>
#include <cstdio>
#include <csignal>
#include <math.h>

inline int maxInt(int a, int b) {
  return (a<b) ? b : a;
}

typedef int index_type;
using namespace std;

template <class TYPE,int A_NDIM>
class ArrayNd {

 protected:
  TYPE*	      data;	           // Array elements
  index_type  nsize[A_NDIM];       // Number of elements
  index_type  nelem[A_NDIM];       // products used for indexing
  bool ext_workspace;

 public:

  //  Constructor:
  ArrayNd() {
    data=0;
    //data=new TYPE[A_NDIM + 1]; // if nothing else touches this, data will point to a valid location of some sort w/ buffer space
    ext_workspace = false;
    for (int i=0; i < A_NDIM; i++) {
      nsize[i]=0;
      nelem[i]=0;
    }
  }


  void build_ArrayNd(const index_type n_in[], TYPE x=0, TYPE *workspace=0) 
    {
      for (index_type i=0; i<A_NDIM; i++) nsize[i]=n_in[i];
      check_dimensions();

      // nelem contains the number of elements between sucessive elements
      nelem[A_NDIM-1]=nsize[A_NDIM-1];
      for (int id = A_NDIM-2; id>=0; id--)  nelem[id]=nsize[id]*nelem[id+1];

      if (workspace==0) {
	data=new TYPE[nelem[0]];
	ext_workspace=false;
	} else {
	//It is the user's responsibility to insure that 
	// workspace is sufficiently large.
	data=workspace;
	ext_workspace=true;
      }
      for (index_type ip=0;ip<nelem[0];ip++) data[ip]=x;
    }

  ArrayNd(const index_type n_in[], TYPE x=0, TYPE *workspace=0) 
    {
      build_ArrayNd(n_in, x, workspace);
    }

  /* Obsolete:
  ArrayNd(const int n_in[], TYPE x=0, TYPE *workspace=0) 
    {
      index_type n_in2[A_NDIM];
      for (index_type i=0; i<A_NDIM; i++) n_in2[i]= n_in[i];

      build_ArrayNd(n_in2, x, workspace);
    }
  */


  ArrayNd(index_type n0_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 1, 
		  " More than 1 argument needed in ArrayNd<T,%d>", A_NDIM);
      
      index_type n_in2[]= {n0_in};

      build_ArrayNd(n_in2, 0, workspace);
      
    }

  ArrayNd(index_type n0_in, index_type n1_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 2, 
		  " More than 2 argument needed in ArrayNd<T,%d>",A_NDIM);
      index_type n_in2[A_NDIM]={n0_in,n1_in};
      build_ArrayNd(n_in2, 0, workspace);
    }

  ArrayNd(index_type n0_in, index_type n1_in, index_type n2_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 3, 
		  " More than 3 argument needed in ArrayNd<T,%d>",A_NDIM);
      index_type n_in2[A_NDIM]={n0_in,n1_in,n2_in};
      build_ArrayNd(n_in2, 0, workspace);
    }

  ArrayNd(index_type n0_in, index_type n1_in, index_type n2_in, 
	  index_type n3_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 4, 
		  " More than 4 argument needed in ArrayNd<T,%d>", A_NDIM);
      index_type n_in2[A_NDIM]={n0_in,n1_in,n2_in,n3_in};
      build_ArrayNd(n_in2, 0, workspace);
    }

  ArrayNd(index_type n0_in, index_type n1_in, index_type n2_in, 
	  index_type n3_in, index_type n4_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 5,
		  " More than 4 argument needed in ArrayNd<T,%d>", A_NDIM);
      index_type n_in2[A_NDIM]={n0_in,n1_in,n2_in,n3_in,n4_in};
      build_ArrayNd(n_in2, 0, workspace);
    }

  ArrayNd(index_type n0_in, index_type n1_in, index_type n2_in,
	  index_type n3_in, index_type n4_in, index_type n5_in, TYPE *workspace=0) 
    {
      check_error(A_NDIM > 6, 
		  " More than 5 argument needed in ArrayNd<T,%d>", A_NDIM);
      index_type n_in2[A_NDIM]={n0_in,n1_in,n2_in,n3_in,n4_in,n5_in};
      build_ArrayNd(n_in2, 0, workspace);
    }

  
  // Copy Constructor
  ArrayNd(const ArrayNd &in) {
      build_ArrayNd(in.nsize);
      for (index_type id=0;id<nelem[0];id++) data[id]=in.data[id];
  };

  
  //  Destructor
  ~ArrayNd(){
    if (!ext_workspace) {
      delete [] data;
      data=0;
    }
  }
  
  //  Access and setting -- Speed is everything
  inline TYPE& operator() (index_type n0){   
    // One can always access any dimension ArrayNd with a single index
#ifdef DEBUG
    check_error( n0>=nelem[0] && nelem[0]>0, "n0 = %d out of range (0,%d]", n0,nelem[0]);
#endif
    return data[n0];
  }

  /* Obsolete:
  //  This one is necessary to avoid ambiguity
  inline TYPE& operator() (int n0){   
    // One can always access any dimension ArrayNd with a single index
#ifdef DEBUG
    check_error( n0>=nelem[0], "n0 = %d out of range (0,%d]", n0,nelem[0]);
#endif
    return data[n0];
    } */


  inline TYPE& operator[] (index_type n0){   
    // One can always access any dimension ArrayNd with a single index
#ifdef DEBUG
    check_error( n0>=nelem[0] && nelem[0]>0, "n0 = %d out of range (0,%d]", n0,nelem[0]);
#endif
    return data[n0];
  }

  inline TYPE& operator() (index_type n0, index_type n1){   
    // 2-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 2,
		"operator(index_type n0, index_type n1) cannot be applied for A_NDIM<2");
    check_error( n0>=nsize[0] ,"n0 = %d out of range (0,%d]", n0, nsize[0]);
    check_error( n1>=nsize[1] ,"n1 = %d out of range (0,%d]", n1, nsize[1]);
#endif
    return data[n0*nelem[1]+n1];
  }

  inline TYPE& operator() (index_type n0, index_type n1, index_type n2){   
    // 3-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 3,
		"operator(index_type n0, index_type n1, index_type n2) cannot be applied for A_NDIM<3");
    check_error( n0>=nsize[0] ,"n0 = %d out of range (0,%d]", n0,nsize[0]);
    check_error( n1>=nsize[1] ,"n1 = %d out of range (0,%d]", n1,nsize[1]);
    check_error( n2>=nsize[2] ,"n2 = %d out of range (0,%d]", n2,nsize[2]);
#endif
    return data[n0*nelem[1]+n1*nelem[2]+n2];
  }

  inline TYPE& operator() (index_type n0, index_type n1, index_type n2, index_type n3){   
    // 4-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 4,
		"operator(index_type n0, index_type n1, index_type n2, index_type n3) cannot be applied for A_NDIM<4");
    check_error( n0>=nsize[0], "n0 = %d out of range (0,%d]", n0, nsize[0]);
    check_error( n1>=nsize[1] ,"n1 = %d out of range (0,%d]", n1, nsize[1]);
    check_error( n2>=nsize[2] ,"n2 = %d out of range (0,%d]", n2, nsize[2]);
    check_error( n3>=nsize[3] ,"n3 = %d out of range (0,%d]", n3, nsize[3]);
#endif
    return data[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3];
  }

  inline TYPE& operator() (index_type n0, index_type n1, index_type n2, index_type n3, index_type n4){   
    // 5-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 5,
		"operator(index_type n0, index_type n1, index_type n2, index_type n3, index_type n4) cannot be applied for A_NDIM<5");
    check_error( n0>=nsize[0] ,"n0 = %d out of range (0,%d]", n0, nsize[0]);
    check_error( n1>=nsize[1] ,"n1 = %d out of range (0,%d]", n1, nsize[1]);
    check_error( n2>=nsize[2] ,"n2 = %d out of range (0,%d]", n2, nsize[2]);
    check_error( n3>=nsize[3] ,"n3 = %d out of range (0,%d]", n3, nsize[3]);
    check_error( n4>=nsize[4] ,"n4 = %d out of range (0,%d]", n4, nsize[4]);
#endif
    return data[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3*nelem[4]+n4];
  }
  
  inline TYPE& operator() (index_type n0, index_type n1, index_type n2, index_type n3, index_type n4, index_type n5){   
    // 6-D array access:
#ifdef DEBUG
    check_error(A_NDIM < 6,
		"operator(index_type n0, index_type n1, index_type n2, index_type n3, index_type n4, index_type n5) cannot be applied for A_NDIM<6");
    check_error( n0>=nsize[0] ,"n0 = %d out of range (0,%d]", n0, nsize[0]);
    check_error( n1>=nsize[1] ,"n1 = %d out of range (0,%d]", n1, nsize[1]);
    check_error( n2>=nsize[2] ,"n2 = %d out of range (0,%d]", n2, nsize[2]);
    check_error( n3>=nsize[3] ,"n3 = %d out of range (0,%d]", n3, nsize[3]);
    check_error( n4>=nsize[4] ,"n4 = %d out of range (0,%d]", n4, nsize[4]);
    check_error( n5>=nsize[5] ,"n5 = %d out of range (0,%d]", n5, nsize[5]);
#endif
    return data[n0*nelem[1]+n1*nelem[2]+n2*nelem[3]+n3*nelem[4]+n4*nelem[5]+n5];
  }
  // Generalize access operator with and index argument
  inline TYPE& operator() (index_type ind[]) {
    index_type index=ind[A_NDIM-1];  // The last term need not be multiplied.
    for (index_type i=0; i<A_NDIM-1; i++) index += ind[i]*nelem[i+1];
    return operator[](index);
  }

  // To pass matricies to canned routines one sometimes needs the address:
  TYPE *address() {return data;}
  TYPE *address(int i) {return data+i;}
  TYPE *address(unsigned int i) {return data+i;}
  TYPE* address(int const  ind[]) { // For references in a vector
    index_type index=ind[A_NDIM-1];  
    for (index_type i=0; i<A_NDIM-1; i++) index += ind[i]*nelem[i+1];
    return data+index;
  }

  
  inline index_type length() const {return nelem[0];};
  inline index_type size() const {return nelem[0];};
  inline index_type size(index_type i) const {return nsize[i];};
  //  inline const int size(int i) const {return static_cast<int>(nsize[i]);};

  inline index_type xsize() const {return nsize[0];};
  inline index_type ysize() const {return nsize[1];};
  inline index_type zsize() const {return nsize[2];};

  // Starting index: while this is trivial for this class, it will make it more generally useful for the inhereted ranged class:

  inline index_type xstart() const {return 0;};
  inline index_type ystart() const {return 0;};
  inline index_type zstart() const {return 0;};

  //  Assignment
  inline const ArrayNd& operator= (TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i]=x;
    return *this;
  };
  
  const ArrayNd& operator= (const ArrayNd &in) {
    // If they're the same size copy the data
    //otherwise, reinitialize:
    index_type old_nelem=nelem[0];
    for (index_type i=0; i<A_NDIM; i++) {
      nsize[i]=in.nsize[i];
      nelem[i]=in.nelem[i];
    }
    if (old_nelem != in.nelem[0]) {
      delete [] data;
      if  (nelem[0] > 0) data=new TYPE[nelem[0]];
      else data=0;
    }
    for (index_type id=0;id<nelem[0];id++) data[id]=in.data[id];
    return *this;
  };

//  Arithmetic:
  const ArrayNd& operator+= (const ArrayNd &x) {
    check_error (nelem[0] != x.nelem[0],
		 "ArrayNd += operator not given same size input");
    for (index_type i=0;i<nelem[0];i++) data[i] += x.data[i];
    return *this;
  };
  
  const ArrayNd& operator-= (const ArrayNd &x) {
    check_error (nelem[0] != x.nelem[0],
	       "ArrayNd -= operator not given same size input");
    for (index_type i=0;i<nelem[0];i++) data[i] -= x.data[i];
    return *this;
  };
  
  const ArrayNd& operator*= (const ArrayNd &x) {
    check_error (nelem[0] != x.nelem[0],
	       "ArrayNd *= operator not given same size input");
    for (index_type i=0;i<nelem[0];i++) {
      data[i] *= x.data[i];
    }
    return *this;
  };

  const ArrayNd& operator/= (const ArrayNd &x) {
    check_error (nelem[0] != x.nelem[0],
		 "ArrayNd /= operator not given same size input");
    for (index_type i=0;i<nelem[0];i++) data[i] /= x.data[i];
    return *this;
  };

  ArrayNd& operator*= (const TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i] *= x;
    return *this;
  };

  ArrayNd& operator/= (const TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i] /= x;
    return *this;
  };

  ArrayNd& operator+= (const TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i] += x;
    return *this;
  };

  ArrayNd& operator-= (const TYPE x) {
    for (index_type i=0;i<nelem[0];i++) data[i] -= x;
    return *this;
  };


  //    ArrayNd operator+ (const ArrayNd &);
  //    ArrayNd operator- (const ArrayNd &);
  //    ArrayNd operator* (const ArrayNd &);
  //    ArrayNd operator/ (const ArrayNd &);

  // Take the values from a smaller array and incorporate them into this array:
  /*  void assign_values(ArrayNd &in, TYPE val_default) {
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
    }*/



  TYPE max() {
    TYPE value=data[0];
    for (index_type i=1;i<nelem[0];i++) if (data[i] > value) value=data[i];
    return value;
  }

  double absmax() {
    double value=abs(data[0]);
    for (index_type i=1;i<nelem[0];i++) if (abs(data[i]) > value) value=abs(data[i]);
    return value;
  }

  TYPE min() {
    TYPE value=data[0];
    for (index_type i=1;i<nelem[0];i++) if (data[i] < value) value=data[i];
    return value;
  }

  TYPE sum() {
    TYPE value=data[0];
    for (index_type i=1;i<nelem[0];i++) value += data[i];
    return value;
  }

  TYPE sumsq() {
    TYPE value=data[0]*data[0];
    for (index_type i=1;i<nelem[0];i++) value += data[i]*data[i];
    return value;
  }

 protected:

  void check_error(bool cond, const char *complaint) const{
    if (cond) {
      cout << "\nError: Problem in ArrayNd\n"
	   << "\t" << complaint << endl;
      cerr << "\nError: Problem in ArrayNd\n"
	   << "\t" << complaint << endl;

      raise(SIGSEGV);
      exit(1);
    }
  };
  void check_error(bool cond, const char *complaint, index_type i1) const{
    if (cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1);
      check_error(0==0,errmess);
    }
  };
  void check_error(bool cond, const char *complaint, index_type i1, index_type i2) const{
    if (cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2);
      check_error(0==0, errmess);
    }
  };
  void check_error(bool cond, const char *complaint, index_type i1, index_type i2, index_type i3) 
    const {
    if (cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2, i3);
      check_error(0==0, errmess);
    }
  };
  void check_error(bool cond, const char *complaint, 
		   index_type i1, index_type i2, index_type i3, index_type i4) const {
    if (cond) {
      char errmess[128];
      sprintf(errmess, complaint, i1, i2, i3, i4);
      check_error(0==0, errmess);
    }
  };
  
  void check_dimensions() const {
    for (index_type i=0; i<A_NDIM; i++) {
      check_error(nsize[i] <= 0,
		  "ArrayND<T,%d> with index %d equal to %d not permitted",
		  A_NDIM, i, nsize[i]);
    }
  }
public:  

  friend ostream &operator<<(ostream& ostr, const ArrayNd& x){

    ostr << "Array( \t";
    for (int nd=0; nd<A_NDIM; nd++) {
      ostr << x.nsize[nd];
      if (nd < A_NDIM-1) ostr << "\t"; else ostr << ")=" << endl;
    }
    
    for (int i=0; i<x.nelem[0]; i++) {
      ostr << i << " \t" ;
      for (int nd=0; nd<A_NDIM-1; nd++){
	index_type ic=(i%x.nelem[nd])/x.nelem[nd+1];
	ostr << ic << " \t" ;
      }
      if (A_NDIM > 1) 	ostr <<  i%x.nelem[A_NDIM-1] << " \t" ;
      ostr << x.data[i] << endl;
    }
    
    return ostr;
    
  };

  /* A simple output to a file: */
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
    int nx; in >> nx;
    int ny; in >> ny;
    //    cout << "About to build arrays...\n";
    if (in.peek() == ')') {
      index_type n_in2[2]={nx,ny};
      build_ArrayNd(n_in2, 0, workspace);
    } else { /* assume 3D */
      int nz; in >> nz;
      index_type n_in3[3]={nx,ny,nz};
      build_ArrayNd(n_in3, 0, workspace);
    }
    //    cout << "in pos: " << in.tellg() << "\n";
    in.ignore(2);
    //    cout << "in pos: " << in.tellg() << "\n";
    //    cout << "number of elements " << nelem[0] << "\n";
    for (int idata=0; idata<nelem[0]; idata++) {
      int count,ix,iy;
      in >> count;
      //      cout << count << "\t";
      //      cout << "in pos: " << in.tellg() << "\n";
      in >> ix; 
      //cout << ix << "\t";
      //cout << "in pos: " << in.tellg() << "\n";
      in  >> iy; 
      //cout << iy << "\t"; 
      //cout << "in pos: " << in.tellg() << "\n";
      in >> data[idata];
      //      cout << data[idata] << "\n";
      //cout << "in pos: " << in.tellg() << " idata " << idata << "\n";

    }
    in.close();
    return 0; 
  }; 

  void bin_output(FILE* fname, int nout_avg=1){
    // Open a file, fname, if needed, and output an array in single precision 
    // binary, averaging over nout points.

    int array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix]=1;
    int array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix]=1;

    for (index_type i=0; i<A_NDIM; i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i+1];

    if ( fname != 0 ) {
      for (int ix=0; ix<array_size[0]; ix += nout_avg)
	for (int iy=0; iy<array_size[1]; iy += nout_avg) 
	  for (int iz=0; iz<array_size[2]; iz += nout_avg) 
	    for (int iu=0; iu<array_size[3]; iu += nout_avg) 
	      for (int iv=0; iv<array_size[4]; iv += nout_avg) 
		for (int iw=0; iw<array_size[5]; iw += nout_avg) {
		  float f=0.;
		  int navg=0;
		  for (int ix2=ix; ix2<ix+nout_avg && 
			 ix2 < array_size[0]; ix2++) 
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

  void bin_input(FILE* fname){
    // Open a file, fname, if needed, and output an array in single precision 
    // binary, averaging over nout points.

    int array_size[6]; for(int ix=0;ix<6;ix++) array_size[ix]=1;
    int array_step[6]; for(int ix=0;ix<6;ix++) array_step[ix]=1;

    for (index_type i=0; i<A_NDIM; i++) array_size[i] = nsize[i];
    for (index_type i=0; i<A_NDIM-1; i++) array_step[i] = nelem[i];

    if ( fname != 0 ) {
      for (int ix=0; ix<array_size[0]; ix ++)
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

#endif // ArrayNd_H 
