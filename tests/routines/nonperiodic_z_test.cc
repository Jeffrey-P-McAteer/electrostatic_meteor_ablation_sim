#include <stdio.h>
#include "tridiag.h"

void print_vectors(ArrayNd<FTYPE,1> a,ArrayNd<FTYPE,1> b,ArrayNd<FTYPE,1> c,
                   ArrayNd<FTYPE,1> rhs, ArrayNd<FTYPE,1> lhs, int arrsize){
    for (int i=0; i<arrsize; i++){
      if (i == 0){
        printf("[ %1.2f %1.2f 0 .............. ] [ %1.2f ]   [ %f ]\n", b(i), c(i), lhs(i), rhs(i));
      } else if (i == 1){
        printf("[ %1.2f %1.2f %1.2f 0 ......... ] [ %1.2f ]   [ %f ]\n", a(i), b(i), c(i), lhs(i), rhs(i));
      } else if (i == arrsize-1){
        printf("[            ... 0 %1.2f %1.2f ] [ %1.2f ]   [ %f ]\n", a(i), b(i), lhs(i), rhs(i));
      } else if (i == arrsize-2){
        printf("[    ...  0  %1.2f  %1.2f %1.2f ] [ %1.2f ]   [ %f ]\n", a(i), b(i), c(i), lhs(i), rhs(i));
      } else {
        printf("[ ... 0 %1.2f %1.2f %1.2f 0 ... ] [ %1.2f ] = [ %f ]\n",a(i),b(i),c(i), lhs(i), rhs(i));
      }
    }
};

void print_vectors(std::complex<FTYPE> *a, int arrsize){
  for (int i=0; i<arrsize; i++){
    printf("[ %1.2f + i%1.2f ]\n", a[i].real(), a[i].imag());
  }
}

void print_vectors(std::complex<FTYPE> *a, std::complex<FTYPE> *b, int arrsize){
  for (int i=0; i<arrsize; i++){
    printf("[ %1.2f + i%1.2f ]  [ %1.2f + i%1.2f ]\n", 
           a[i].real(), a[i].imag(), b[i].real(), b[i].imag());
  }
}

int main(int argc, char* argv[])
{
  printf("Inside nonperiodic_z_test.cc\n");
  // Test tridiag_solver
  int arrsize = 7;
  
  int testtype = 4;

  if (testtype == 0){
    ArrayNd<FTYPE,1> a = ArrayNd<FTYPE,1>(arrsize);
    ArrayNd<FTYPE,1> b = ArrayNd<FTYPE,1>(arrsize);
    ArrayNd<FTYPE,1> c = ArrayNd<FTYPE,1>(arrsize);
    ArrayNd<FTYPE,1> rhs1 = ArrayNd<FTYPE,1>(arrsize);
    ArrayNd<FTYPE,1> lhs = ArrayNd<FTYPE,1>(arrsize);
    for (int i=0; i<arrsize; i++){
      a(i) = 1;
      b(i) = -2;
      c(i) = 1;
      rhs1(i) = i;
    }
    print_vectors(a,b,c,rhs1,lhs,arrsize);
    
    tridiag_solver(a,b,c,rhs1,lhs);
    printf("After tridiag_solver...\n");
    print_vectors(a,b,c,rhs1,lhs,arrsize);
  }  

  // Test dirichlet-dirichlet
  if (testtype ==1){
    FTYPE boundary_condition_low = 2.0;
    FTYPE boundary_condition_high = 3.0;
    FTYPE b = -2.0;
    std::complex<FTYPE> rhs[arrsize];
    std::complex<FTYPE> lhs[arrsize];
    int n_unknown = arrsize;
    
    for (int i=0; i<arrsize; i++){
      rhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
      lhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
    }
    // print the complex vectors
    for(int i=0; i<arrsize; i++){
      printf("rhs[%d] = %1.3f+i%1.3f, lhs[%d] = %1.3f+i%1.3f\n",i,
             rhs[i].real(), rhs[i].imag(), i, lhs[i].real(), lhs[i].imag());
    }
    
    print_vectors(lhs, rhs, arrsize);
    
    
    printf("Using dirichlet solver\n");
    tridiag_solver_complex_dirichlet_dirichlet(boundary_condition_low,
                                               boundary_condition_high,
                                               b, rhs, lhs, arrsize);
  /*
    lhs[0] = 0.0;
    lhs[arrsize-1] = 0.0;
    printf("Using tridiagsolvercomplex\n");
    void tridiag_solver_complex(FTYPE b,
    std::complex<FTYPE> *rhs,
    std::complex<FTYPE> *lhs,
    int n_unknown);
    tridiag_solver_complex(b, &rhs[1], &lhs[1], arrsize-2);
  */
    
    printf("After dirichlet\n");
    print_vectors(lhs, rhs, arrsize);
  }
  
  // Test Dirichlet-Nuemann
  if (testtype == 2){
    printf("Testing Dirichlet-Neumann\n");
    FTYPE boundary_condition_low = 2.0;
    FTYPE boundary_condition_high = 3.0;
    FTYPE b = -2.0;
    std::complex<FTYPE> rhs[arrsize];
    std::complex<FTYPE> lhs[arrsize];
    int n_unknown = arrsize;
    
    for (int i=0; i<arrsize; i++){
      rhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
      lhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
    }
    // print the complex vectors
    for(int i=0; i<arrsize; i++){
      printf("rhs[%d] = %1.3f+i%1.3f, lhs[%d] = %1.3f+i%1.3f\n",i,
             rhs[i].real(), rhs[i].imag(), i, lhs[i].real(), lhs[i].imag());
    }
    
    print_vectors(lhs, rhs, arrsize);
    
    printf("Using dirichlet-neumann solver\n");
    tridiag_solver_complex_dirichlet_neumann(boundary_condition_low,
                                             boundary_condition_high,
                                             b, rhs, lhs, arrsize);
    printf("After dirichlet-neumann\n");
    print_vectors(lhs, rhs, arrsize);
  }

  if (testtype == 3){
    printf("Testing Neumann-Dirichlet\n");
    FTYPE boundary_condition_low = 2.0;
    FTYPE boundary_condition_high = 3.0;
    FTYPE b = -2.0;
    std::complex<FTYPE> rhs[arrsize];
    std::complex<FTYPE> lhs[arrsize];
    int n_unknown = arrsize;
    
    for (int i=0; i<arrsize; i++){
      rhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
      lhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
    }
    // print the complex vectors
    for(int i=0; i<arrsize; i++){
      printf("rhs[%d] = %1.3f+i%1.3f, lhs[%d] = %1.3f+i%1.3f\n",i,
             rhs[i].real(), rhs[i].imag(), i, lhs[i].real(), lhs[i].imag());
    }
    
    print_vectors(lhs, rhs, arrsize);
    
    printf("Using neumann-dirichlet solver\n");
    tridiag_solver_complex_neumann_dirichlet(boundary_condition_low,
                                             boundary_condition_high,
                                             b, rhs, lhs, arrsize);
    printf("After neumann dirichlet\n");
    print_vectors(lhs, rhs, arrsize);
  }

  if (testtype == 4){
    printf("Testing Neumann-Neumann\n");
    FTYPE boundary_condition_low = 2.0;
    FTYPE boundary_condition_high = 3.0;
    FTYPE b = -2.0+0.5;
    //FTYPE b = -2.0;
    boundary_condition_low = 0.0;
    boundary_condition_high = 0.0;
    std::complex<FTYPE> rhs[arrsize];
    std::complex<FTYPE> lhs[arrsize];
    int n_unknown = arrsize;
    
    for (int i=0; i<arrsize; i++){
      //rhs[i] = std::complex<FTYPE>((FTYPE) i-(arrsize-1.0)/2.0, (FTYPE)-i+(arrsize-1.0)/2.0);
      //lhs[i] = std::complex<FTYPE>((FTYPE) i-(arrsize-1.0)/2.0, (FTYPE)-i+(arrsize-1.0)/2.0);
      rhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
      lhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
    }
    // print the complex vectors
    for(int i=0; i<arrsize; i++){
      printf("rhs[%d] = %1.3f+i%1.3f, lhs[%d] = %1.3f+i%1.3f\n",i,
             rhs[i].real(), rhs[i].imag(), i, lhs[i].real(), lhs[i].imag());
    }    
    print_vectors(lhs, rhs, arrsize);
    
    printf("Using neumann-neumann solver\n");
    tridiag_solver_complex_neumann_neumann(boundary_condition_low,
                                           boundary_condition_high,
                                           b, rhs, lhs, arrsize);
    printf("After neumann-neumann\n");
    print_vectors(lhs, rhs, arrsize);
  }  

  if (testtype == 5){
    printf("Testing Open-Open\n");
    FTYPE m=0;
    FTYPE b=-2.0*(1+2.0*pow(sin(3.14159*m/arrsize),2));
    FTYPE r = -b/2.0 + sqrt(b*b/4.0-1);
    FTYPE rinv = 1.0/r;
    std::complex<FTYPE> rhs[arrsize];
    std::complex<FTYPE> lhs[arrsize];
    int n_unknown = arrsize;
    
    for (int i=0; i<arrsize; i++){
      //rhs[i] = std::complex<FTYPE>((FTYPE) i-(arrsize-1.0)/2.0, (FTYPE)-i+(arrsize-1.0)/2.0);
      //lhs[i] = std::complex<FTYPE>((FTYPE) i-(arrsize-1.0)/2.0, (FTYPE)-i+(arrsize-1.0)/2.0);
      rhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
      lhs[i] = std::complex<FTYPE>((FTYPE) i, (FTYPE)i+1);
    }
    // print the complex vectors
    for(int i=0; i<arrsize; i++){
      printf("rhs[%d] = %1.3f+i%1.3f, lhs[%d] = %1.3f+i%1.3f\n",i,
             rhs[i].real(), rhs[i].imag(), i, lhs[i].real(), lhs[i].imag());
    }    

    print_vectors(lhs,rhs,arrsize);
    printf("Using open-open solver\n");
    tridiag_solver_complex_open_open(b, rinv, rhs, lhs, arrsize);
    printf("After open-open\n");
    print_vectors(lhs,rhs,arrsize);
  }
  return 0;
}
