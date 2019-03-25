

#include "eppic-types.h"
#include "tridiag.h"

/* 
   Algorithm taken from Press et al., Numerical Recipies.
   Solves for x in Ax=b, where A is a tridiagonal matrix.
   INPUTS:
   a -- ArrayNd<FTYPE,1> -- lower stripe of tridiag matirx
   b -- ArrayNd<FTYPE,1> -- diagonal, middle stripe of tridiag matrix
   c -- ArrayNd<FTYPE,1> -- upper stripe of tridiag matrix
   rhs -- ArrayNd<FTYPE,1> -- corresponds to 'b' in Ax=b 
   lhs -- ArrayNd<FTYPE,1> -- solution and corresponds to 'x' in Ax=b

   Note: all inputs must be predefined and have the same dimensions. There 
   is no test done to verify this.
*/

void tridiag_solver(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
		    ArrayNd<FTYPE,1> &c,
		    ArrayNd<FTYPE,1>&rhs, ArrayNd<FTYPE,1> &lhs) 
{
  int n_unknown = rhs.size(0);  
  int j;
  FTYPE bet=b[0]; 
  FTYPEAVec gam=FTYPEAVec(n_unknown)=0;
  if (b[0] == 0.0)
    throw("Error 1 in tridiag_solver");

  lhs[0] = rhs[0]/bet;

  /* decomposition and forward substitution */
  for (j=1;j<n_unknown;j++) {
    gam[j] = c[j-1]/bet;
    bet=b[j]-a[j]*gam[j];
#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
#endif    
    lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
  }
  /* backsubstitution */
  for (j=(n_unknown-2);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];  
}


void tridiag_solver_complex_neumann_dirichlet(FTYPE boundary_condition_low,
                                              FTYPE boundary_condition_high,
                                              FTYPE b,
                                              std::complex<FTYPE> *rhs,
                                              std::complex<FTYPE> *lhs,
                                              int n_unknown)
{
  int j;
  std::complex<FTYPE> bet=b; 
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");
  
  // Set the first element to the low boundary condition (neumann)
  lhs[0] = (rhs[0]+2.0*boundary_condition_low)/bet;
  // Deal with the second element and c[0]=2 instead of c[0]=1
  gam[1] = 2.0/bet;
  bet = b-gam[1];
  lhs[1] = (rhs[1]-lhs[0])/bet;
  // decomposition and forward substitution 
  for (j=2;j<n_unknown-2;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
#endif   
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Take care of cell before dirichlet boundary
  gam[n_unknown-2] = 1.0/bet;
  bet = b-gam[n_unknown-2];
  lhs[n_unknown-2] = ((rhs[n_unknown-2]-boundary_condition_high)-lhs[n_unknown-3])/bet;
  // Take care of dirichlet boundary
  gam[n_unknown-1] = 1.0/bet;
  bet = b-gam[n_unknown-1];
  // Fill in the last element
  lhs[n_unknown-1] = boundary_condition_high;
  
  // backsubstitution (normally starts at j=n_unknown-2 but now it starts at 
  // j=n_unknown-3 because lhs[n_unknown-1] is constant
  for (j=(n_unknown-3);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];
}


void tridiag_solver_complex_dirichlet_neumann(FTYPE boundary_condition_low,
                                              FTYPE boundary_condition_high,
                                              FTYPE b,
                                              std::complex<FTYPE> *rhs,
                                              std::complex<FTYPE> *lhs,
                                              int n_unknown)
{
  int j;
  std::complex<FTYPE> bet=b; 
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  gam[1] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");
  
  // Set the first element to the low boundary condition (dirichlet)
  lhs[0] = boundary_condition_low;
  lhs[1] = (rhs[1]-lhs[0])/bet;
  // decomposition and forward substitution 
  for (j=2;j<n_unknown-1;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
#endif   
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Deal with the last element
  gam[n_unknown-1] = 1.0/bet;
  bet = b-(2.0*gam[n_unknown-1]);
  lhs[n_unknown-1] = ((rhs[n_unknown-1]-2.0*boundary_condition_high)-2.0*lhs[j-1])/bet;
  
  // backsubstitution 
  for (j=(n_unknown-2);j>=1;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];

/* Old method
  int j;
  //std::complex<FTYPE> bet=b[0]; 
  std::complex<FTYPE> bet=b; 
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");
  
  // Set the first element to the low boundary condition
  lhs[0] = boundary_condition_low;
  lhs[1] = (rhs[1]-lhs[0])/bet;
  // decomposition and forward substitution 
  for (j=2;j<n_unknown-2;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
    //#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
    //#endif   
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Deal with the second to last element
  gam[n_unknown-2] = 1.0/bet;
  bet = (b+4.0/3.0)-(gam[n_unknown-2]*2.0/3.0);
  lhs[n_unknown-2] = (rhs[n_unknown-2]-boundary_condition_high*2.0/3.0-lhs[n_unknown-3])/bet;
  
  // backsubstitution 
  for (j=(n_unknown-3);j>=1;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];

  // Solve lhs[n_unknown-1] enforcing Neumann condition
  lhs[n_unknown-1] = (2.0*boundary_condition_high+4.0*lhs[n_unknown-2]-lhs[n_unknown-3])/3.0;
*/
}


void tridiag_solver_complex_neumann_neumann(FTYPE boundary_condition_low,
                                            FTYPE boundary_condition_high,
                                            FTYPE b,
                                            std::complex<FTYPE> *rhs,
                                            std::complex<FTYPE> *lhs,
                                            int n_unknown)
{
  // Note that neumann-neumann has a singular matrix.  To fix this 
  // problem, we set the last element=0

  // make sure that b is not 2, if it is then we will have a singular matrix
  // if there is a singular matrix, eliminate the first row and set solution to 0

  int j;
  std::complex<FTYPE> bet=b; 
  std::complex<FTYPE> gam[n_unknown];
  // Test to make sure b is not -2 (singular)
  if (b == -2.0){
    gam[0] = 0.0;
    if (b == 0.0)
      throw("Error 1 in tridiag_solver");
    
    // Set the first element to the low boundary condition (neumann)
    lhs[0] = (rhs[0]+2.0*boundary_condition_low)/bet;
    // Deal with the second element and c[0]=2 instead of c[0]=1
    gam[1] = 2.0/bet;
    bet = b-gam[1];
    lhs[1] = (rhs[1]-lhs[0])/bet;
    // decomposition and forward substitution 
    for (j=2;j<n_unknown-2;j++) {
      //gam[j] = c[j-1]/bet;
      gam[j] = 1.0/bet;
      //bet=b[j]-a[j]*gam[j];
      bet=b-gam[j];
#if DEBUG
      if (bet==0.0)
        throw("Error 2 in tridiag_solver");
#endif   
      //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
      lhs[j]=(rhs[j]-lhs[j-1])/bet;
    }
    
    // Deal with the adjusted final (really 2nd to last) equation
    // 2nd to last because we throw away the redundant final equation
    gam[n_unknown-2] = 1.0/bet;
    bet = (b-2.0) - gam[n_unknown-2];
    lhs[n_unknown-2] = ((rhs[n_unknown-2]-rhs[n_unknown-1]+2.0*boundary_condition_high)
                        -lhs[n_unknown-3])/bet;
    // Set last value equal to 0
    lhs[n_unknown-1] = 0;
    
    // backsubstitution
    for (j=(n_unknown-3);j>=0;j--) 
      lhs[j] -= gam[j+1]*lhs[j+1];
  } else {
    // The matrix is not singular, so use normal algorithm
    gam[0] = 0.0;
    if (b == 0.0)
      throw("Error 1 in tridiag_solver");
    
    // Set the first element to the low boundary condition (neumann)
    // rhs[0] = rhs[0] + 2*bchigh
    lhs[0] = (rhs[0]+2.0*boundary_condition_low)/bet;
    // Deal with the second element and c[0]=2 instead of c[0]=1
    gam[1] = 2.0/bet;
    bet = b-gam[1];
    lhs[1] = (rhs[1]-lhs[0])/bet;
    // decomposition and forward substitution 
    for (j=2;j<n_unknown-1;j++) {
      //gam[j] = c[j-1]/bet;
      gam[j] = 1.0/bet;
      //bet=b[j]-a[j]*gam[j];
      bet=b-gam[j];
#if DEBUG
      if (bet==0.0)
        throw("Error 2 in tridiag_solver");
#endif   
      //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
      lhs[j]=(rhs[j]-lhs[j-1])/bet;
    }
    // Deal with the adjusted final equation
    gam[n_unknown-1] = 1.0/bet;
    bet = b - (2.0*gam[n_unknown-1]);
    lhs[n_unknown-1] = ((rhs[n_unknown-1]-2.0*boundary_condition_high)
                        -2.0*lhs[n_unknown-2])/bet;
    
    // backsubstitution
    for (j=(n_unknown-2); j>=0; j--){
      lhs[j] -= gam[j+1]*lhs[j+1];
    }
  }
}


// ASSUMES a and c = 1, b is a constant except at boundaries where it is a function of r
void tridiag_solver_complex_open_open(FTYPE b, FTYPE rinv,
                                      std::complex<FTYPE> *rhs,
                                      std::complex<FTYPE> *lhs,
                                      int n_unknown)
{
  // Test to make sure b is not -2 (singluar)
  if (b == -2.0){
    // printf("Singular!\n");
    // Solve for phi_-1 and phi_0 using equation 19 on page 322 of Birdsall Langdon
    std::complex<FTYPE> phi_im1=0.0;
    std::complex<FTYPE> phi_i=0.0;
    // Declare a variable to store the rhs previous value in case rhs,lhs are the same pointer
    std::complex<FTYPE> rhs_m1 = 0.0;
    std::complex<FTYPE> rhs_store = 0.0;
/* // Calculate phi_-1 and phi_N
    std::complex<FTYPE> gam[n_unknown];
    std::complex<FTYPE> bet=b;
    for (int i=0; i<n_unknown; i++){
      phi_im1 += abs((FTYPE) i + 1.0)*rhs[i];
      phi_i   += abs((FTYPE) i - n_unknown)*rhs[i];
    }
    phi_im1 /= 2.0;
    phi_i /= 2.0;
    int j;
    lhs[0] = (rhs[0]-phi_im1)/bet;
    for (j=1; j<n_unknown-1;j++){
      //gam[j] = c[j-1]/bet;
      gam[j] = 1.0/bet;
      //bet=b[j]-a[j]*gam[j];
      bet=b-gam[j];
      //#if DEBUG
      if (bet==0.0) 
        throw("Error 2 in tridiag_solver");
      //#endif 
      //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
      lhs[j]=(rhs[j]-lhs[j-1])/bet;
    }
    // Take care of last index
    gam[n_unknown-1] = 1.0/bet;
    bet = b-gam[n_unknown-1];
    lhs[n_unknown-1] = ((rhs[n_unknown-1]-phi_i)-lhs[n_unknown-2])/bet;
    // Backsubstitution
    for (j=(n_unknown-2);j>=0;j--){
      lhs[j] -= gam[j+1]*lhs[j+1];
    }

*/
    // Calculate phi_-1 and phi_0
    for (int i=0; i<n_unknown; i++){
      phi_im1 += abs((FTYPE) i + 1.0)*rhs[i];
      phi_i   += abs((FTYPE) i)*rhs[i];
    }
    phi_im1 /= 2.0;
    phi_i /= 2.0;
    // Sweep through all z indicies and solve for phi using eq 2 from Birdsall pg 319
    rhs_store = rhs[0];
    lhs[0] = phi_i;
    rhs_m1 = rhs_store;
    rhs_store = rhs[1];
    lhs[1] = rhs_m1 - b*lhs[0] - phi_im1;
    for (int i=2; i<n_unknown; i++){
      rhs_m1 = rhs_store;
      rhs_store = rhs[i];
      lhs[i] = rhs_m1 - b*lhs[i-1] - lhs[i-2];
    }

  } else {
    // The matrix is not singular, so use algorithm described in Birdsall pg 321
    std::complex<FTYPE> psi[n_unknown];
    psi[n_unknown-1] = rinv*rhs[n_unknown];
    for (int i=n_unknown-2; i>=0; i--){
      psi[i] = rinv*(rhs[i]+psi[i+1]);
    }
    
    lhs[0] = psi[0]/(1-rinv*rinv);
    for (int i=1; i<n_unknown; i++){
      lhs[i] = psi[i] + rinv*lhs[i-1];
    }
  }
  /* Old method
  int j;
  // Set bet to the first row b value (b+1/r)
  //std::complex<FTYPE> bet=b[0]; 
  std::complex<FTYPE> bet=b+rinv; 
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");

  lhs[0] = rhs[0]/bet;
  // decomposition and forward substitution 
  for (j=1;j<n_unknown-1;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
    #if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
    #endif   
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Take care of the last row
  gam[n_unknown-1] = 1.0/bet;
  bet = (b+rinv)-gam[n_unknown-1];
  lhs[n_unknown-1] - (rhs[n_unknown]-lhs[n_unknown-1])/bet;
  // backsubstitution 
  for (j=(n_unknown-2);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];  
  */
}


void tridiag_solver_complex_dirichlet_dirichlet(FTYPE boundary_condition_low,
                                                FTYPE boundary_condition_high,
                                                FTYPE b,
                                                std::complex<FTYPE> *rhs,
                                                std::complex<FTYPE> *lhs,
                                                int n_unknown)
{
  int j;
  //std::complex<FTYPE> bet=b[0]; 
  std::complex<FTYPE> bet=b;
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  gam[1] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");
  
  // Set the first element to the low boundary condition
  lhs[0] = boundary_condition_low;
  lhs[1] = (rhs[1]-lhs[0])/bet;
  // decomposition and forward substitution
  for (j=2;j<n_unknown-2;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
    #if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
    #endif 
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Take care of cell before dirichlet boundary
  gam[n_unknown-2] = 1.0/bet;
  bet = b-gam[n_unknown-2];
  lhs[n_unknown-2] = ((rhs[n_unknown-2]-boundary_condition_high)-lhs[n_unknown-3])/bet;
  // Take care of dirichlet boundary
  gam[n_unknown-1] = 1.0/bet;
  bet = b-gam[n_unknown-1];
  // Fill in the last element
  lhs[n_unknown-1] = boundary_condition_high;
  // backsubstitution 
  for (j=(n_unknown-3);j>=1;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];

/* OLD METHOD
  int j;
  //std::complex<FTYPE> bet=b[0]; 
  std::complex<FTYPE> bet=b;
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");
  
  // Set the first element to the low boundary condition
  lhs[0] = boundary_condition_low;
  
  lhs[1] = (rhs[1]-lhs[0])/bet;
  // decomposition and forward substitution 
  for (j=2;j<n_unknown-2;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
    //#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
    //#endif 
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  // Take care of last element before boundaryand dirichlet boundary
  gam[n_unknown-2] = 1.0/bet;
  bet = b-gam[n_unknown-2];
  lhs[n_unknown-2] = (rhs[n_unknown-2]-boundary_condition_high-lhs[n_unknown-2])/bet;
  // Fill in the last element
  lhs[n_unknown-1] = boundary_condition_high;
  // backsubstitution
  for (j=(n_unknown-3);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];
*/
}


// ASSUMES a and c = 1, b is a constant
void tridiag_solver_complex(FTYPE b,
                            std::complex<FTYPE> *rhs,
                            std::complex<FTYPE> *lhs,
                            int n_unknown)
{
  int j;
  //std::complex<FTYPE> bet=b[0]; 
  std::complex<FTYPE> bet=b; 
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  //if (b[0] == 0.0)
  if (b == 0.0)
    throw("Error 1 in tridiag_solver");

  lhs[0] = rhs[0]/bet;
  /* decomposition and forward substitution */
  for (j=1;j<n_unknown;j++) {
    //gam[j] = c[j-1]/bet;
    gam[j] = 1.0/bet;
    //bet=b[j]-a[j]*gam[j];
    bet=b-gam[j];
    #if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
    #endif   
    //lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
    lhs[j]=(rhs[j]-lhs[j-1])/bet;
  }
  /* backsubstitution */
  for (j=(n_unknown-2);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];  
}

void tridiag_solver_complex(std::complex<FTYPE> *a,
                            std::complex<FTYPE> *b,
                            std::complex<FTYPE> *c,
                            std::complex<FTYPE> *rhs,
                            std::complex<FTYPE> *lhs,
                            int n_unknown)
{
  int j;
  std::complex<FTYPE> bet=b[0]; 
  //std::complex<FTYPE> gam[n_unknown] = 0;
  std::complex<FTYPE> gam[n_unknown];
  gam[0] = 0;
  if (b[0] == 0.0)
    throw("Error 1 in tridiag_solver");

  lhs[0] = rhs[0]/bet;

  /* decomposition and forward substitution */
  for (j=1;j<n_unknown;j++) {
    gam[j] = c[j-1]/bet;
    bet=b[j]-a[j]*gam[j];
#if DEBUG
    if (bet==0.0) 
      throw("Error 2 in tridiag_solver");
#endif   
    lhs[j]=(rhs[j]-a[j]*lhs[j-1])/bet;
  }
  /* backsubstitution */
  for (j=(n_unknown-2);j>=0;j--) 
    lhs[j] -= gam[j+1]*lhs[j+1];  
}
