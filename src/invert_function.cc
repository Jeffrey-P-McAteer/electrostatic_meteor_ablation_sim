/** This takes a vector F(x), calculated for a fixed step size in x and 
    inverts the function along a fixed step size in F, producing 
    Inv_F(F).
**/


#include "eppic-types.h"
#include "eppic-mpi.h"


void invert_function(FTYPEAVec &F, FTYPEAVec &Inv_F, int nsteps) {
  
  if (Inv_F.size()!=nsteps)
    Inv_F = FTYPEAVec(nsteps);
  int nx = F.size();
  FTYPE delta_inv_F=F(nx-1)/FTYPE(nsteps-1);
  FTYPE curr_delta = 0;
  // assume dx=1
  // loop over ix until F_inv(if)+ix(+) > 
  int ix=1;
  Inv_F(0)=0;
  int iinv=1;
  while((iinv<nsteps)&&(ix<nx)) {
    curr_delta = F(ix)-(iinv-1)*delta_inv_F;
    if (curr_delta>delta_inv_F) {
      FTYPE slope=(F(ix)-F(ix-1));
      FTYPE ix_fraction = (delta_inv_F*iinv-F(ix-1))/slope;
      Inv_F(iinv) = ix-1+ix_fraction;
      iinv++;
    } else {
      ix++;
    }
  }
  Inv_F(nsteps-1)=nx-1;
}
