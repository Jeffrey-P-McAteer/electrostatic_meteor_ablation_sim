

#include "efield_inject_tridiag.h" /* to know about field, types, fftw stuff */
#include "eppic-mpi.h"

void efield_inject_tridiag_init(FArrayND &rho,
				field &Efield,
				rfftw_many &working_array)
{

  //    rho dc calculation:
  //    This is a tricky part of the calculation. The k=0 component of rho,
  //    the charge density must be zero for specific cases. It can be set to 
  //    zero for (I think) all cases, but this would change the meaning of the 
  //    solution. The foundation for deciding when to set the k=0 values (ie
  //    the average of rho) to zero came from reading Birdsall and Langdon 
  //    (First Edition), section 4-11: Boundary conditions for charge, current,
  //    field and potential.

  //    Here they start by  describing the periodic case, where E(x=0)=E(x=Nx)
  //    as is also the case for the potential and the charge density. Equation
  //    (1) shows that
  //    $\int_{x}^{x+L}\frac{\partial E}{\partial x} dx = E(x+L)-E(x)=
  //    L/\epsilon <\rho>$ and since E is periodic, then the average of 
  //    the charge density ($<\rho>$) must also be zero. 

  //    This analysis may also be extended to the other possible boundary 
  //    conditions. For example, Neumann boundary conditions set E=constant 
  //    on the boundaries, and currently we only allow E=0 on the boundaries. 
  //    With this condition, the average charge density should also be zero. 

  //    For Dirichlet and Open (aka Natural) boundary conditions, there is no
  //    constraint on the electric field on the boundaries, only the potential.
  //    There for the average should be left untouched, and it's value will
  //    play a part in a meaningful solution. 
     
  //    For mixed boundary conditions (in x), the problem becomes a bit more
  //    complicated. Obviously you can't have a mixed boundary condition for 
  //    periodic cases, but you could have a Neumann boundary condition on 
  //    one boundary with either an Open or Dirichlet boundary on the other. 
  //    While the Neumann boundary on one end fixes E, the other end should be
  //    determined only by the solution to Poisson's equation. Take for example
  //    E(x) = 0, for a Neumann boundary condition on the LHS, then 
  //    E(x+L)-E(x) = E(x+L) = L/eps <rho>. 

  //    NOTE: This mixed boundary case has not been thoroughly tested, and as a
  //    precaution <rho> is set to zero when ever Neumann is chosen for either 
  //    boundary condition. 


  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

  FTYPE rho_dc=0;
  if (Efield.zero_dc) {
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<ny;iy++) 
	for(int iz=0;iz<nz;iz++) 
	  rho_dc+=rho(INDICIES(ix,iy,iz));
    rho_dc/=(nx*ny*nz);

    if (nsubdomains>1)  {
      FTYPE rho_dc_tmp=rho_dc;
      MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
		    subdomain.neighbor_comm);
      rho_dc/=nsubdomains;
    }
  }


  int iy_start=0;
  int iy_end=ny;
#if NDIM==3
  iy_start=working_array.x_start;
  iy_end=iy_start+working_array.local_nx;
  if (! working_array.allocatable()) { 
    return;
  }
#endif
  // see headers for INDICIES and RFFTW_INDICIES defs
  for (int ix=0;ix<nx;ix++)
    for(int iy=iy_start;iy<iy_end;iy++) 
      for(int iz=0;iz<nz;iz++) 
	working_array.data(RFFTW_INDICIES(ix,iy,iz)) = 
	  rho(INDICIES(ix,iy,iz))-rho_dc;
}
