#include <cmath>
#include "eppic.h"
#include "eppic-mpi.h"

void background_dust(FArrayND &qnden, FTYPE n0)
{

  FTYPE invSigma,invSigma_x,invSigma_y;
  FTYPE gArg,gArgx,gArgy,gArgxL,gArgxR;

  switch (dust_shape) {

  case 1: // 1-D Gaussian
    invSigma = 1.0/(ny*dust_sigma);
    // FTYPE gArg;
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  gArg = (iy-ny*dust_mid)*invSigma;
	  qnden(INDICIES(ix,iy,iz)) += dust_charge*dust_den
	    *n0*exp(-0.5*gArg*gArg);
	}
      }
    }
    break;
  case 2: // 2-D symmetric Gaussian
    invSigma = 1.0/(ny*dust_sigma);
    // FTYPE gArgx,gArgy;
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  gArgx = (ix-nx*dust_mid)*invSigma;
	  gArgy = (iy-ny*dust_mid)*invSigma;
	  qnden(INDICIES(ix,iy,iz)) += dust_charge*dust_den
	    *n0*exp(-0.5*(gArgx*gArgx + gArgy*gArgy));
	}
      }
    }
    break;
  case 3: // 2 x 2-D symmetric Guassian
    invSigma = 1.0/(ny*dust_sigma);
    // FTYPE gArgxL,gArgxR;
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  gArgxL = (ix-nx*0.49)*invSigma; // TEMP HARDCODE
	  gArgxR = (ix-nx*0.51)*invSigma; // TEMP HARDCODE
	  gArgy = (iy-ny*dust_mid)*invSigma;
	  qnden(INDICIES(ix,iy,iz)) += 
	    dust_charge*dust_den*
	    n0*exp(-0.5*(gArgxL*gArgxL + gArgy*gArgy))+
	    dust_charge*dust_den*
	    n0*exp(-0.5*(gArgxR*gArgxR + gArgy*gArgy));
	}
      }
    }
    break;
  case 4: // 2-D asymmetric Gaussian across entire box
    invSigma_x = 1.0/(nx*nsubdomains*dust_sigma_x);
    invSigma_y = 1.0/(ny*dust_sigma_y);
    // FTYPE gArgx,gArgy;
    int xshift = nx*subdomain.id_number;
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  gArgx = (ix+xshift-nx*nsubdomains*dust_mid)*invSigma_x;
	  gArgy = (iy-ny*dust_mid)*invSigma_y;
	  qnden(INDICIES(ix,iy,iz)) += dust_charge*dust_den
	    *n0*exp(-0.5*(gArgx*gArgx + gArgy*gArgy));
	}
      }
    }
    break;

  } // END case

}
