

#include "eppic-velsample.h"
#include "eppic-math.h"
#include "eppic-misc.h"


/* Indefinite integral of v*f(v)dv
   f(v) = 1/sqrt(2*PI*vth**2)*exp(-.5*[(v-vnot)/vth]**2)
*/
inline FTYPE cdffunc(FTYPE v,FTYPE vth,FTYPE v0) {
  return (-vth/sqrt(2*PI)*exp(-Sqr((v-v0)/vth)/2)
	  +v0/2*erf(1/sqrt(2)*(v-v0)/vth));
    }

inline FTYPE cdffunc_vpos(FTYPE v,FTYPE vth,FTYPE v0) {
  /*
  FTYPE erf_nov = erf(v0/(sqrt(2)*vth));
  FTYPE normalization = (v0+exp(-v0*v0/(2*vth*vth))*sqrt(2/M_PI)*vth+
                         v0*erf_nov)/2;
  FTYPE vmv0 = v-v0;
  return ((erf_nov+erf(vmv0/(sqrt(2)*vth)))*v0/2+
          vth*(exp(-v0*v0/(2*vth*vth))-exp(-vmv0*vmv0/(2*vth*vth)))/(sqrt(2)*M_PI))/
    normalization;
  */
  FTYPE erf_nov = erf(v0/(sqrt(2.0)*vth));
  FTYPE normalization = vth*(v0*sqrt(M_PI/2.0)*(1+erf_nov) + vth*exp(-v0*v0/(2.0*vth*vth)));
  FTYPE vmv0 = v-v0;
  return vth*(vth*(exp(-v0*v0/(2.0*vth*vth))-exp(-vmv0*vmv0/(2.0*vth*vth))) + v0*sqrt(M_PI/2.0)*(erf_nov-erf(-vmv0/(sqrt(2.0)*vth))))/normalization;
}

inline FTYPE cdffunc_vneg(FTYPE v,FTYPE vth,FTYPE v0) {
  /*
  FTYPE erf_nov = erf(v0/(sqrt(2)*vth));
  FTYPE normalization = v0*(1-erf_nov)/2-exp(-v0*v0/(2*vth*vth))*vth/(sqrt(2*M_PI));
  FTYPE vmv0 = v-v0;
  //return (v0*(1+erf(vmv0/(sqrt(2)*vth))) - exp(-vmv0*vmv0/(2*vth*vth))*sqrt(2/M_PI))/normalization;
  return (v0*(1+erf(vmv0/(sqrt(2)*vth))) - exp(-vmv0*vmv0/(2*vth*vth))*sqrt(2/M_PI))/
    (normalization*2);
  */
  /*
  FTYPE erf_nov = erf(-v0/(sqrt(2)*vth));
  //FTYPE normalization = v0*sqrt(M_PI/2)*(1-erf_nov)-exp(-v0*v0/(2*vth*vth));
  FTYPE normalization = vth*(-v0*sqrt(M_PI/2)*(1+erf_nov) + vth*exp(-v0*v0/(2*vth*vth)));
  FTYPE vmv0 = v-v0;
  return ((v0*sqrt(M_PI/2)*(1+erf(vmv0/(sqrt(2)*vth))) - exp(-vmv0*vmv0/(2*vth*vth))))/normalization;
  */

  FTYPE erf_nov = erf(v0/(sqrt(2.0)*vth));
  FTYPE vmv0 = v-v0;
  FTYPE normalization = sqrt(M_PI/2.0)*vth*v0*(1-erf_nov) - exp(-v0*v0/(2.0*vth*vth))*vth*vth;
  return (vth*vth*(exp(-vmv0*vmv0/(2.0*vth*vth))-exp(-v0*v0/(2.0*vth*vth))) + 
          sqrt(M_PI/2.0)*v0*vth*(erf(-vmv0/(sqrt(2.0)*vth))-erf_nov))/normalization;
}

void init_velsample(velsample &veldev, PTYPE vth, PTYPE v0, 
		    int precision, int direction, int seed) 
{

  veldev.ran3Seed = seed;
  PTYPE &dv = veldev.dv;
  ArrayNd<FTYPE,1> &inv_cdf = veldev.inv_cdf;
  FTYPE maxSpeed = 10*vth;
  if (direction>0) {
    direction = 1;
    maxSpeed+=v0;
    if (maxSpeed<0) maxSpeed = vth;
  } else {
    direction = -1;
    maxSpeed-=v0;
    if (direction*maxSpeed>0) maxSpeed = vth;
  }

  dv=pow(10,log10(maxSpeed)-precision)*direction;
  int ncdf = static_cast<int>(pow(10,precision));
  ArrayNd<FTYPE,1> cdf= ArrayNd<FTYPE,1>(ncdf);
  if (dv>0) // 0 domain, inject to the left
    {
      // CHANGED BY GLENN TO TEST A CDF OBTAINED VIA MATHEMATICA
      // Remove the FvZero and FvNorm terms because they are confusing
      // Constain entire cdf inthe cdf function
      /*
      // build dist specfic pos and neg cdf from analytic form 
      FTYPE FvZero = -vth/sqrt(2*PI)*exp(-Sqr(v0/vth)/2)-
	v0/2*erf(1/sqrt(2)*v0/vth);
      FTYPE FvNorm = cdffunc(ncdf*dv,vth,v0)-FvZero;
      for (int i=0;i<ncdf;i++) {
	cdf(i)=(cdffunc(i*dv,vth,v0)-FvZero)/FvNorm;
      }
      */
      cdf(0) = 0.0;
      for (int i=1;i<ncdf;i++) {
        cdf(i) = cdffunc_vpos(i*dv,vth,v0);
      }

      // invert the cdfs for the dist specfic number of injected particles 
      invert_function(cdf,inv_cdf,ncdf);
    }
  else // ndomain - 1 domain, inject to the right 
    {
      // CHANGED BY GLENN TO TEST A CDF OBTAINED VIA MATHEMATICA
      /*
      // build dist specfic pos and neg cdf from analytic form 
     FTYPE FvZero = -vth/sqrt(2*PI)*exp(-Sqr(v0/vth)/2)-
	v0/2*erf(1/sqrt(2)*v0/vth);
     FTYPE FvNorm = cdffunc(ncdf*dv,vth,v0)-FvZero; // note the -vx0d, diff from above
      for (int i=0;i<ncdf;i++) {
	cdf(i)=(cdffunc(i*dv,vth,v0)-FvZero)/FvNorm;
      }
      */
      cdf(0) = 0.0;
      for (int i=1;i<ncdf;i++) {
        //cdf(i) = cdffunc_vneg((ncdf-i)*dv,vth,v0);
        //cdf(i) = 1-cdffunc_vneg(i*dv,vth,v0);
        cdf(i) = cdffunc_vneg(i*dv,vth,v0);
      }      
      // invert the cdfs for the dist specfic number of injected particles 
      invert_function(cdf,inv_cdf,ncdf);
    } // 0 or nsubdomains -1 if statement
}
