#include "eppic-velsample.h"
#include "eppic-math.h"
#include "eppic-misc.h"


/* Indefinite integral of v*f(v)dv
   f(v) = 1/(2*PI)^1.5*vth^3)*exp(-.5*[v/vth]^2)
*/
inline FTYPE cdffunc_ablating(FTYPE v,FTYPE vth,FTYPE v0) {
  //return 1/(4*PI*pow(vth,2))*(erf(v/(sqrt(2)*vth))+1);
  //  return 1/(4*PI*pow(vth*sqrt(2),2))*(erf(v/(2*vth))+1);
  //  return 1/(sqrt(2*PI)*vth)*(erf(v/(sqrt(2)*vth))+1);
  //  return 1/(sqrt(2*PI)*vth)*(erf(v/(sqrt(2)*vth))+1);
  //return (1/(4*PI*vth*vth))*(1+erf(v/(sqrt(2)*vth)));
  //return erf(v/(sqrt(2.0)*vth)) - (sqrt(2.0/PI)*v*exp(-v*v/(2*vth*vth)))/vth;
  return 1 - (v*v/(2.0*vth*vth) + 1)*exp(-v*v/(2*vth*vth)); //v^3/(2*vth^4)*exp(-v^2/2*vth^2) pdf
  //return erf(v/(sqrt(2.0)*vth)) - sqrt(2.0/PI)*v*exp(-v*v/(2*vth*vth))/vth; //v^2 exp(-v^2/2*vth^2)
}

void init_velsample_ablating(velsample &veldev, PTYPE vth, PTYPE v0, 
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
      // build dist specfic pos and neg cdf from analytic form 
      for (int i=0;i<ncdf;i++) {
	cdf(i)=cdffunc_ablating(i*dv,vth,v0);
      }
      // invert the cdfs for the dist specfic number of injected particles 
      invert_function(cdf,inv_cdf,ncdf);
    }
  else // ndomain - 1 domain, inject to the right 
    {
      // build dist specfic pos and neg cdf from analytic form 
      for (int i=0;i<ncdf;i++) {
	cdf(i)=cdffunc_ablating(i*dv,vth,v0);
      }
      // invert the cdfs for the dist specfic number of injected particles 
      invert_function(cdf,inv_cdf,ncdf);
    } // 0 or nsubdomains -1 if statement
}
