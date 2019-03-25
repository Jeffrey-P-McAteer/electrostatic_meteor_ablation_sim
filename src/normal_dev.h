
// Function to define a normal random number generator
// using the routines Ranfib (ran3_inline) and Normaldev in the 3rd edition
// of Numerical Recipes



#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC ((double) 1.0/MBIG)

inline double ran3_inline(int *idum)
{
  static int inext,inextp;
  static int ma[56];
  static int iff=0;
  int mj,mk;
  int i,ii,k;
  double x;

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  x=mj*FAC;
  return x;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

// mu is mean, sigma is standard deviation
inline double normal_dev(double mu, double sig, int *idum) { 
 
    double u,v,x,y,q;
    
    do {
      u = ran3_inline(idum);
      v = 1.7156*(ran3_inline(idum)-0.5);
      x = u - 0.449871;
      y = fabs(v) + 0.386595;
      q = x*x + y*(0.19600*y - 0.25472*x);
    } while (q > 0.27597 && (q > 0.27846 || v*v > -4.0*log(u)*u*u));
    return mu+sig*v/u;
}

