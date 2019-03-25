//----------------------------------------------------

// A collection of random number generators inside of C++ structures.  
// These can be included as inlined functions.

// Meers Oppenheim June 2009
// 

//----------------------------------------------------

#ifndef RANDOM_H
#define RANDOM_H

struct RanFast {
    // Based on Ranq1 in NR.
    unsigned long long int v;
    RanFast(unsigned long long int j) : v(4101842887655102017LL) {
	v ^= j;
	v = int64();
    }
    inline unsigned long long int int64() {
	v ^= v >> 21; 
	v ^= v << 35; 
	v ^= v >> 4;
	return v * 2685821657736338717LL;
    }
    inline double dbl() { 
	return 5.42101086242752217E-20 * int64(); 
    }
    inline unsigned int int32() { 
	return (unsigned int)int64(); 
    }
};

struct GaussDev : RanFast {
    // Based on the NR 3.0 Normaldev routine
    GaussDev(unsigned long long int i)
	: RanFast(i){}
    inline double dev() {
	double u,v,x,y,q;
	do {
	    u = dbl();
	    v = 1.7156*(dbl()-0.5);
	    x = u - 0.449871;
	    y = fabs(v) + 0.386595;
	    q = x*x + y*(0.19600*y-0.25472*x);
	} while (q > 0.27597
		 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
	return v/u;
    }
};

#endif
