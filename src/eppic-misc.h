


#ifndef EPPIC_MISC_H
#define EPPIC_MISC_H

#include "eppic-types.h"

/** This inverts a function F(x) at nsteps discrete points to produce Inv_F */
void invert_function(FTYPEAVec &F, FTYPEAVec &Inv_F, int nsteps);

/** Calculates the Cumulative Distribution Function for a particular 
    probability distribution function.
**/

void calc_cdf(FTYPEAVec &dist, FTYPEAVec &cdf);



/** Calculated the max value and max local average of a 1D function
    across domains. 
**/

void prob_dist_1D_ddecomposed(FTYPEAVec &local_dist,
			      FTYPE &global_max,
			      FTYPE &global_avg_max, 
			      FTYPE &right_total, 
			      FTYPE &left_total);

#endif
