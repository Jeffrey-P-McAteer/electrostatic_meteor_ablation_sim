
#include "eppic-velsample.h"
#include "eppic-math.h"

PTYPE sampleVel(velsample &veldev)
{
  // the -1 in the next line is not exactly right
  FTYPE sample = ran3(&veldev.ran3Seed)*(veldev.inv_cdf.size()-1);
  int inv_idx = static_cast<int>(sample);
  return veldev.dv*(veldev.inv_cdf(inv_idx) + 
		    (veldev.inv_cdf(inv_idx+1)-
		     veldev.inv_cdf(inv_idx))*(sample-inv_idx)
		    );

}
