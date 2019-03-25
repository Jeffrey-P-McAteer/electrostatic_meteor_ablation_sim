


#ifdef EPPIC_FFTW_USE_D_PREFIX
#include "drfftw.h"

#else

#ifdef use_fftw3
#include "fftw3.h"
#define rfftwnd_plan fftw_plan
#define rfftwnd_create_plan rfftwnd_create_plan_1
#else
#include "rfftw.h"
#endif

#endif
