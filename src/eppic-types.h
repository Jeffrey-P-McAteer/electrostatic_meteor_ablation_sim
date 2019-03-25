


#ifndef SRCCONFIG_H
#define SRCCONFIG_H
#include "srcconfig.h"
#endif

#ifndef EPPIC_TYPES_H
#define EPPIC_TYPES_H


// This file is used to define basic types throughout the code:

// Define the fundamental field type: (default for most calculations)
// One defines the appropriate flag,  F_PRECISION = 1,2 and P_PRECISION = 1,2
// the type is then defined (1 is for single precision, 2 for double)

#include "srcconfig.h"
#include <string>
#include <vector>

#define stringVec std::vector<std::string>

#if F_PRECISION == 2
#define FTYPE double
#else

#if F_PRECISION == 1
#define FTYPE float
#else
// Default - Make sure you switch both terms
#define F_PRECISION 2
#define FTYPE double
#endif

#endif

// Define the particle type: the huge particle arrays 
//                           do not need high presision

#if P_PRECISION == 2
#define PTYPE double
#else

#if P_PRECISION == 1
#define PTYPE float
#else
// Default - Make sure you switch both terms
#define P_PRECISION 1
#define PTYPE float
#endif
#endif

// Define an output floating type - usually the lowest precision available
#define OTYPE float
#include "ArrayNd.h"

#define intAVec ArrayNd<int,1>

#if P_PRECISION == 1
#define PTYPEAVec ArrayNd<float,1>
#else
#define PTYPEAVec ArrayNd<double,1>
#endif

#if F_PRECISION == 1
#define FTYPEAVec ArrayNd<float,1>
#else
#define FTYPEAVec ArrayNd<double,1>
#endif

typedef ArrayNd<FTYPE,NDIM> FArrayND;
#ifdef USE_DOMAINS
#include "ArrayNd_ranged.h"
typedef ArrayNd_ranged<FTYPE,NDIM> FArrayND_ranged;
#else
typedef ArrayNd<FTYPE,NDIM> FArrayND_ranged;
#endif


#define MAX_NDIM 3

#ifndef NDIM
#define NDIM 2
#endif

#if NDIM == 1
#define INDICIES(ix,iy,iz) ix
#elif NDIM == 2
#define INDICIES(ix,iy,iz) ix, iy
#else
#define INDICIES(ix,iy,iz) ix, iy, iz
#endif

#if NDIM == 1
#define DONDIM(ix,iy,iz) ix;
#elif NDIM == 2
#define DONDIM(ix,iy,iz) ix; iy;
#else
#define DONDIM(ix,iy,iz) ix; iy; iz;
#endif

#define DO_VELDIM(veldim,ix,iy,iz) if(veldim==1){ix;}	\
  elseif(veldim==2){iy;}				\
  else{iz;}


#endif
