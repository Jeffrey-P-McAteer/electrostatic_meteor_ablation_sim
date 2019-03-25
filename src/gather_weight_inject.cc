/* Multiple gather-scale routines below for 1-D, 2-D 4-D and 6-D */

#include "eppic.h"
#include<cmath>   
//***************************************************************************

/* Calculate the 1-D distribution function on a grid by gathering particle 
   characteristics and scaling the solution by scale.  
   xmax and xmin is the range which fits into the nmesh points.
   (points outside nmesh are discarded) except:
   iwrap specifies the integer values at which the mesh wraps one grid cell
     from the max of the mesh to the bottom

   The incoming array is not zeroed!

*/

const int n_dist1=1;

void gather_weight_inject(ArrayNd<OTYPE, n_dist1> &den, particle_dist &pic, 
		   PTYPEAVec &weight, FTYPE *xmin, FTYPE *xmax, 
			  int *nmesh, int *iwrap, FTYPE scale, FTYPE nxmin)
{
  
  // Local variables 
  FTYPE w;
  const int np=pic.np;
  int ix[n_dist1][2];
  FTYPE wx[n_dist1][2], xscale[n_dist1];
  //Define a pointer to access the position arrays
  PTYPEAVec *data[NDIM]={INDICIES(&(pic.x), &(pic.y), &(pic.z))};

  for (int id = 0; id<n_dist1; ++id) {
    FTYPE range=(xmax[id]-xmin[id]);
    if (range >0) xscale[id]=nmesh[id]/range;
    else xscale[id]=0;
  }

  /* 
     nxmin is typically -1 for inject @ LHS, but if averaging, should be
     -.5 for navg = 2, etc, which is accounted for in xscale.
  */
  nxmin*=xscale[0];

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {
    // Define the four nearest grid points: 
    for (int id = 0; id<n_dist1; ++id) {
      FTYPE xx=((*(data[id]))(i)-xmin[id]) * xscale[id];

      if (xx >= 0) ix[id][0] = (int) xx;
      else ix[id][0] = (int) xx - 1;
      ix[id][1] = ix[id][0] + 1;

      // Wrap periodic dimensions:  
      // This routine only wraps the top element to the 0 element:
      // (If one nmesh[dim]=1, this can also be used to obtain all info
      // for that point)
      if (ix[id][1] == iwrap[id]) ix[id][1] = 0;

      // and the corresponding linear weighting factors: 
      wx[id][1] = xx - ix[id][0];
      wx[id][0] = 1. - wx[id][1];
    }

    // Loop through the mesh assigning particle fractions.
    // Skip out-of-bounds particles
    for (int i0= 0; i0<2; i0++)
      if (ix[0][i0] >= 0 && ix[0][i0] < nmesh[0])
	{
	  w=wx[0][i0];
	  den(ix[0][i0]) += w*weight[i]*scale;
	}
      else if (i0 == 0 && ix[0][i0] >= nxmin &&  ix[0][i0] < nmesh[0]) 
	{
	  w=fabs(wx[0][1]);
	  den(ix[0][1]) += w*weight[i]*scale;
	}
  } // end for (i = 0; i < np; ++i) 

  
} // End 1-D gather 
//***************************************************************************

/* Calculate the 2-D distribution function on a grid by gathering particle 
   characteristics and scaling the solution by scale.  
   xmax and xmin is the range which fits into the nmesh points.
   (points outside nmesh are discarded) except:
   iwrap specifies the integer values at which the mesh wraps one grid cell
     from the max of the mesh to the bottom

   The incoming array is not zeroed!

*/

const int n_dist2=2;

void gather_weight_inject(ArrayNd<OTYPE, n_dist2> &den, particle_dist &pic, 
		   PTYPEAVec &weight,  FTYPE *xmin, FTYPE *xmax, 
			  int *nmesh, int *iwrap, FTYPE scale, FTYPE nxmin)
{
  
  // Local variables 
  FTYPE w;
  const int np=pic.np;
  int ix[n_dist2][2];
  FTYPE wx[n_dist2][2], xscale[n_dist2];
  //Define a pointer to access the position arrays
  PTYPEAVec *data[NDIM]={INDICIES(&pic.x, &pic.y, &pic.z)};

  for (int id = 0; id<n_dist2; ++id) {
    FTYPE range=(xmax[id]-xmin[id]);
    if (range >0) xscale[id]=nmesh[id]/range;
    else xscale[id]=0;
  }

  /* 
     nxmin is typically -1 for inject @ LHS, but if averaging, should be
     -.5 for navg = 2, etc, which is accounted for in xscale.
  */
  nxmin*=xscale[0];

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {
    // Define the four nearest grid points: 
    for (int id = 0; id<n_dist2; ++id) {
      FTYPE xx=((*data[id])(i)-xmin[id]) * xscale[id];

      if (xx >= 0) ix[id][0] = (int) xx;
      else ix[id][0] = (int) xx - 1;
      ix[id][1] = ix[id][0] + 1;

      // Wrap periodic dimensions:  
      // This routine only wraps the top element to the 0 element:
      // (If one nmesh[dim]=1, this can also be used to obtain all info
      // for that point)
      if (ix[id][1] == iwrap[id]) ix[id][1] = 0;

      // and the corresponding linear weighting factors: 
      wx[id][1] = xx - ix[id][0];
      wx[id][0] = 1. - wx[id][1];
    }

    // Loop through the mesh assigning particle fractions.
    // Skip out-of-bounds particles
    for (int i0= 0; i0<2; i0++)
      if (ix[0][i0] >= 0 && ix[0][i0] < nmesh[0]) {
	
	for (int i1= 0; i1<2; i1++)
	  if (ix[1][i1] >= 0 && ix[1][i1] < nmesh[1])
	    {
	      w=wx[0][i0]*wx[1][i1];
	      den(ix[0][i0],ix[1][i1]) += w*weight[i]*scale;
	    }
      }
      else if (i0 == 0 && ix[0][i0] >= nxmin &&  ix[0][i0] < nmesh[0]) 
	{
	for (int i1= 0; i1<2; i1++)
	  if (ix[1][i1] >= 0 && ix[1][i1] < nmesh[1])
	    {
	      w=fabs(wx[0][i0]*wx[1][i1]);
	      den(ix[0][i0],ix[1][i1]) += w*weight[i]*scale;
	    }
	}
  } // end for (i = 0; i < np; ++i) 

  
} // End 2-D gather 
   
//***************************************************************************


/* Calculate the 3-D distribution function on a grid by gathering particle 
   characteristics and scaling the solution by scale.  
   xmax and xmin is the range which fits into the nmesh points.
   (points outside nmesh are discarded) except:
   iwrap specifies the integer values at which the mesh wraps one grid cell
     from the max of the mesh to the bottom

   The incoming array is not zeroed!

*/

const int n_dist3=3;

void gather_weight_inject(ArrayNd<OTYPE, n_dist3> &den, particle_dist &pic, 
		   PTYPEAVec &weight,  FTYPE *xmin, FTYPE *xmax, 
			  int *nmesh, int *iwrap, FTYPE scale, FTYPE nxmin)
{
  
  // Local variables 
  FTYPE w;
  const int np=pic.np;
  int ix[n_dist3][2];
  FTYPE wx[n_dist3][2], xscale[n_dist3];
  //Define a pointer to access the position arrays
  PTYPEAVec *data[NDIM]={INDICIES(&pic.x, &pic.y, &pic.z)};

  for (int id = 0; id<n_dist3; ++id) {
    FTYPE range=(xmax[id]-xmin[id]);
    if (range >0) xscale[id]=nmesh[id]/range;
    else xscale[id]=0;
  }

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {
    // Define the four nearest grid points: 
    for (int id = 0; id<n_dist3; ++id) {
      FTYPE xx=((*data[id])(i)-xmin[id]) * xscale[id];

      if (xx >= 0) ix[id][0] = (int) xx;
      else ix[id][0] = (int) xx - 1;
      ix[id][1] = ix[id][0] + 1;

      // Wrap periodic dimensions:  
      // This routine only wraps the top element to the 0 element:
      // (If one nmesh[dim]=1, this can also be used to obtain all info
      // for that point)
      if (ix[id][1] == iwrap[id]) ix[id][1] = 0;

      // and the corresponding linear weighting factors: 
      wx[id][1] = xx - ix[id][0];
      wx[id][0] = 1. - wx[id][1];
    }

    // Loop through the mesh assigning particle fractions.
    // Skip out-of-bounds particles
    for (int i0= 0; i0<2; i0++)
      if (ix[0][i0] >= 0 && ix[0][i0] < nmesh[0])
	
	for (int i1= 0; i1<2; i1++)
	  if (ix[1][i1] >= 0 && ix[1][i1] < nmesh[1])
	    
	    for (int i2= 0; i2<2; i2++)
	      if (ix[2][i2] >= 0 && ix[2][i2] < nmesh[2])
		{
		  w=wx[0][i0]*wx[1][i1]*wx[2][i2];
		  den(ix[0][i0],ix[1][i1],ix[2][i2]) += w*weight[i]*scale;
		}

  } // end for (i = 0; i < np; ++i) 

  
} // End 3-D gather 

/* Calculate the 4-D distribution function on a grid by gathering particle 
   characteristics and scaling the solution by scale.  
   xmax and xmin is the range which fits into the nmesh points.
   (points outside nmesh are discarded) except:
   iwrap specifies the integer values at which the mesh wraps one grid cell
     from the max of the mesh to the bottom

   The incoming array is not zeroed!

*/

const int n_dist4=4;

void gather_weight_inject(ArrayNd<OTYPE, n_dist4> &den, particle_dist &pic, 
		   PTYPEAVec &weight,  FTYPE *xmin, FTYPE *xmax, 
		   int *nmesh, int *iwrap, FTYPE scale, FTYPE nxmin)
{
  
  // Local variables 
  FTYPE w;
  const int np=pic.np;
  int ix[n_dist4][2];
  FTYPE wx[n_dist4][2], xscale[n_dist4];
  //Define a pointer to access the position arrays
  PTYPEAVec *data[NDIM]={INDICIES(&pic.x, &pic.y, &pic.z)};

  for (int id = 0; id<n_dist4; ++id) {
    FTYPE range=(xmax[id]-xmin[id]);
    if (range >0) xscale[id]=nmesh[id]/range;
    else xscale[id]=0;
  }

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {
    // Define the four nearest grid points: 
    for (int id = 0; id<n_dist4; ++id) {
      FTYPE xx=((*data[id])(i)-xmin[id]) * xscale[id];

      if (xx >= 0) ix[id][0] = (int) xx;
      else ix[id][0] = (int) xx - 1;
      ix[id][1] = ix[id][0] + 1;

      // Wrap periodic dimensions:  
      // This routine only wraps the top element to the 0 element:
      // (If one nmesh[dim]=1, this can also be used to obtain all info
      // for that point)
      if (ix[id][1] == iwrap[id]) ix[id][1] = 0;

      // and the corresponding linear weighting factors: 
      wx[id][1] = xx - ix[id][0];
      wx[id][0] = 1. - wx[id][1];
    }

    // Loop through the mesh assigning particle fractions.
    // Skip out-of-bounds particles
    for (int i0= 0; i0<2; i0++)
      if (ix[0][i0] >= 0 && ix[0][i0] < nmesh[0])
	
	for (int i1= 0; i1<2; i1++)
	  if (ix[1][i1] >= 0 && ix[1][i1] < nmesh[1])
	    
	    for (int i2= 0; i2<2; i2++)
	      if (ix[2][i2] >= 0 && ix[2][i2] < nmesh[2])

		for (int i3= 0; i3<2; i3++)
		  if (ix[3][i3] >= 0 && ix[3][i3] < nmesh[3])
		{
		  w=wx[0][i0]*wx[1][i1]*wx[2][i2]*wx[3][i3];
		  den(ix[0][i0],ix[1][i1],ix[2][i2],ix[3][i3]) += 
		    w*weight[i]*scale;
		}

  } // end for (i = 0; i < np; ++i) 

  
} // End 4-D gather 


/* Calculate the 6-D distribution function on a grid by gathering particle 
   characteristics and scaling the solution by scale.  
   xmax and xmin is the range which fits into the nmesh points.
   (points outside nmesh are discarded) except:
   iwrap specifies the integer values at which the mesh wraps one grid cell
     from the max of the mesh to the bottom

   The incoming array is not zeroed!

*/

const int n_dist6=6;

void gather_weight_inject(ArrayNd<OTYPE, n_dist6> &den, particle_dist &pic, 
		   PTYPEAVec &weight,  FTYPE *xmin, FTYPE *xmax, 
		   int *nmesh, int *iwrap, FTYPE scale, FTYPE nxmin)
{
  
  // Local variables 
  FTYPE w;
  const int np=pic.np;
  int ix[n_dist6][2];
  FTYPE wx[n_dist6][2], xscale[n_dist6];
  //Define a pointer to access the position arrays
  PTYPEAVec *data[NDIM]={INDICIES(&pic.x, &pic.y, &pic.z)};

  for (int id = 0; id<n_dist6; ++id) {
    FTYPE range=(xmax[id]-xmin[id]);
    if (range >0) xscale[id]=nmesh[id]/range;
    else xscale[id]=0;
  }

  //  For each particle ... 
  for (int i = 0; i < np; ++i) {
    // Define the four nearest grid points: 
    for (int id = 0; id<n_dist6; ++id) {
      FTYPE xx=((*data[id])(i)-xmin[id]) * xscale[id];

      if (xx >= 0) ix[id][0] = (int) xx;
      else ix[id][0] = (int) xx - 1;
      ix[id][1] = ix[id][0] + 1;

      // Wrap periodic dimensions:  
      // This routine only wraps the top element to the 0 element:
      // (If one nmesh[dim]=1, this can also be used to obtain all info
      // for that point)
      if (ix[id][1] == iwrap[id]) ix[id][1] = 0;

      // and the corresponding linear weighting factors: 
      wx[id][1] = xx - ix[id][0];
      wx[id][0] = 1. - wx[id][1];
    }

    // Loop through the mesh assigning particle fractions.
    // Skip out-of-bounds particles
    for (int i0= 0; i0<2; i0++)
    if (ix[0][i0] >= 0 && ix[0][i0] < nmesh[0])
	
      for (int i1= 0; i1<2; i1++)
      if (ix[1][i1] >= 0 && ix[1][i1] < nmesh[1])
	    
	for (int i2= 0; i2<2; i2++)
	if (ix[2][i2] >= 0 && ix[2][i2] < nmesh[2])

	  for (int i3= 0; i3<2; i3++)
	  if (ix[3][i3] >= 0 && ix[3][i3] < nmesh[3])
		    
	    for (int i4= 0; i4<2; i4++)
	    if (ix[4][i4] >= 0 && ix[4][i4] < nmesh[4])

	      for (int i5= 0; i5<2; i5++)
	      if (ix[5][i5] >= 0 && ix[5][i5] < nmesh[5])
		{
		  w=wx[0][i0]*wx[1][i1]*wx[2][i2]*wx[3][i3]*
		    wx[4][i4]*wx[5][i5];
		  den(ix[0][i0], ix[1][i1], ix[2][i2], ix[3][i3],
		      ix[4][i4], ix[5][i5]) += w*weight[i]*scale;
		}

  } // end for (i = 0; i < np; ++i) 

  
} // End 6-D gather 
