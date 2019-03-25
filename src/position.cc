#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

void position(particle_dist &dist, int id, int &iran)
{
  
  // This routine is for customized particle distributions.



  int i;
  int prime(int);
  FTYPE linear(FTYPE x);
  void gasuniform(PTYPEAVec &, PTYPEAVec &, int, int, int);
  void shuffle (PTYPEAVec &x, int &seed);
  void pwluniform(PTYPEAVec &vec, FTYPE x0, FTYPE x1, int &,
		  FTYPE (den_func)(FTYPE x));
  FTYPE gasdev(int &);
  int x_start=nx*subdomain.id_number;

  // A distribution with a gassian bump in the middle:
  // param3 is magnitude of the bump above dist.n0avg
  // param4 is the width of the bump 
  /*  
  //use to be in x direction,here's the code:
  FTYPE fraction_in_bump=param3[id]*dist.n0avg*param4[id]*sqrt(3.14159)/
    (param3[id]*dist.n0avg*param4[id]*sqrt(3.14159)+dist.n0avg*nx*dx);
  for (i = 0; i < dist.np*fraction_in_bump; ++i) 
    dist.x[i] = fmod(gasdev(iran) * param4[id] / dx + nx/2., FTYPE(nx));
  for (i = i; i < dist.np; ++i) 
    dist.x[i] = ran3(&iran) * nx;
  if (ndim >= 2)
    for (i = 0; i < dist.np ; ++i)  dist.y[i] = ran3(&iran) * ny;
  */
  
  // make gaussian in y direction

  // make an error if param3[id] < 1
  //  FTYPE fraction_in_bump=param3[id]/(param3[id]+1.0);
   // old
  // n_gauss: number of ions under guassian
  //          ie area of gaussian*n0avg*param3[id]
  //          area of gaussian = param4[id]*dy*SQRT2PI
  //          because variance = SQR(param4[id]*dy)
  FTYPE n_gauss = param3[id]*dist.n0avg*param4[id]*SQRT2PI;
  FTYPE fraction_in_bump=n_gauss/(n_gauss+dist.n0avg*FTYPE(ny)*dy);
  dist.n0avg = n_gauss/(dy*FTYPE(ny))+dist.n0avg;
  if (mpi_rank == 0) {
    printf("\tUPDATE: Using Guassian initial conditions for PIC\n\tParticles represented by single PIC particle: %g\n",
	   (dist.n0avg/(FTYPE(dist.np*mpi_np))*nx*nsubdomains*dx*ny*dy*nz*dz));
    printf("\tNew N0: %g\n\n\n",
	   dist.n0avg);
  }
  
  // for (i = 0; i < dist.np ; ++i)  dist.x[i] = gasdev(iran) * param4[id]/dx + nx/2.;
  // if (ndim >= 2) {
  //   for (i = 0; i < fraction_in_bump*dist.np; ++i) {
  //     dist.y[i] = gasdev(iran) * param4[id]/dy + ny/2.;
  //     if (dist.y[i] > FTYPE(ny) || dist.y[i] < 0.0) 
  // 	dist.y[i] = ran3(&iran) * ny;
  // 	//      dist.y[i] = -1.0;
  //   }
  //       for (i = i; i < dist.np; ++i) 
  //     dist.y[i] = ran3(&iran) * ny;
  // }
  /*for (i=0; i<dist.np; i++) dist.x[i] = ran3(&iran)*nx;
  if (ndim >= 2) {
    for (i = 0; i < fraction_in_bump*dist.np; ++i) {
      dist.y[i] = gasdev(iran) * param4[id]/dy + ny/2.;
      if (dist.y[i] > FTYPE(ny) || dist.y[i] < 0.0) 
  	dist.y[i] = ran3(&iran) * ny;
  	//      dist.y[i] = -1.0;
    }
        for (i = i; i < dist.np; ++i) 
      dist.y[i] = ran3(&iran) * ny;
  }
  */

  for (i=0; i<dist.np; i++) dist.x[i] = ran3(&iran)*nx;
  if (ndim >= 2) {
    for (i=0; i<fraction_in_bump*dist.np; ++i) {
      dist.y[i] = gasdev(iran)*param4[id]/dy + ny/2.;
      if (dist.y[i] > FTYPE(ny) || dist.y[i] < 0.0)
	dist.y[i] = ran3(&iran)*ny;
    }
    for (i=i; i<dist.np; ++i) dist.y[i] = ran3(&iran)*ny;
  }
  
  // Add sinusoidal perturbation to peak location
  // FTYPE Hx=1.0/nx,Px=1.0,Ax=32;
  // for (i=0; i<dist.np; i++) {
  //   dist.y[i] += -1.0*Ax*sin(2*PI*Px*dist.x[i]*Hx);
  //   if (dist.y[i] < 0) dist.y[i] += ny;
  //   if (dist.y[i] > ny-1) dist.y[i] -= ny;
  // }

  // Add a sinusoidal perturbation (PPIC3D copy)
  /* -->Requires update to infile.cc
  if (denSeed_param1[id] > 0.0) {
    FTYPE Hx=1.0/(nx*nsubdomains);
    // if (mpi_rank == 0) printf("%s:%d\n\tHx = %g\n\tdenSeed_param1[%d] = %g\n\tdenSeed_param2[%d] = %g\n\n",
    // 			      __func__,__LINE__,Hx,id,denSeed_param1[id],id,denSeed_param2[id]);
    for (i=0; i<dist.np; i++) {
      dist.y[i] += -1.0*denSeed_param1[id]*sin(2*PI*
					       denSeed_param2[id]*
					       (dist.x[i]+x_start)*Hx);
      if (dist.y[i] < 0) dist.y[i] += ny;
      if (dist.y[i] > ny-1) dist.y[i] -= ny;
    }
  }
  */

  // Along x, let's have a linearly dropping distribution, as a test:
  //  pwluniform(dist.x, 0, nx, iran, linear);
  //  shuffle(dist.x, iran);

      /*
  if (ndim >= 2) {
    int prime_seed=prime(id*6+3);
    for (i = 0; i < dist.np ; ++i)  
      dist.y[i] = reverse(i+dist.np*mpi_rank, prime_seed) * ny;
  }
  
  if (ndim == 3) {
    int prime_seed=prime(id*6+4);
    for (i = 0; i < dist.np ; ++i)  
      dist.z[i] = reverse(i+dist.np*mpi_rank, prime_seed) * nz;
  }
  
  if (ndim == 1)
    gasuniform(dist.vx, dist.vx, dist.np*mpi_rank,
		       prime(id*6+5), prime(id*6+6) );
  else
    gasuniform(dist.vx, dist.vy, dist.np*mpi_rank,
		       prime(id*6+5), prime(id*6+6) );

  dist.vx *= vxthd[id] * (dt / dx);
  if (vx0d[id] !=0.0) dist.vx += vx0d[id] * (dt / dx);
  // shuffle(dist.vx,iran);

  if (ndim >=2) {
    dist.vy *= vythd[id] * (dt / dy);
    if (vy0d[id] !=0.0) dist.vy += vy0d[id]  * (dt / dy);
    // shuffle(dist.vy,iran);
  }

  if (vel_dim[id] == 3) {
    gasuniform(dist.vz, dist.vz, dist.np*mpi_rank,
		       prime(id*6+7), prime(id*6+8) );
    dist.vz *= vzthd[id] * (dt / dz);
    if (vz0d[id] !=0.0) dist.vz += vz0d[id] * (dt / dz);
    // shuffle(dist.vz,iran);
  }
  */

}

FTYPE linear(FTYPE x)
{
  
  
return (FTYPE) 1-.3*exp(-Sqr((x-nx/2.)/(nx/16.)));
}

