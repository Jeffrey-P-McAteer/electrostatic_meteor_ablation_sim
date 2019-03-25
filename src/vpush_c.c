/* The following is a clone of vpush written in straight C - for the purposes 
   of evaluating the efficiency of the c compiler. */
/* Advance the particle velocities to the next time step. */

extern double dt;     /* time step in seconds */
extern int nt;        /* Number of timesteps to evaluate */
extern double dx, dy; /* grid step size in meters */
extern int nx, ny;    /* number of grid cells */

void vpush_c(float *adv_vx, float *adv_vy, float *old_vx, float *old_vy, 
	     float *cur_x, float *cur_y, 
	     double **Efield_Ex, double **Efield_Ey, 
	     double qmrat, double alpha, int np)
{
  double Exscale, Exiscale;
  double Eyscale, Eyiscale;
  int i, ix, iy;
  int ixl, ixh;
  int iyl, iyh;
  double wxh;
  double wyh;
  double Exll, Exlh, Exhl, Exhh, Exl, Exh, Exi;
  double Eyll, Eylh, Eyhl, Eyhh, Eyl, Eyh, Eyi;

  /* Rescale the electric field so that it may simply be added to
  particle velocities in the particle loop. */
  Exscale =  alpha * qmrat * (dt*dt/dx) ;
  Exiscale = 1. / Exscale;
  Eyscale =  alpha * qmrat * (dt*dt/dy) ;
  Eyiscale = 1. / Eyscale;
  for (ix = 0; ix<nx; ix++) 
    for (iy = 0; iy<ny; iy++) {
      Efield_Ex[ix][iy] *= Exscale;
      Efield_Ey[ix][iy] *= Eyscale;
    }
  
 
  /* For each particle, update the velocity ... */
  for (i = 0; i < np; ++i) {

    /* Define the nearest grid points: */
    ixl = (int) cur_x[i];
    ixh = ixl + 1;
    if (ixh == nx) ixh = 0;

    iyl = (int) cur_y[i];
    iyh = iyl + 1;
    if (iyh == ny) iyh = 0;

    /* And the corresponding linear weighting factors: */
    wxh = cur_x[i] - ixl;
    wyh = cur_y[i] - iyl;

    /* Calculate Exi, the linearly-weighted Ex at the particle: */
    Exll = Efield_Ex[ixl][iyl];
    Exlh = Efield_Ex[ixl][iyh];
    Exhl = Efield_Ex[ixh][iyl];
    Exhh = Efield_Ex[ixh][iyh];
    Exl = Exll + wyh*(Exlh-Exll);
    Exh = Exhl + wyh*(Exhh-Exhl);
    Exi = Exl + wxh*(Exh-Exl);

    /* Calculate Ezi, the linearly-weighted Ez at the particle: */
    Eyll = Efield_Ey[ixl][iyl];
    Eylh = Efield_Ey[ixl][iyh];
    Eyhl = Efield_Ey[ixh][iyl];
    Eyhh = Efield_Ey[ixh][iyh];
    Eyl = Eyll + wyh*(Eylh-Eyll);
    Eyh = Eyhl + wyh*(Eyhh-Eyhl);
    Eyi = Eyl + wxh*(Eyh-Eyl);

    /* Advance the particle velocities: */
    adv_vx[i] = old_vx[i] + Exi;
    adv_vy[i] = old_vy[i] + Eyi;

    } /* End for */

  /*  Restore original values to the electric field array: */
  
  for (ix = 0; ix<nx; ix++) 
    for (iy = 0; iy<ny; iy++) {
      Efield_Ex[ix][iy] *= Exiscale;
      Efield_Ey[ix][iy] *= Eyiscale;
    }

    return ;
} /* vpush */
