
inline PTYPE hole(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
  // This function returns a density with a constant gradient length scale:
  // n/Grad(n)
  // Adjust x for the position with respect to the entire domain 
  // decomposed system.
  // this function is only for 2d, ie is uniform in 3d; however, the uniform direction is x
  PTYPE den;
  x=x+subdomain.id_number*nx;
  PTYPE rsqr=Sqr(x-nx*nsubdomains/2)+Sqr(y-ny/2);
#if NDIM == 3 
  rsqr=Sqr(y-ny/2)+Sqr(z-nz/2);
#endif
  PTYPE gradn=1./param3[id];
  PTYPE routside_hole_sqr=param4[id]*param4[id];
  PTYPE rinside_hole_sqr=param5[id]*param5[id]; 
  PTYPE routside_hole=param4[id];
  PTYPE rinside_hole=param5[id];
  PTYPE den_min = 1./exp(gradn*(routside_hole-rinside_hole));

  if (rsqr>routside_hole_sqr) den=1.0;
  else if (rsqr<routside_hole_sqr && rsqr>rinside_hole_sqr) {
	 den=den_min*exp(gradn*(sqrt(rsqr)-rinside_hole));
  }
  else den=den_min;
  return den;
};
