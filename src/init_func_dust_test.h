
inline PTYPE dust_test_init(int id, INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
  PTYPE den;
  x=x+subdomain.id_number*nx;

  /*
  FTYPE s = 1.0/(ny*dust_sigma);
  FTYPE a = (y-ny*dust_mid)*s;
  FTYPE g = -0.5*dust_charge*dust_den*n0d[id]*exp(-0.5*a*a);
  return g*(1+sqrt(1+(n0d[id]/g)*(n0d[id]/g)));
  */

  /*
  FTYPE s = 1.0/(ny*dust_sigma);
  FTYPE a = (y-ny*dust_mid)*s;
  FTYPE g = -0.5*dust_charge*dust_den*n0d[id]*exp(-0.5*a*a);
  den = g*(1+sqrt(1+(n0d[id]/g)*(n0d[id]/g)));
  den /= n0d[id];
  */

  FTYPE s = 1.0/(ny*dust_sigma);
  FTYPE a = (y-ny*dust_mid)*s;
  FTYPE g = -1.0*dust_charge*dust_den*n0d[id]*exp(-0.5*a*a);
  den = g/n0d[id];

  return den + 1.0;
}
