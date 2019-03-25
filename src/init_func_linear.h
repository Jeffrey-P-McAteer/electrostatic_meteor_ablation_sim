inline PTYPE linear_gradient_init(int id,INDICIES(PTYPE x,PTYPE y,PTYPE z))
{ 
  x += subdomain.id_number*nx;
  /* return param3[id]*x + param8[id]; */
  PTYPE slope=(n0rhsd[id]-n0lhsd[id])/nx/nsubdomains;
  return slope*x + n0lhsd[id];
}
