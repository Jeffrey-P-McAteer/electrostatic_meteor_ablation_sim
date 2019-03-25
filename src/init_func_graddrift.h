

inline PTYPE smooth_heaviside(PTYPE x, PTYPE xnot, PTYPE k) {
  return 1.0/(1.0+exp(-2*k*(x-xnot)));
}


inline PTYPE init_func_graddrift(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
  /*
  FTYPE thisy=y-param7[id];//fmod(y,float(ny-0.4));
  FTYPE m=dy/param4[id];
    //param4[id]/float(ny-param7[id]);
  FTYPE k=param5[id];
  FTYPE y0=param6[id];
  return 1.0+m*thisy*smooth_heaviside(y,y0+param7[id],k)
    *smooth_heaviside(y,ny-y0,-1*k);
  */
  x+=subdomain.id_number*nx;
  return param3[id]*sin(2*PI*param4[id]*x)+param5[id]*sin(2*PI*param6[id]*y)+ param7[id]*fabs(y-ny/2)+ 1.0;    
}
