/*
  Contains two functions for initializing sinusoidal density distributions.
*/

inline PTYPE sine_init(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
  PTYPE den;
  x=x+subdomain.id_number*nx;
#if NDIM == 1
  den = param3[id]*sin(param4[id]*x)
    + 1.0;    
#endif 
#if NDIM == 2
  den = param3[id]*sin(2*PI*param4[id]*x)+param5[id]*sin(2*PI*param6[id]*y)
    + 1.0;    
#endif 
#if NDIM == 3
  den = param3[id]*(sin(param4[id]*x)*sin(param5[id]*y)*sin(param6[id]*z))
    + 1.0;    
#endif 

  return den;
};

inline PTYPE cosine_init(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
  PTYPE den;
  x=x+subdomain.id_number*nx;
#if NDIM == 1
  den = param3[id]*cos(param4[id]*x)
    + 1.0;    
#endif 
#if NDIM == 2
  den = param3[id]*cos(2*PI*param4[id]*x)+param5[id]*cos(2*PI*param6[id]*y)
    + 1.0;    
#endif 
#if NDIM == 3
  den = param3[id]*(cos(param4[id]*x)*cos(param5[id]*y)*cos(param6[id]*z))
    + 1.0;    
#endif 

  return den;
};

