
inline PTYPE gaussian_axis(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z))
{
// U is a unit vector defining the direction of column
// A is a vector point from the "origin", ie the middle of the box to x,y,z
// B is the vector from the nearest center of the column to x,y,z
// Param3-5 = Ux,Uy,Uz (need to be normalized ahead of time)
// Param6 = height of column maximum in units of N0
// Param7 = variance (in m)
//          this gets converted in init_particles into the use below, ie
//          param7[id]=Sqr(dx/param7[id]);
// Param8 = background height, typically 1 or 0 (default)
  x=x+subdomain.id_number*nx;
  PTYPE Ux=param3[id],Uy=param4[id],Uz=param5[id];
  PTYPE den;  
  PTYPE Ax=x-PTYPE(nx*nsubdomains/2),Ay=0,Az=0;
#if NDIM == 2
  Ay=y-PTYPE(ny/2);
#elif NDIM == 3
  Ay=y-PTYPE(ny/2);
  Az=z-PTYPE(nz/2);
#endif
//   static int first_entry=0;
//   if (mpi_rank == 0) {
//     if (first_entry < 100) {
//       cout << "Entry: " << first_entry;
//       cout << " X = " << x;
//       cout << " y = " << y;
//       cout << " z = " << z;
//       cout << " AX = " << Ax;
//       cout << " Ay = " << Ay;
//       cout << " Az = " << Az;
//       cout << " UX = " << Ux;
//       cout << " Uy = " << Uy;
//       cout << " Uz = " << Uz << "\n";
//       first_entry++;
//     }
//   }
  PTYPE L=Ax*Ux+Ay*Uy+Az*Uz;
  PTYPE Bx=Ax-L*Ux,By=Ay-L*Uy,Bz=Az-L*Uz;
  PTYPE distsqr=Sqr(Bx)+Sqr(By)+Sqr(Bz);

  if (start_col[id] == 0.0 && stop_col[id] == 1.0) den = param6[id]*exp(-distsqr*param7[id])+param8[id]; 
  else {
    // Finite trail: cut off meteor at the edges of the box (i.e. restrict the column to between 
    // x=startc and x=stopc) - LKT 12/1/16
    // Update: now viable for a meteor tilted in any direction - LKT 06/28/17

    //calculate length of desired meteor trail & distance of current location from center
    float crit_angle = atan((float)nz/(float)(nx*nsubdomains));
    float met_angle = atan((float)Uz/(float)Ux);
    int box_dist;
    if (met_angle > crit_angle) {
      float len = nz/2;
      box_dist = 2*len/sin(met_angle);
    } else if (met_angle < crit_angle) {
      float len = nx*nsubdomains/2;
      box_dist = 2*len/cos(met_angle);
    } else box_dist = 2*sqrt((nx*nsubdomains/2)*(nx*nsubdomains/2) + (nz/2)*(nz/2));
    float met_length = (stop_col[id]-start_col[id])*box_dist;
    float dist_from_origin = sqrt(Bx*Bx + By*By + Bz*Bz);

    /*
    // Exponential fall off at ends:
    //    use base distribution without adding in background, overlay exponential fall off,
    //    add background
    den = param6[id]*exp(-distsqr*param7[id]);
    if (dist_from_origin > met_length/2.) {
      float frac=exp(-(dist_from_origin-met_length/2.)/met_length/2.);
      den = den*frac + param8[id];
    } else {
      den = den + param8[id];
      }*/

    // just for x-direction meteors
    // convert start_col & stop_col fractions to grid locations
    int startc,stopc;
    startc = start_col[id]*nsubdomains*nx;
    stopc = stop_col[id]*nsubdomains*nx;
    // coefficients of 4/startc and 4/stopc are arbitrary - a higher coefficient will make the ends        fall off more quickly
    den = param6[id]*exp(-distsqr*param7[id]);
    if (x < startc) {
      den = den*exp(4*(x-startc)/startc) + param8[id];
    } else if (x > stopc) {
      den = den*exp(-4*(x-stopc)/startc) + param8[id];
    } else {
      den = den + param8[id];
    }
    
  }
  
  return den;
};

