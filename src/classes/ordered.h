/* This class generates a repeatable bit reversed sequence of base, b, for 
   a more uniform distribution than true random numbers would generate.  
   n should be prime */

struct ordered {

  unsigned long long int n;
  double n_inv;
  ordered(int n_in) : n(n_in), n_inv(1.0/n_in) {}

  double reverse(long long int istart) {
    int inum=istart;
    double rev=0;
    double power=1.0;
    do {
      unsigned long long int iquot = inum / n;
      unsigned long long int irem = inum - n * iquot;
      power *= n_inv;
      rev += irem * power;
      inum = iquot;
      //      cout << "\ny\t" << iquot << "\t" << irem << "\t" << power << "\t" << rev << "\n";
    } while (inum > 0);
    return rev;
  }
};

// Return the next 2 points for a uniformly filled gaussian distribution
/*inline double GaussDev(int i, double &vx,double &vy) {
  static ordered part(39);
  double vmag=part.reverse(2*i);
  double theta = (2.*M_PI) * part.reverse(2*i+1);//Problem - can't use reverse twice!
  double fac = sqrt ((-2.)*log(vmag));
  vx = fac * cos(theta);
  vy = fac * sin(theta);
  return vx;
  }*/


  
    
