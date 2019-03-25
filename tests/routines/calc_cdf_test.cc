
#include "eppic-misc.h"
#include "eppic-math.h"

int main(int argc, char* argv[])
{

  int nx=200;
  FTYPEAVec F=FTYPEAVec(nx);

  F(0) = 0;
  FTYPE ftotal=0;
  for (int ix=1;ix<nx;ix++) {
    F(ix)=sin(3.14*10*ix/(nx+1.0))+1;
  }
  FTYPEAVec cdf;
  calc_cdf(F,cdf);
  FTYPEAVec Inv_cdf;
  int ninv = nx*10;
  invert_function(cdf,Inv_cdf,ninv);

  char buffer[80];
  sprintf(buffer,"calc_cdf_test_F.dat\0");
  std::ofstream gnuplot_out(buffer);
  for (int i=0;i<nx;i++) {
    sprintf(buffer,"%8d\t%16.8e\n",i,F(i));
    gnuplot_out << buffer;
  }
  gnuplot_out.close();

  sprintf(buffer,"calc_cdf_test_cdf.dat\0");
  gnuplot_out.open(buffer);
  for (int i=0;i<nx;i++) {
    sprintf(buffer,"%8d\t%16.8e\n",i,cdf(i));
    gnuplot_out << buffer;
  }
  gnuplot_out.close();

  sprintf(buffer,"calc_cdf_test_Inv_cdf.dat\0");
  gnuplot_out.open(buffer);
  for (int i=0;i<ninv;i++) {
    sprintf(buffer,"%16.8e\t%16.8e\n",i*cdf(nx-1)/(ninv-1),Inv_cdf(i));
    gnuplot_out << buffer;
  } 
  gnuplot_out.close();


}
