
#include "eppic-misc.h"
#include "eppic-math.h"

int main(int argc, char* argv[])
{

  int nx=200;
  FTYPEAVec F=FTYPEAVec(nx);

  F(0) = 0;
  for (int ix=1;ix<nx;ix++) {
    F(ix)=sqrt(ix);
  }
  FTYPEAVec Inv_F;
  int ninv = 8;
  invert_function(F,Inv_F,ninv);

  char buffer[80];
  sprintf(buffer,"invert_function_test_F.dat\0");
  std::ofstream gnuplot_out(buffer);
  for (int i=0;i<nx;i++) {
    sprintf(buffer,"%8d\t%16.8e\n",i,F(i));
    gnuplot_out << buffer;
  }
  gnuplot_out.close();

  sprintf(buffer,"invert_function_test_Inv_F.dat\0");
  gnuplot_out.open(buffer);
  for (int i=0;i<ninv;i++) {
    sprintf(buffer,"%16.8e\t%16.8e\n",i*F(nx-1)/(ninv-1),Inv_F(i));
    gnuplot_out << buffer;
  } 
  gnuplot_out.close();


}
