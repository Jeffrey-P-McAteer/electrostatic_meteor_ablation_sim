
#include "ArrayNd.h"

#define NDIM 2

typedef ArrayNd<double,NDIM> FArrayND;
typedef ArrayNd<int,NDIM> IArrayND;

main(int argc, char* argv[])
{
  cout << "....... Starting ... arraynd_test\n";

  /* construction */
  FArrayND arr1;  
  arr1 = FArrayND(5,3);
  FArrayND arr2 = FArrayND(5,3);
  FArrayND arr3;
  FArrayND arr4; /* use: formated i/o test */
  FArrayND arr5 = FArrayND(5,3); /* use: binary i/o test */

  /* assignment */
  arr1 = -3.0;
  /* assignment: lhs access */
  arr1(1,1) = 1.0;
  arr1(2,2) = 2.0;

  /* copy */
  arr3=arr1;

  /* rhs access */
  cout << "....... For arr1, the value of 2,2 is " << arr1(2,2) << "\n";
  cout << "....... For arr2, the value of 2,2 is " << arr2(2,2) << "\n";
  cout << "....... For arr3, the value of 2,2 is " << arr3(2,2) << "\n";

  /* arithmetic */

  /* other functions */

  /* output */
  arr1.output("arr1data.out");
  FILE *fp;

  if ((fp = fopen("binary_arr1data.out","wb")) != NULL) {
    arr1.bin_output(fp);
    fclose(fp);
  }

  /* construction: read in */
  arr4.input("arr1data.out");
  if ((fp = fopen("binary_arr1data.out","rb")) != NULL) {
    arr5 = 10;
    arr5.bin_input(fp);
    fclose(fp);
  }

  /* C++ style output */
  cout << "....... arr5: \n" << arr5;

  /* C++ style output to file */
  arr4.output("arr4data.out");
  arr5.output("arr5data.out");

  /* Test destroy method using  mkArray() */
  void mkArray();
  for (int i=0;i<5;i++) {
    mkArray();
  }

  cout << "....... Ending ... arraynd_test\n";
}


void mkArray()
{
  /* make an array and see if destroy method works */
  FArrayND arr5 = FArrayND(5,3); 
  arr5 = 5;
}
