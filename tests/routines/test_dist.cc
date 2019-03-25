// Program to test my distribution creating programs.


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;
#include "Random.h"
#include "ordered.h"
#include "eppic-types.h"

int mpi_rank=0;
//#include "../create_shell_particle.cc"

int main() {

  //Test random Normal Deviate routine:
  cout << "Test random Normal Deviate routine: \n";
  GaussDev test_dist(133);
  for (int i=0; i<10; i++) {
    cout << test_dist.dev()<< endl;
  }

  //Test uniform distribution creation
  cout << "Test uniform distribution creation \n";
  FTYPE reverse(int num, int n);

  ordered two(2), five(5);

  for (int i=0; i<100; i++) {
    cout <<i << "\t" << reverse(i,2) <<"\t" << two.reverse(i) <<"\t" 
	 << reverse(i,5) <<"\t" << five.reverse(i) 
	 << endl;
  }

  //Test random ring generation
  ofstream ring_file;
  ring_file.open("ring.dat");
  PTYPE vradius=10., vth=1., vx, vy, vz;
  void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz);
  for (int i=0; i<1000000; i++) {
    create_shell_particle(vradius, vth, vx, vy, vz);
    ring_file << vx << "\t" << vy  << "\t" << vz << endl;
  }
  ring_file.close();
  exit(0);
}

FTYPE reverse(int num, int n)
{

  /* Local variables */
  int irem, inum;
    FTYPE power;
    int iquot;
    FTYPE rev;
    
    
    /*     Function to reverse the digits of num in base n. */
    
    rev = 0.;
    inum = num;
    power = 1.;

    do {
      iquot = inum / n;
      irem = inum - n * iquot;
	power /= n;
	rev += irem * power;
	inum = iquot;
	//	cout << "\nx\t" << iquot << "\t" << irem << "\t" << power << "\t" << rev << "\n";
    } while(inum > 0);
    
    return rev;
} /* End reverse */
