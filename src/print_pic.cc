

#include "eppic.h"
#include "eppic-mpi.h"
#include <iostream>
#include <fstream>

void print_pic(char *fname, particle_dist &pic_dist)
{
  
  filebuf fb;
  fb.open(fname,ios::out);
  ostream fout(&fb);

  char line[160];

  for (int i=0;i<pic_dist.np;i++) {
    if (pic_dist.x(i)>-1) {
      sprintf(line,"%4d\t%.6e\t%.6e\t%.6e\t%.6e\0",mpi_rank,pic_dist.x(i)+subdomain.id_number*nx,pic_dist.y(i),pic_dist.vx(i),pic_dist.vy(i));
      fout << line << endl;
    }
  }
  fb.close();
}
