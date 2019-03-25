
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "eppic-types.h"
#include "eppic-mpi.h"

/* void read_domains(ArrayNd_ranged<FTYPE,2> &array,MPI_Comm across_comm,
		  char *path_name,int nghost_pts[2]) {
  int mpi_rank = 0;
  int mpi_err=MPI_Comm_rank(across_comm,&mpi_rank);

  // replace '*' in filename and path with mpi_rank 
  char buff[3];
  sprintf(buff,"%03d",mpi_rank);
  string replaceStr(buff);
  string filename(path_name);
  string::size_type pos = 0;
  while ( (pos = filename.find("*", pos)) != string::npos ) {
        filename.replace( pos, 1, replaceStr );
        pos++;
  }

  / open the file, read in and pass guard cells 
  FILE* fopensafe(char* filename, char* mode);
  char file_name[256];
  strcpy(file_name,filename.c_str());
  FILE *fp = fopensafe(file_name,"rb");
  array.bin_input_ghost(fp,nghost_pts);
  fclose(fp);

// pass the guard cells 
  void pass_guards(ArrayNd_ranged<FTYPE,2> &in, const int nx_guard[]);
  pass_guards(array, nghost_pts);

}*/


void read_domains(ArrayNd_ranged<FTYPE,NDIM> &array,MPI_Comm across_comm,
		  char *path_name,int nghost_pts[2]) {
  int mpi_rank = 0;
  int mpi_err=MPI_Comm_rank(across_comm,&mpi_rank);

  /* replace '*' in filename and path with mpi_rank */
  char buff[3];
  sprintf(buff,"%03d",mpi_rank);
  string replaceStr(buff);
  string filename(path_name);
  string::size_type pos = 0;
  while ( (pos = filename.find("*", pos)) != string::npos ) {
        filename.replace( pos, 1, replaceStr );
        pos++;
  }

  /* open the file, read in and pass guard cells */
  FILE* fopensafe(char* filename, char* mode);
  char file_name[256];
  strcpy(file_name,filename.c_str());
  FILE *fp = fopensafe(file_name,"rb");
  if (fp != NULL) cout << "file_name = " << file_name << " is open\n";
  array.bin_input_ghost(fp,nghost_pts);
  cout << "ARRAY read domain: " << array;
  fclose(fp);


  /* pass the guard cells */
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  pass_guards(array, nghost_pts);


}
