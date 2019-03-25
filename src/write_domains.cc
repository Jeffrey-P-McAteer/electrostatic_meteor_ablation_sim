
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "eppic-types.h"
#include "eppic-mpi.h"

void write_domains(ArrayNd_ranged<FTYPE,2> &array,MPI_Comm across_comm,
		   char *path_name,int nghost_pts[2],bool non_periodic_x) {
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
  /* open the file and write */
  FILE* fopensafe(char* filename, char* mode);
  char file_name[256];
  strcpy(file_name,filename.c_str());
  FILE *fp = fopensafe(file_name,"wb");
  int nghost_pts_tmp[2]={0,0};
  nghost_pts_tmp[0]=nghost_pts[0];
  nghost_pts_tmp[1]=nghost_pts[1];
  if (non_periodic_x) {
    if (subdomain.id_number==0) 
      nghost_pts_tmp[0]=0;
    if (subdomain.id_number==nsubdomains-1)
      nghost_pts_tmp[1]=0;
  }

  array.bin_output_ghost(fp,nghost_pts_tmp);
  fclose(fp);

}

void write_domains(ArrayNd_ranged<FTYPE,3> &array,MPI_Comm across_comm,
		  char *path_name,int nghost_pts[2],bool non_periodic_x) {
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
  /* open the file and write */
  FILE* fopensafe(char* filename, char* mode);
  char file_name[256];
  strcpy(file_name,filename.c_str());
  FILE *fp = fopensafe(file_name,"wb");
  int nghost_pts_tmp[2]={0,0};
  nghost_pts_tmp[0]=nghost_pts[0];
  nghost_pts_tmp[1]=nghost_pts[1];
  if (non_periodic_x) {
    if (subdomain.id_number==0) 
      nghost_pts_tmp[0]=0;
    if (subdomain.id_number==nsubdomains-1)
      nghost_pts_tmp[1]=0;
  }

  array.bin_output_ghost(fp,nghost_pts_tmp);
  fclose(fp);

}
