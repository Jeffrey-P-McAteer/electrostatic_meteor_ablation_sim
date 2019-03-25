// We do not want MPI C++ bindings
#define MPI_NO_CPPBIND
#define MPI_SKIP_MPICXX
#define MPICH_IGNORE_CXX_SEEK
#include <mpi.h>

#include "eppic-types.h"


int mpi_rank=0;
int mpi_np=1;
int proc_west,proc_east;

#if NDIM == 2

main(int argc, char* argv[])
{
  cout << "....... Starting ... arraynd_test\n";
  /* setup MPI */
  MPI_Errhandler mpi_err_hand;

  /* setup Petsc */
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if (mpi_rank == 0) {
    MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
  }




  /* construction */
  int xrange=6;
  int yrange=4;
  if (mpi_np>1) {
    xrange/=mpi_np;
    proc_east=(mpi_rank + 1)%mpi_np;
    proc_west=(mpi_rank + mpi_np-1)%mpi_np;
  }


  FArrayND arr1;  
  FArrayND_ranged arr1_ranged; /* for copying from arr1 */
  const int range[]={xrange+2,yrange};
  int start[]={-1,0};
  /* for single assign */
  FArrayND_ranged arr2_ranged=FArrayND_ranged(start,range); 
  int start_and_range[2][2]={{-1,xrange+2},{-1,yrange}};
  /* for assign every point */
  FArrayND_ranged arr3_ranged=FArrayND_ranged(start_and_range);
  /* for reading in */
  FArrayND_ranged arr4_ranged; //=FArrayND_ranged(start,range);
  FArrayND_ranged arr5_ranged; //=FArrayND_ranged(start,range);


  arr1 = FArrayND(xrange+2,yrange);
  
  /* copy */
  /* does not work*/
  /* arr1_ranged = arr1; */
  arr1_ranged = arr2_ranged;
  
  /* assignment */
  arr3_ranged = 4;

  /* assignment: lhs access */
  arr2_ranged(-1,0) = -1024;
  arr1_ranged(1,0) = 1024;

  /* rhs access */
  //  if (mpi_rank==0) cout << "....... For arr2_ranged, the value of -1,0 is " 
  //       << arr2_ranged(-1,0) << "\n";


  /* arithmetic */

  /* other functions */

  /* output */
  /* C++ style output */
  //  if (mpi_rank==0) cout << arr2_ranged;

  /* C++ style output to file */
  arr2_ranged.output("arr2_rangeddata.out");

  /* binary write to file */
  FILE *fp;
  if ((fp = fopen("binary_arr2_rangeddata.out","wb")) != NULL) {
    arr2_ranged.bin_output(fp);
    fclose(fp);
  }
  
  /* construction: read in */
  arr4_ranged.input("arr2_rangeddata.out");
  //  if (mpi_rank==0) cout << "....... For arr4_ranged, the value of -1,0 is " 
  //       << arr4_ranged(-1,0) << "\n";
  if ((fp = fopen("binary_arr2_rangeddata.out","rb")) != NULL) {
    arr5_ranged = arr2_ranged;
    arr5_ranged = 10;
    //    if (mpi_rank==0) cout << "....... BEFORE: For arr5_ranged, the value of -1,0 is " 
    //       << arr5_ranged(-1,0) << "\n";
    arr5_ranged.bin_input(fp);
    //    if (mpi_rank==0) cout << "....... AFTER: For arr5_ranged, the value of -1,0 is " 
    //       << arr5_ranged(-1,0) << "\n";
    fclose(fp);
  }

  /* if multiple processors, read/write domain decomposed data */

  if (mpi_np>1) {
    void read_domains(ArrayNd_ranged<FTYPE,NDIM> &array,MPI_Comm across_comm,
		    char *path_name,int *nghost_pts=0);

    /* for domain domposed input */
    int ghost_cells[]={1,1};
    FArrayND_ranged arr6_ranged=FArrayND_ranged(start,range);
    arr6_ranged = 4;
    read_domains(arr6_ranged,MPI_COMM_WORLD,
		 "./arraynd_ranged_test*.dat",ghost_cells);
    // 
    void write_domains(FArrayND_ranged &array,MPI_Comm across_comm,
		       char *path_name,int nghost_pts[2]);

    
    for (int ix=0; ix<xrange; ix++) {
      if (mpi_rank == 0) {
	for (int iy=0; iy<yrange; iy++) {
	  //cout << arr6_ranged(ix,iy)<<"\t";
	}
	flush(cout);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 1) {
	for (int iy=0; iy<yrange; iy++) {
	  //cout << arr6_ranged(ix,iy)<<"\t";
	}
	flush(cout);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      if (mpi_rank==0) {cout << endl;	flush(cout);}
      MPI_Barrier(MPI_COMM_WORLD);
    }


    //write_domains(arr6_ranged,MPI_COMM_WORLD,"./arraynd_ranged_test*.out",ghost_cells);

  }


  /* Test destroy method using  mkArray() */
  void mkArray();
  for (int i=0;i<5;i++) {
    mkArray();
  }

  /* Test larger arrays */
  /* this is a work in progress, so commented out for now */
  /*
  if (mpi_np>1) {
    void read_domains(ArrayNd_ranged<FTYPE,3> &array,MPI_Comm across_comm,
		      char *path_name,int *nghost_pts=0);

    void write_domains(ArrayNd_ranged<FTYPE,3> &array,MPI_Comm across_comm,
		       char *path_name,int nghost_pts[2]);

    ArrayNd_ranged<FTYPE,3> vel_dist;
    const int range2[]={xrange+2,yrange,yrange};
    int start2[]={-1,0,0};
    
    vel_dist = ArrayNd_ranged<FTYPE,3>(start2,range2);
    
    int ghost_cells[]={1,1};
    
    if (mpi_rank==0) cout << "about to read vel_dist\n";
    read_domains(vel_dist,MPI_COMM_WORLD,
		 "./vel_dist*.dat",ghost_cells);


    write_domains(vel_dist,MPI_COMM_WORLD,
		 "./vel_dist*.out",ghost_cells);



    vel_dist(xrange,2,2) = -1024;
    if (mpi_rank == 0) cout << vel_dist;
  }
  */


  MPI_Finalize();
  cout << "....... Ending ... arraynd_test\n";
}

void mkArray()
{
  /* make an array and see if destroy method works */
  const int range[]={5,3};
  int start[]={-1,0};
  FArrayND_ranged arr2_ranged=FArrayND_ranged(start,range);
  
}


#else


main(int argc, char* argv[]) {
  cout << "3D test not implemented!" << endl;
  return 1;
}

#endif
