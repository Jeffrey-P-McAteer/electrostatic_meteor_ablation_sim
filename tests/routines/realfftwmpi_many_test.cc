

// a test program of the realfftwmpi_many class


// things that I need to check:
// Can I change the ntransform counter from being the fastest to the slowest
// counter?
// Does the answer to the above question change for an MPI call?
// I think I can do it if I set istride to 2, idist =1;

// For the 1D case, I need to set things up using the same structure as the
// 2/3D cases, but I will manually do the "many" transforms.


#include <iostream>

#include "realfftwmpi_many.h"
#include "eppic-math.h"
#include "eppic-mpi.h"

int mpi_rank,mpi_np;
int main(int argc, char* argv[])
{

  //   realfftwmpi<FTYPE,1> oneD;
//   realfftwmpi<FTYPE,3> threeD;


  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if (mpi_rank == 0) {
    MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
  }


  int nx=32,ny=32,nz=32,ntransforms=4;
  int size_oneD[]={nx};
  int size_twoD[]={nx,ny};
  int size_threeD[]={nx,ny,nz};

  
  realfftwmpi_many<FTYPE,2> twoD=
    realfftwmpi_many<FTYPE,2>(MPI_COMM_WORLD,size_twoD,ntransforms);
  
  realfftwmpi_many<FTYPE,3> threeD=
    realfftwmpi_many<FTYPE,3>(MPI_COMM_WORLD,size_threeD,ntransforms);

  realfftwmpi_many<FTYPE,1> oneD=
    realfftwmpi_many<FTYPE,1>(MPI_COMM_WORLD,size_oneD,ntransforms);

  int x;

  fftw_real *data, *work;
  fftw_complex *cdata;
  
  /// non_mpi
  //int fftw_plan_type = FFTW_ESTIMATE|FFTW_IN_PLACE;
  cout << mpi_rank << ": not broken at 1\n";
  //  rfftw_plan plan_nompi, iplan_nompi;
//   plan_nompi=
//     rfftw_create_plan(nx,FFTW_FORWARD,fftw_plan_type);
//   iplan_nompi=
//     rfftw_create_plan(nx,FFTW_BACKWARD,fftw_plan_type);
//   rfftwnd_plan plan_nompi, iplan_nompi;
//   plan_nompi=
//     rfftwnd_create_plan(1,&nx,FFTW_FORWARD,fftw_plan_type);

//   iplan_nompi=
//     rfftwnd_create_plan(1,&nx,FFTW_BACKWARD,fftw_plan_type);

  data =new fftw_real[2*(nx/2+1)*ntransforms];
  work =new fftw_real[2*(nx/2+1)*ntransforms];

  int realsize[2] ={2*(nx/2+1),ntransforms};
  ArrayNd<fftw_real,2> oneDdata;
  oneDdata.build_ArrayNd(realsize,0,data);
  int complexsize[2] = {nx/2+1,ntransforms};
  cdata = (fftw_complex *)(data);


  fftw_complex czero;
  czero.re = 0.;
  czero.im = 0.;
  ArrayNd<fftw_complex,2> oneDcdata; 
  oneDcdata.build_ArrayNd(complexsize,czero,cdata);
//   for (int it=0;it<ntransforms;it++)
//     for (x=0;x<nx;x++) 
//       //      oneDdata(x,it) = sin(2*PI*x/nx)*(it+1);
//       oneD(x,it) = sin(2*PI*x/nx)*(it+1);
  //data[x] = sin(2*PI*x/nx);
  cout << mpi_rank << ": not broken at 2\n";
//   rfftw(plan_nompi, 1,
// 	data, 1, 0,
// 	work , 1, 0);
  //rfftw_one(plan_nompi,data,	work);

  if (mpi_rank==0) {
    ofstream quickout;
    quickout.open("1dtest_before.dat");
    for (x=0;x<nx;x++) 
      //quickout << x << "\t" << oneDdata(x,0) << "\t" << oneDdata(x,1) << endl;
      quickout << x << "\t" << oneD.data(x,0) << "\t" << oneD.data(x,1) << endl;
  }


  //  for (int itransform=0;itransform<ntransforms;itransform++)
//     rfftwnd_real_to_complex(plan_nompi,ntransforms,
// 			    data,//+2*(nx/2+1)*itransform,
// 			    ntransforms,1
// 			    ,NULL,
// 			    //,oneDcdata.address(0),//+(nx/2+1)*itransform,
// 			    ntransforms,1);
			    
  oneD.transform();

  //  *data = *work;
			    
  
  cout << mpi_rank << ": not broken at 3\n";
  oneD.cdata(2,0)=complex<FTYPE>(oneD.cdata(2,0).real(),1);
  oneD.cdata(1,1)=complex<FTYPE>(oneD.cdata(2,0).real(),1);
  //work[4]=1;
  //  cdata[3].re = 1;
  //cdata[2].re = 1;
  //cdata[nx/2+1+2].re = 1;
  //oneDcdata(1,0).re=1;
  //cdata[1].re = 1;
  //oneDcdata(2,1).re=2;
  //data[1] = 1.0;


  if (mpi_rank==0) {
    ofstream quickout;
    quickout.open("1dtest_mid.dat");
    for (x=0;x<nx/2+1;x++) {
      quickout << x; 
      for (int itransform=0;itransform<ntransforms;itransform++)
// 	quickout << "\t" << oneDcdata(x,itransform).re
// 		 << "\t" << oneDcdata(x,itransform).im 
// 		 << "\t" << oneDdata(x,itransform);
	quickout << "\t" << oneD.cdata(x,itransform).real()
		 << "\t" << oneD.cdata(x,itransform).imag()
		 << "\t" << oneD.data(x,itransform);
      quickout << endl;
      
    }
  }


  cout << mpi_rank << ": not broken at 4\n";
	
  //  for (int itransform=0;itransform<ntransforms;itransform++)
//     rfftwnd_complex_to_real(iplan_nompi,ntransforms,
// 			    oneDcdata.address(0),//+(nx/2+1)*itransform,
// 			    ntransforms,1,
// 			    NULL,
// 			    //data,//+2*(nx/2+1)*itransform,
// 			    ntransforms,1);

  oneD.invtransform();
  //   rfftw(iplan_nompi, 1,
// 	data, 1,0,
// 	work, 1,0);

  //*data = *work;

  cout << mpi_rank << ": not broken at 5\n";


  if (mpi_rank==0) {
    ofstream quickout;
    quickout.open("1dtest_after.dat");
    for (x=0;x<nx;x++) 
//       quickout << x 
// 	       << "\t" << oneDdata(x,0) << "\t" << oneDdata(x,1) 
// 	       << "\t" << data[x*ntransforms] 
// 	       << "\t" << data[x*ntransforms+1] 
      quickout << x 
	       << "\t" << oneD.data(x,0) << "\t" << oneD.data(x,1) 
	       << endl;
  }

  cout << mpi_rank << ": not broken at 6\n";
  delete [] data;
  delete [] work;

  cout << mpi_rank << ": not broken at 7\n";
//   int itransform=0;
//   for (itransform=0;itransform<ntransforms;itransform++) 
//     for (x = 0; x < twoD.local_nx; ++x)
//       for (y = 0; y < ny; ++y) 
// 	twoD.data(x,y,itransform) 
// 	  = cos(2*PI*(x+twoD.x_start)/nx)*cos(2*PI*y/ny)*
// 	  (itransform+1);
//   for (itransform=0;itransform<ntransforms;itransform++) 
//     for (z=0; z < nz; ++z) 
//       for (y = 0; y < ny; ++y) 
// 	for (x = 0; x < threeD.local_nx; ++x){
// 	  threeD.data(x,y,z,itransform) 
// 	    = cos(2*PI*(x+threeD.x_start)/nx)*cos(2*PI*y/ny)*
// 	    cos(2*PI*z/nz)*
// 	    (itransform+1);
// 	}


//   char filename[256];

//   void output_vtkdata(char*filename,ArrayNd_ranged<FTYPE,4> &data, 
// 		      int itransform, int* index, int xstart);
//   int index[3] = {threeD.local_nx,ny,nz};

//   sprintf(filename,"before_%d_%d.vts\0",0,mpi_rank);
//   output_vtkdata(filename,threeD.data,0,index,threeD.x_start);
//   sprintf(filename,"before_%d_%d.vts\0",1,mpi_rank);
//   output_vtkdata(filename,threeD.data,1,index,threeD.x_start);
		 
//   cout << mpi_rank << ": not broken at 1\n";


//   twoD.transform();
//   threeD.transform();



//   cout << mpi_rank << ": not broken at 2\n";
//   twoD.invtransform();
//   threeD.invtransform();


//   for (itransform=0;itransform<ntransforms;itransform++) 
//     for (z=0; z < nz; ++z) 
//       for (y = 0; y < ny; ++y) 
// 	for (x = 0; x < threeD.local_nx; ++x)
// 	  threeD.data(x,y,z,itransform) /=(nx*ny*nz);

//   cout << mpi_rank << ": not broken at 4\n";
    
//   sprintf(filename,"after_%d_%d.vts\0",0,mpi_rank);
//   output_vtkdata(filename,threeD.data,0,index,threeD.x_start);
//   sprintf(filename,"after_%d_%d.vts\0",1,mpi_rank);
//   output_vtkdata(filename,threeD.data,1,index,threeD.x_start);

  // 2D case testing ordering...

  
//   rfftwnd_mpi_plan plan, iplan;
//   plan = rfftw2d_mpi_create_plan(MPI_COMM_WORLD, nx, ny,
// 				 FFTW_REAL_TO_COMPLEX, fftw_plan_type);
  
//   iplan = rfftw2d_mpi_create_plan(MPI_COMM_WORLD, nx, ny,
// 				  FFTW_COMPLEX_TO_REAL, fftw_plan_type);

//   ArrayNd_ranged<fftw_real,3> datampi;
//   ArrayNd_ranged<std::complex<fftw_real>,3> cdatampi;


//   int local_nx,local_x_start,local_ny_after_transpose,
//     local_y_start_after_transpose,total_local_size;

//   rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
//                             &local_ny_after_transpose,
//                             &local_y_start_after_transpose,
//                             &total_local_size);


//   //// ----- COPIED FROM CLASS ----- //////
//   /// need to clean up and get working for my test....
  
//   index_type nsize[]={ntransforms,nx,2*(ny/2+1)};
//   index_type memory_after_0=local_nx;

//   // Set the starting index for the array:
//   int  nstart[]={0,0,0};

//   //Our array supplies sufficient memory:
//   int local_size=memory_after_0*ntransforms;
//   datampi=ArrayNd_ranged<fftw_real,3>(nstart,nsize,x);

	
//   // transposed array
//   // Generally {0,0,0}
//   //Dimension of complex array
//   int index_zero[]={0,0,0};
//   fftw_real *start_add=datampi.address(index_zero);
//   index_type nstart_transpose[]={0,0,0};
//   index_type ncsize_transpose[]={ntransforms,local_ny_after_transpose,nx/2+1};
//   cdatampi.build_ArrayNd_ranged(nstart_transpose, ncsize_transpose, 
// 			     std::complex<fftw_real>(0.,0.), 
// 			     (std::complex<fftw_real> *)(start_add));



//   //Fourier Transform Array:
//   fftw_real *workmpi=new fftw_real[total_local_size];
	      
//   /* Now, compute the forward transform: */
//   start_add=datampi.address(index_zero);
//   rfftwnd_mpi(plan, ntransforms, start_add, workmpi, FFTW_TRANSPOSED_ORDER);
  
//   //Inverse Fourier Transform Array:
//   rfftwnd_mpi(iplan, ntransforms,  data.address(index_zero), work, 
// 		    FFTW_TRANSPOSED_ORDER);
//       else {
// 	// fix me!
// 	cerr << "1D many not implemented yet!!" << endl;
//       }

//       ordering=NOT_TRANSPOSED;

//       if (work_defined) {
// 	  delete [] work;
// 	  nwork=0;
//       }


//       return nwork; // Return the size of the work array that remains
    
//     }



  MPI_Finalize();

}


void output_vtkdata(char*filename,ArrayNd_ranged<FTYPE,4> &data, int itransform
		    , int* index, int xstart) {
  ofstream fout;
  fout.open(filename);

  int nx = index[0];
  int ny = index[1];
  int nz = index[2];
  cout << nx 
       << " - " << ny 
       << " - " << nz
       << endl;

  fout << "# vtk DataFile Version 1.4.2" << endl;
  fout << "Data created by eppic" << endl;
  fout << "ASCII" << endl;
  fout << "DATASET STRUCTURED_POINTS" << endl;
  fout << "DIMENSIONS " 
       << nx << " "
       << ny << " "
       << nz << endl;
  fout << "ORIGIN " << xstart << " 0 0" << endl;
  fout << "SPACING 1 1 1" << endl;
  fout << "POINT_DATA " << nx*ny*nz << endl;
  fout << "SCALARS data double 1" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  int newline_count=0;

  cout << "Sizes: " 
       << data.size(0) << " "
       << data.size(1) << " "
       << data.size(2) << " "
       << data.size(3) << endl;
  cout << "Starts: " 
       << data.xstart() << " "
       << data.ystart() << " "
       << data.zstart() << " "
       << endl;


  for (int z=0;z<nz;z++) 
    for (int y =0;y<ny;y++)
      for (int x =0;x<nx;x++) {
	//fout << data(x,y,z,itransform) << " ";
	fout << x + y +z ;
	newline_count++;
	if(newline_count>8) {
	  fout << endl;
	  newline_count=0;
	}
      }

  fout.close();
}
