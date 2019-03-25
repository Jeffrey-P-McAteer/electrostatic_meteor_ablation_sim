


//    A test program for the efield_inject_tridiag routine found in the 
//    libeppic.a library
//    if compiled with debugging, the program tests the accuracy, otherwise
//    the program tests the efficiency, by running the routine multiple times
//    and reporting various timing information. 

#include "efield_inject_tridiag_test.h"


// -------------------------------------------------------------------------
// Subroutines 
// -------------------------------------------------------------------------

// -------------------
// setup mpi variables
// -------------------
void mpiSetup(int &mpi_rank, int &mpi_int,int argc, char* argv[]){

   // setup MPI
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);

}

// -------------------
// setup eppic-system
// -------------------

void setupSystem(eppic_system &sys) 
{

  sys.ndim = NDIM;
  
  sys.dx=.9;
  sys.dy=.9;
  sys.dz=.9;
  
  sys.eps=8.8542e-12;
  

#ifdef DEBUG
   sys.nx = 132;
   sys.ny = 4;
#if NDIM==2
   sys.nz = 1;
#else
   sys.nz = 16;
#endif

   if (mpi_np > 1) {
     if (mpi_np > 4) 
       nsubdomains = 4;
     else 
       nsubdomains = mpi_np;
     nsubdomains=1;
     sys.nx/=nsubdomains;
   }

#endif
   domain_decomp();
}

// --------------------------------
// deal with command line arguments
// --------------------------------
void cmdLine(int argc, char* argv[],bc_type* boundaries,eppic_system &sys,
	     int &Nit) 
{

  sys.boundary_type = 1; // assume non-periodic

 #ifdef DEBUG  
   if (argc>1) {
     if (atoi(argv[1])>-2)
       boundaries[0] = static_cast<bc_type>(atoi(argv[1]));
       if (atoi(argv[2])<3)
       boundaries[1] = static_cast<bc_type>(atoi(argv[2]));
       if (boundaries[0]==periodic) {
	 boundaries[1]==periodic;
	 sys.boundary_type = 0;
       }
   }


 #else
   if (mpi_rank==0) {
     if (argc<7) {
       cout << "Not enough arguments" <<endl;
       cout << "NX NY NZ NIT NSUBDOMAINS LHS_BOUNDARY RHS_BOUNDARY" << endl;
       MPI_Finalize();
       return 1;
     }
   }
   // initialize variables
   // SYSTEM defining parameters
   sys.nx = atoi(argv[1]);
   sys.ny = atoi(argv[2]);
   sys.nz = atoi(argv[3]);

   Nit=atoi(argv[4]); /* Number of times to run efield (for scaling) */

   nsubdomains = atoi(argv[5]);  
   boundaries[2]={static_cast<bc_type>(atoi(argv[6])),
		  static_cast<bc_type>(atoi(argv[7]))}; // natural, ie open
   if (boundaries[0]==periodic) {
     boundaries[1] = periodic;
     sys.boundary_type = 0;
   }

#endif
}

// ----------------------------
// a routine to initalize phi0 
// ----------------------------

void initPhi0(FArrayND_ranged &phi0, eppic_system sys)
{


  int phisize[NDIM] = {INDICIES(sys.nx+3,sys.ny,sys.nz)};
  int phistart[NDIM] = {INDICIES(-1,0,0)};
  
  phi0 = FArrayND_ranged(phistart,phisize);
  
  for (int iy=0;iy<sys.ny;iy++)
    for (int iz=0;iz<sys.nz;iz++) {
      FTYPE x1,xnxm1;
      x1 = sys.dx;
      xnxm1=-sys.dx;

      FTYPE slope =sys.nx*nsubdomains-1;
      for (int ix=-1;
	   ix<sys.nx+2;
	   ix++) {
	FTYPE xval = (ix+sys.nx*subdomain.id_number)/slope;
	phi0(INDICIES(ix,iy,iz)) = 
	  //	  x1
	  // 	  +slope*
   	  //val-sys.nx*nsubdomains/2)*
   	  //val-sys.nx*nsubdomains/2)
	  //  	  +sin(-10*PI*xval/(sys.nx*nsubdomains-5))
	  sin(-4*(iy+1)*PI*xval)//*(xval)
	  
	  //     	  (xval-sys.nx*nsubdomains/2)
	  //     	  +(xval)/slope
	  //+(xval)//*(xval)*(xval)
	  // 	  +slope*
	  // 	  (xval-sys.nx*nsubdomains/2)*
// 	  (xval-sys.nx*nsubdomains/2)*
// 	  (xval-sys.nx*nsubdomains/2)
	  ;
      }
//       if (subdomain.id_number == 0) {
// //  	phi0(INDICIES(0,iy,iz)) = 
// //  	  2*phi0(INDICIES(1,iy,iz))-phi0(INDICIES(2,iy,iz));
// //  	phi0(INDICIES(-1,iy,iz)) = 
// //  	  2*phi0(INDICIES(0,iy,iz))-phi0(INDICIES(1,iy,iz));
//       }
//       if (subdomain.id_number == nsubdomains-1) {
// // 	phi0(INDICIES(sys.nx,iy,iz)) = 
// // 	  2*phi0(INDICIES(sys.nx-1,iy,iz)) - phi0(INDICIES(sys.nx-2,iy,iz));
// 	phi0(INDICIES(sys.nx+1,iy,iz)) = 
// 	  2*phi0(INDICIES(sys.nx,iy,iz))-phi0(INDICIES(sys.nx-1,iy,iz));
//       }
      if (subdomain.id_number==0)
	phi0(INDICIES(-1,iy,iz)) = 0.0;
      if (subdomain.id_number==nsubdomains-1){
	phi0(INDICIES(sys.nx,iy,iz)) = phi0(INDICIES(sys.nx-1,iy,iz));
	phi0(INDICIES(sys.nx+1,iy,iz)) = phi0(INDICIES(sys.nx-1,iy,iz));
      }

    }
  FTYPE phiavg=0;
  for (int ix=-1;ix<sys.nx+2;ix++)
    for (int iy=0;iy<sys.ny;iy++)
      for (int iz=0;iz<sys.nz;iz++)
	phiavg+=phi0(INDICIES(ix,iy,iz));
  phiavg=0;
  phi0-=phiavg/(sys.nx*nsubdomains+3)/sys.ny/sys.nz;


}

// ---------------------------------------------
// a routine to return rho given a specific phi0
// ---------------------------------------------

void buildRhoFromPhi0(FArrayND_ranged phi0,FArrayND &rho,eppic_system sys)
{
  
  FTYPE dx2 = sys.dx*sys.dx;
  FTYPE dy2 = sys.dy*sys.dy;
  FTYPE dz2 = sys.dz*sys.dz;
  for(int ix=0; ix<sys.nx+1;ix++)
    for (int iy=0; iy<sys.ny; iy++)
      for (int iz=0; iz<sys.nz; iz++) {
	int iyp1=(iy+1)%sys.ny;
	int iym1=(iy-1+sys.ny)%sys.ny;
	int izp1=(iz+1)%sys.nz;
	int izm1=(iz-1+sys.nz)%sys.nz;
	rho(INDICIES(ix,iy,iz)) = 
	  -sys.eps*(
		    ((phi0(INDICIES(ix+1,iy  ,iz  ))+
		      phi0(INDICIES(ix-1,iy  ,iz  ))-
		      phi0(INDICIES(ix  ,iy  ,iz  ))*2)/dx2
		     )+
		    ((phi0(INDICIES(ix  ,iyp1,iz  ))+
		      phi0(INDICIES(ix  ,iym1,iz  ))-
		      phi0(INDICIES(ix  ,iy  ,iz  ))*2)/dy2
		     )+
		    ((phi0(INDICIES(ix  ,iy  ,izp1))+
		      phi0(INDICIES(ix  ,iy  ,izm1))-
		      phi0(INDICIES(ix  ,iy  ,iz  ))*2)/dz2
		     ));
      }

}

// -------------------------
// do all the efield solving
// -------------------------

void efield_inject_orig_noglobal(field &Efield,FArrayND &rho, int bc);

void solveForPhi(FArrayND rho, FArrayND_ranged &phi1, field &Efield) 
{

  Efield.phi_rho = 0.0;
  
  if (Efield.algorithm == 5) {
    cout << "solving with inject tridiag\n";
    efield_inject_tridiag(Efield,rho);
  } else if (Efield.algorithm == 6) {
    cout << "solving with inject orig\n";
    efield_inject_orig_noglobal(Efield,rho,1);
  } else {
    cout << "Algorithm not setup: " 
	 << Efield.algorithm 
	 << endl;
  }

  phi1 = Efield.phi_rho;

} 

// ----------------------------------------------
// test if phi solution makes sense for given rho
// ----------------------------------------------

void testPhiVsRho(FArrayND_ranged phi, FArrayND rho, eppic_system sys)
{

  FTYPE totalErr=0;
  FTYPE dx2=sys.dx*sys.dx;
  FTYPE dy2=sys.dy*sys.dy;
  FTYPE dz2=sys.dz*sys.dz;
  
  for (int ix=1;ix<sys.nx;ix++)
    for (int iy=0;iy<sys.ny;iy++)
      for (int iz=0;iz<sys.nz;iz++) {
	int iyp1=(iy+1)%sys.ny;
	int iym1=(iy-1+sys.ny)%sys.ny;
	int izp1=(iz+1)%sys.nz;
	int izm1=(iz-1+sys.nz)%sys.nz;
	FTYPE rhoxyz = -sys.eps*((phi(INDICIES(ix+1,iy  ,iz  ))+
				  phi(INDICIES(ix-1,iy  ,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dx2
				 +
				 (phi(INDICIES(ix  ,iyp1,iz  ))+
				  phi(INDICIES(ix  ,iym1,iz  ))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dy2
				 +
				 (phi(INDICIES(ix  ,iy  ,izp1))+
				  phi(INDICIES(ix  ,iy  ,izm1))-
				  phi(INDICIES(ix  ,iy  ,iz  ))*2)/dz2
				 );
	FTYPE err = abs(rho(INDICIES(ix,iy,iz))-rhoxyz);
	totalErr+=err;
      }

  if (totalErr > 1e-8)
    cout << "Average error for phi-rho comparison at domain "
	 << subdomain.id_number
	 << " is "
	 << totalErr/(sys.ny*sys.nz*(sys.nx-1)) 
	 << endl;

}

// ------------------------------
// test phi compared to given phi
// ------------------------------

void testTwoPhis(FArrayND phi0, FArrayND phi1)
{

  // note, phi0,phi1 should be FArrayND_ranged, but they are declared
  // FArrayND, so they lose their 'start' information. Scaning from 0 to nx
  // is like scaning from nxstart to nx-nxstart. 

  int nx = phi0.size(0);
  int ny = phi0.size(1);
  int nz = 1;
  if (NDIM > 2) nz = phi0.size(2);

  FTYPE totalErr=0;
  for (int ix=0;ix<nx; ix++)
    for (int iy=0;iy<ny; iy++)
      for (int iz=0;iz<nz; iz++) {
	FTYPE err = abs(phi0(INDICIES(ix,iy,iz))-phi1(INDICIES(ix,iy,iz)));
	totalErr+=err;
      }
  if (totalErr/(nx*ny*nz) > 1e-8)
    cout << "Average error for phi-phi comparison at domain "
	 << subdomain.id_number
	 << " is "
	 << totalErr/(nx*ny*nz) 
	 << endl;
  

}


// ------------------------------
// save phis and rhos
// ------------------------------

void writePhiRho(char *suffix, FArrayND_ranged phi0, FArrayND rho)
{

  FILE *Fphi0=0;
  FILE *Frho =0;
  char name[256];
  char *openbtype;
  int nghost[2] = {0,0};
  openbtype="wb\0 ";

  sprintf(name,"phi%s_%03d.bin",suffix,subdomain.id_number);
  Fphi0=fopensafe(name,openbtype);
  phi0.bin_output_ghost(Fphi0,nghost);
  //phi0.bin_output_ghost(Fphi0,nghost);
  fclose(Fphi0);

  sprintf(name,"rho%s_%03d.bin",suffix,subdomain.id_number);
  Frho=fopensafe(name,openbtype);
  rho.bin_output(Frho);
  fclose(Frho);
}


// -------------------------------------------------------------------------
// The test program
// -------------------------------------------------------------------------



int main(int argc, char* argv[]) 
{
  
  // mpi_rank, mpi_np are globals!
  mpiSetup(mpi_rank,mpi_np,argc,argv);

  eppic_system sys; 
  setupSystem(sys);
  
  int Nit=1; /* Number of times to run efield (for scaling only) */
  bc_type boundaries[2]={open,open}; 
  cmdLine(argc,argv,boundaries,sys,Nit);

  FArrayND rho;
  field Efield;
  FTYPE boundary_vals[4] = {0,0,0,0};
  init_field(Efield,rho,
	     sys.nx,sys.ny,sys.nz,
	     sys.dx,sys.dy,sys.dz,
	     sys.eps,boundaries,
	     boundary_vals,0,0,0,0,
	     // 3 = my 'cleaned' up version of marcos's code
	     // 4 = My version of Marcos's code in parallel
	     // 5 = batch parallel
	     // 6 = Marcos's original (probably not exactly)
	     5);
  FArrayND_ranged phi0,phi1;
  int nghost[2]={1,2};
  initPhi0(phi0,sys);
  buildRhoFromPhi0(phi0,rho,sys);
  testPhiVsRho(phi0,rho,sys);
  writePhiRho("old",phi0,rho);
  for (int it=0;it<Nit;it++) {
    if (mpi_rank ==0)
      cout << "Running step " << it << endl;
    solveForPhi(rho,phi1,Efield);
    testPhiVsRho(phi1,rho,sys);
    testTwoPhis(phi0,phi1);
  }
  buildRhoFromPhi0(phi1,rho,sys);
  writePhiRho("new",phi1,rho);

  MPI_Finalize();
}






// old main


//  int main(int argc, char* argv[])
//  {


//    FArrayND rho;
//    FArrayND rho_start;
//    field Efield;

//    // setup MPI

//    MPI_Errhandler mpi_err_hand;
//    void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

//    MPI_Init(&argc,&argv);
//    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
//    if (mpi_rank == 0) {
//      MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
//     }



//  #ifdef DEBUG  

//    // initialize variables
//    // system defining parameters
//    //int boundaries[2]={neumann,open}; // natural, ie open
//    bc_type boundaries[2]={periodic,periodic}; 
//    //bc_type boundaries[2]={dirichlet,dirichlet}; 
//    //bc_type boundaries[2]={open,open}; 

//    if (argc>1) {
//      if (atoi(argv[1])>-2)
//        boundaries[0] = static_cast<bc_type>(atoi(argv[1]));
//        if (atoi(argv[2])<3)
//        boundaries[1] = static_cast<bc_type>(atoi(argv[2]));
//      if (boundaries[0]==periodic) 
//        boundaries[1]==periodic;
//    }

//    int nx = 16;
//    int ny = 16;
// #if NDIM==2
//    int nz = 1;
// #else
//    int nz = 16;
// #endif

//    int Nit=1; /* Number of times to run efield (for scaling only) */


//    if (mpi_np > 1) {
//      if (mpi_np > 4) 
//        nsubdomains = 4;
//      else 
//        nsubdomains = mpi_np;
//        nx/=nsubdomains;
//    }

//  #else


//    if (mpi_rank==0) {
//      if (argc<7) {
//        cout << "Not enough arguments" <<endl;
//        cout << "NX NY NZ NIT NSUBDOMAINS LHS_BOUNDARY RHS_BOUNDARY" << endl;
//        MPI_Finalize();
//        return 1;
//      }
//    }
//    // initialize variables
//    // SYSTEM defining parameters
//    int nx = atoi(argv[1]);
//    int ny = atoi(argv[2]);
//    int nz = atoi(argv[3]);

//    int Nit=atoi(argv[4]); /* Number of times to run efield (for scaling) */

//    nsubdomains = atoi(argv[5]);  
//    bc_type boundaries[2]={static_cast<bc_type>(atoi(argv[6])),
// 			  static_cast<bc_type>(atoi(argv[7]))}; // natural, ie open

//  #endif
//    FTYPE dx=.9,dy=.9,dz=.9;
//    FTYPE eps=1;//8.8542e-12;
//    FTYPE boundary_vals[4] = {0,0,0,0};
//    init_field(Efield,rho,nx,ny,nz,dx,dy,dz,eps,boundaries,boundary_vals);
//    domain_decomp();

//    rho=0.0;
//    int seed=77;
//    FTYPE ran_num=0;
//    for(int ix=0;ix<nx*nsubdomains;ix++) {
// 	 for(int iy=0; iy<ny;iy++) {
// 	   for(int iz=0; iz<nz;iz++) {
// 	     ran_num = 2*ran3(&seed)-1;
// 	       if (ix>=subdomain.id_number*nx) 
// 		 if (ix<subdomain.id_number*nx+nx) 
// 		   rho(INDICIES(ix-nx*subdomain.id_number,iy,iz)) = ran_num;
// 	 //cos(2*PI*1*(ix+mpi_rank*nx)/FTYPE(nx*nsubdomains))*
	 
// 	 //	   	   sin(2*PI*1*iy/FTYPE(ny));
// // 	   sin(2*PI*1*iz/FTYPE(nz))+
// // 	   cos(4*PI*1*(ix+mpi_rank*nx)/FTYPE(nx*nsubdomains))*
	 
// // 	   sin(4*PI*1*iy/FTYPE(ny))*
// // 	   sin(4*PI*1*iz/FTYPE(nz));
// ;
// 	//sin(ix/FTYPE(nx))*sin(iy/FTYPE(ny));
//        }
//      }
//    }

//    //rho dc part calculation
//   FTYPE rho_dc=0.;
//   for(int ix=0;ix<nx;ix++) 
//     for(int iy=0;iy<ny;iy++)
//       for(int iz=0;iz<nz;iz++) rho_dc+=rho(INDICIES(ix,iy,iz));
//   rho_dc/=(nx*ny*nz);
//   FTYPE rho_dc_tmp=rho_dc;
//   MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
// 		subdomain.neighbor_comm);
//   rho_dc/=nsubdomains;
  
//   bool zero_dc=false;
//   zero_dc = (Efield.boundary[0]==periodic)||(Efield.boundary[1]==periodic);
//   zero_dc = (zero_dc||
// 	     ((Efield.boundary[0]==neumann)||
// 	      (Efield.boundary[1]==neumann))
// 	     );
//   // if not periodic or neumann, set rho_dc to 0, ie don't use it
//   if (!zero_dc) rho_dc=0;

//   rho_start=rho;
//   for (int it=0;it<Nit;it++){
//     rho = rho_start;
//     //get phi
//     Efield.phi_rho = 0.0;
//     efield_inject_tridiag(Efield,rho);
//     if (mpi_rank == 0) 
//       if (it==0) 
// 	first_time = efield_time;

//     /* binary write to file */

//     if (it==0) {
//       void write_domains(FArrayND_ranged &array,MPI_Comm across_comm,
// 			 char *path_name,int nghost_pts[2],
// 			 bool non_periodic_x);
      
//       FArrayND_ranged rho_ranged;

//       if (subdomain.rank==0) {
// 	int non_periodic = 1;
// 	if ((Efield.boundary[0]==periodic)||(Efield.boundary[1]==periodic))
// 	  non_periodic = 0;
// 	int nghost_phi[2]={1,2};
// 	write_domains(Efield.phi_rho,
// 		      MPI_COMM_WORLD,"efield_inject_tridiag_phi*.out",
// 		      nghost_phi,non_periodic);
	
// 	int nghost_rho[2]={0,1};
// 	//int nstart_rho[]={INDICIES(0,0,0)};
// 	//rho_ranged= FArrayND_ranged(nstart_rho,nghost_rho);
// 	rho_ranged=rho_start;
// 	write_domains(rho_ranged,
// 		      MPI_COMM_WORLD,"efield_inject_tridiag_rho*.out",
// 		      nghost_rho,0);
//       }
//     }


// #ifdef DEBUG
//    int failed=0;

//     FArrayND error_list=rho;
//     error_list=0;
//     //differentiate phi

//     int ixm1=-1;
//     int ixp1=0;
//     int iz=0;
//     FTYPE errtot=0.;
//     int ixstart=0;
//     if (boundaries[0] != periodic) {
//       ixstart=1;
//     }
//     ixm1=ixstart-1;
//     for (int ix=ixstart; ix<nx; ix++) {
//       ixp1=ix+1;

//       int iym1=ny-1;
//       for (int iy=0; iy<ny; iy++) {
// 	int izm1=nz-1;
// 	for (int iz=0; iz<nz; iz++) {
// 	FTYPE rho_test,rho_now;
	
// 	rho_now = rho_start(INDICIES(ix,iy,iz));
// 	rho_test=0;
// 	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dx));
// 	rho_test += Efield.phi_rho(INDICIES(ixp1,iy,iz))/Sqr(dx);
// 	rho_test += Efield.phi_rho(INDICIES(ixm1,iy,iz))/Sqr(dx);
// 	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dy));

// 	rho_test += Efield.phi_rho(INDICIES(ix,(iy+1)%ny,iz))/Sqr(dy);
// 	rho_test += Efield.phi_rho(INDICIES(ix,iym1,iz))/Sqr(dy);
// #if NDIM == 3
// 	rho_test += Efield.phi_rho(INDICIES(ix,iy,iz))*(-2/Sqr(dz));
// 	rho_test += Efield.phi_rho(INDICIES(ix,iy,(iz+1)%nz))/Sqr(dz);
// 	rho_test += Efield.phi_rho(INDICIES(ix,iy,izm1))/Sqr(dz);
// #endif	
// 	rho_test *= -1.0*eps;
// 	rho_test += rho_dc;

	
// 	if (fabs(rho_now) > 1e-6) {
// 	  if (fabs((rho_test-rho_now)/rho_now)>1e-6){
// 	    failed++;
// 	    error_list(INDICIES(ix,iy,iz)) = (rho_test-rho_now);
// 	  }
// 	}
// 	else {
// 	  if (fabs((rho_test-rho_now))>1e-6){
// 	    failed++;
// 	    error_list(INDICIES(ix,iy,iz)) = (rho_test-rho_now);
// 	  }
// 	}

// 	errtot +=fabs(rho_test-rho_now);

// 	izm1=iz;
// 	}
// 	iym1=iy;
//       }
//       ixm1=ix;
//     }
//     // sync failed
//     int allfailed=0;
//     MPI_Allreduce(&failed,&allfailed,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
//     if (allfailed) {
//       for (int iproc=0;iproc<mpi_np;iproc++) {
// 	if (iproc==mpi_rank) {
// 	  if (failed) {
// 	    cerr << "Error on processor " << mpi_rank << endl;
// 	    cerr << "Subdomain: " << subdomain.id_number 
// 		 << " Rank: " << subdomain.rank << endl;
// 	    cerr.flush();
// 	    char errormsg[300];
// 	    sprintf(errormsg,"%10s %4s %4s %4s %16s\n",
// 		    " ","ix","iy","iz","deviation");
// 	    cerr << errormsg;
// 	    for (int ix=0;ix<nx;ix++) {
// 	      for (int iy=0;iy<ny;iy++) {
// 		for (int iz=0;iz<nz;iz++) {
// 		  if (error_list(INDICIES(ix,iy,iz))!=0) {
// 		    sprintf(errormsg,"%10s %4d %4d %4d %16g\n",
// 			    " ",ix+subdomain.id_number*nx,iy,iz,
// 			    error_list(INDICIES(ix,iy,iz)));
// 		    cerr << errormsg;
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	MPI_Barrier(MPI_COMM_WORLD);
//       }
//       MPI_Finalize();
//       return 1;
//     }
// #endif
//   }

// #ifndef DEBUG
//   /* total time, average time, total less first, average less first */
//   if (mpi_rank == 0) cout <<  nx << "\t";
//   if (mpi_rank == 0) cout <<  ny << "\t";
//   if (mpi_rank == 0) cout <<  nz << "\t";
//   if (mpi_rank == 0) cout <<  mpi_np << "\t";
//   if (mpi_rank == 0) cout <<  nsubdomains << "\t";
//   if (mpi_rank == 0) cout <<  Nit << "\t";
//   if (mpi_rank == 0) cout <<  boundaries[0] << "\t";
//   if (mpi_rank == 0) cout <<  boundaries[1] << "\t";
//   if (mpi_rank == 0) cout <<  efield_time << "\t";
//   if (mpi_rank == 0) cout <<  efield_time/Nit<< "\t";
//   if (mpi_rank == 0) cout <<  efield_time - first_time<< "\t";
//   if (mpi_rank == 0) cout << (efield_time - first_time)/(Nit-1) << endl;
// #endif
//   MPI_Finalize();
//   return 0;
//    }


// // a routine to initialize phi to make open boundary conditions

// void build_phi0(FArrayND &phi0, 
// 		ArrayNd<FTYPE,1> phi_left,
// 		ArrayNd<FTYPE,1> phi_right
// 		)
// {
//   nx = phi0.size(0)-2;
//   ny = phi0.size(1);
//   nz = 1;
// #if NDIM==3
//   nz = phi0.size(2);
// #endif
  


// }

