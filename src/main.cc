 //     A 3-D Parallelized PIC Code to study electron/ion beam kinetics
//     Code written by: 
//     Meers Oppenheim 
//     Started Feb., 1999

#include <cstdio>
#include <iostream>
#include <ctime>
#include <hdf5.h>
#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-efield.h"
#include "eppic-times.h"

//#ifdef HAVE_PETSC  
#if HAVE_PETSC
#include "petsc.h"
#include <fftw3.h>

static char help[] = "EPPIC is a parallel PIC-based plasma simulator.\n\n";
static char petscOpts[] = "petsc.io";
//static char * petscOpts = 0;
#endif

#include "global_defs.h"

int main(int argc, char* argv[])
{
  particle_dist *pic;
  fluid *fspecie;
  field Efield;
  particle_misc *misc;
  FArrayND rho;
  FArrayND_ranged qden,Gx;
#if NDIM > 1
  FArrayND_ranged Gy;
#if NDIM > 2
  FArrayND_ranged Gz;
#endif
#endif
  int it0=1;
  H5_FID = -1;

  extern void infile(char *,int);
  extern void init_misc(char *dir_name_header);
  extern void init_particles(particle_dist *&pic, particle_misc *&misc);
  extern void init_fluid(fluid *&, particle_dist *&pic);

  extern void output(char *outdir, particle_dist *, fluid *, field &, 
  		     FArrayND &, FArrayND_ranged &, 
		     INDICIES(FArrayND_ranged &,
			      FArrayND_ranged &,
			      FArrayND_ranged &),
		     particle_misc *misc, int , int);
  extern void leapadv_subcycle(particle_dist *, fluid *, field &, 
			       FArrayND &, FArrayND_ranged &, 
			       INDICIES(FArrayND_ranged &,
					FArrayND_ranged &,
					FArrayND_ranged &),
			       particle_misc *, int it);
  extern void boundary(particle_dist *pic, field &Efield, 
		       particle_misc *misc, int it);
  extern void check(void);
  extern void restart(char* outdir, particle_dist *&pic, 
		      fluid *&fspecie, field &Efield, int &it0);


  // Start the parallel processor 
#ifdef USE_MPI
  int mpi_err;
  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

  MPI_Init(&argc,&argv);

  mpi_err=MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"rank call in main");
  mpi_err=MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"size call in main");
  if (mpi_rank == 0) {
    mpi_err=MPI_Comm_get_errhandler(MPI_COMM_WORLD,&mpi_err_hand);
    printf("\nMPI STARTING %d PROCESSORS\n MPI Error Handler: %d\n\n", 
	   mpi_np, mpi_err_hand);
  }
#endif

  // Get starting time
  time_t timebin;
  char timestr[128];
  time(&sim_start_time);
  if (mpi_rank==0) {
    time(&timebin);
    strftime(timestr, 128, "%A %x %T",localtime(&timebin));
    printf("EPPIC Starting at: %s\n\n", timestr);
  }

  if (mpi_rank == 0 && subdomain.rank == 0) {
    printf("\n1.0/CLOCKS_PER_SEC = %f\n",(1.0)/CLOCKS_PER_SEC);
  }
  
  // Syncronize the c++ and c print buffers:
  ios::sync_with_stdio();

  //  Read parameters from the input file: 
  infile(argv[1],1);

#if HAVE_PETSC
  /* This can be refactored to avoid creating a subcomm if
     petsc_np == mpi_np */
  MPI_Group world_group,petsc_group;
  MPI_Comm_group(MPI_COMM_WORLD,&world_group);
  if (petsc_np < 0 || petsc_np > mpi_np) petsc_np = mpi_np;
  if (mpi_rank == 0) 
    printf("\nRunning PETSc on %d of %d processors\n\n",petsc_np,mpi_np);
  int petsc_ranks[petsc_np];
  for (int i=0;i<petsc_np;i++) petsc_ranks[i] = i;
  MPI_Group_incl(world_group,petsc_np,petsc_ranks,&petsc_group);
  // MPI_Comm_create(MPI_COMM_WORLD,petsc_group,&petsc_comm);
  // PETSC_COMM_WORLD = petsc_comm;
  MPI_Comm_create(MPI_COMM_WORLD,petsc_group,&PETSC_COMM_WORLD);

  PetscErrorCode perr;
  if (mpi_rank < petsc_np) {
    perr=PetscInitialize(&argc,&argv,petscOpts,help);CHKERRQ(perr);
    // extern PetscErrorCode echo_petsc(char *);
    // perr=echo_petsc(petscOpts);CHKERRQ(perr);

    perr=KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(perr);
    perr=KSPGetPC(ksp,&pc);
    perr=PCSetFromOptions(pc);CHKERRQ(perr);
    perr=KSPSetFromOptions(ksp);CHKERRQ(perr);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Check the input parameters to make sure no obvious errors were made
  check();

  // Organize the domain decomposition communicators:
  domain_decomp();

  //  Initialize the dynamic variables : 
  init_misc(argv[1]);

#if HAVE_PETSC
  /* Initialize PETSc Mat object after defining global_length
     in init_misc() */
  if (mpi_rank < petsc_np) {
    perr=MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(perr);CHKMEMQ;
    perr=MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,
		     global_length,global_length);CHKERRQ(perr);CHKMEMQ;
    perr=MatSetType(A,MATMPIAIJ);CHKERRQ(perr);CHKMEMQ;
    perr=MatSetFromOptions(A);CHKERRQ(perr);CHKMEMQ;
    perr=MatSeqAIJSetPreallocation(A,stencil_size,NULL);CHKERRQ(perr);CHKMEMQ;
    perr=MatMPIAIJSetPreallocation(A,stencil_size,NULL,
				   stencil_size,NULL);CHKERRQ(perr);CHKMEMQ;
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  // Initialize distributions
  // -->There must be an elegant way to initialize only the necessary 
  //    number of distributions of each type (i.e. NOT 
  //    "pic=new particle_dist [ndist]" when not all the distributions
  //    are PIC). 11Jan2017 (may)
  init_particles(pic, misc);
  init_fluid(fspecie, pic);
  bc_type boundary_default[2]={bc_type(field_boundary_type[0][0]),
                               bc_type(field_boundary_type[0][1])};
  if (mpi_rank==0) {
    for (int idim=0; idim<ndim; idim++){
      if (boundary_type[idim] != PERIODIC) {
        // Changed by Glenn to deal with open boundaries
        //if (boundary_type = INJECT) {
        cout << "; running with non-periodic ";
        if (idim==0) cout << "x ";
        if (idim==1) cout << "y ";
        if (idim==2) cout << "z ";
        cout << "particle boundary conditions" << endl;
        cout << ";    LHS field boundary = ";
        switch (field_boundary_type[idim][0]) {
        case open:
          cout << "open ";
          break;
        case periodic:
          cout << "periodic ";
          break;  
        case neumann:
          cout << "neumann ";
          break;
        case dirichlet:
          cout << "dirichlet ";
          break;
        }
        cout << endl;
        cout << ";    RHS field boundary = ";
        switch (field_boundary_type[idim][1]) {
        case open:
          cout << "open ";
          break;
        case periodic:
          cout << "periodic ";
          break;  
        case neumann:
          cout << "neumann ";
          break;
        case dirichlet:
          cout << "dirichlet ";
          break;
        }
        cout << endl;
      }
    }
  }
  FTYPE boundary_vals[4] = {0,0,0,0};
  if (boundary_type[1] != PERIODIC){
    yguard_size = 1;
  }
  if (boundary_type[2] != PERIODIC){
    zguard_size = 1;
  }

  if (efield_algorithm == 2 )
    // init_field(Efield,qden,divG,nx,ny,nz,dx,dy,dz,eps,
    // 	       boundary_default,boundary_vals,
    // 	       Ex0_external,Ey0_external,Ez0_external,
    // 	       Ex0_rate,Ey0_rate,Ez0_rate,
    // 	       efield_algorithm,
    // 	       efield_zero_dc);
    init_field(Efield,qden,
	       INDICIES(Gx,Gy,Gz),nx,ny,nz,dx,dy,dz,eps,
    	       boundary_default,boundary_vals,
    	       Ex0_external,Ey0_external,Ez0_external,
    	       Ex0_rate,Ey0_rate,Ez0_rate,
    	       efield_algorithm,
    	       efield_zero_dc);
  else
    init_field(Efield,rho,nx,ny,nz,dx,dy,dz,eps,
	       boundary_default,boundary_vals,
	       Ex0_external,Ey0_external,Ez0_external,
	       Ex0_rate,Ey0_rate,Ez0_rate,
	       efield_algorithm,
	       efield_zero_dc);

  // Restart from dump file if requested 
  if (iread != 0) {
    restart(outdir, pic, fspecie, Efield, it0);
  }

  //  Calculate the charges and currents on the grid. 
  if (efield_algorithm == 2)
    // quasineutral_den_flux(qden, divG, pic, fspecie, 0);
    quasineutral_den_flux(qden, INDICIES(Gx,Gy,Gz), pic, fspecie, 0);
  else
    charges(rho, pic, fspecie, 0);

  //  Find the electric field on the grid at t=n:
  // efield_wrapper(Efield, rho, qden, divG, 0);
  efield_wrapper(Efield, rho, qden, INDICIES(Gx,Gy,Gz), 0);

  // output any initial diagnostics: 
  // if (iread == 0) 
  //   output (argv[1], pic, fspecie, Efield, rho, qden, divG, misc, it, it0);
  if (iread == 0) 
    output (argv[1], pic, fspecie, Efield, rho, qden, 
	    INDICIES(Gx,Gy,Gz), misc, it, it0);

  if(mpi_rank==0) cout << "; --- All initialization complete: starting iterations" << endl;

  // Main timestep loop: 
  for (it = it0; it <= nt; it++) {

    // Apply the standard leapfrog method 
    // leapadv_subcycle(pic, fspecie, Efield, rho, qden, divG, 
    // 		     misc, it);
    leapadv_subcycle(pic, fspecie, Efield, rho, qden, 
		     INDICIES(Gx,Gy,Gz), 
		     misc, it);

    // Deal with any Boundary condition issues
    boundary(pic, Efield, misc, it);

    // Output data, diagnostics and restart: 
    if (it%nout == 0) {
      // DEBUGGING
      /*if (subdomain.rank==0) {
        for (int id=0; id < ndist; id++)
          {
            int vpart = 0; 
            PTYPE vx = 0;
            PTYPE vy = 0;
            PTYPE vz = 0;
            PTYPE vvar = 0;;
 
            for (int ip=0; ip < pic[id].np; ip++){
              if (pic[id].x[ip] > fabsent2){
                vpart++;
                vx += pic[id].vx[ip]*dx/dt;
                vy += pic[id].vy[ip]*dy/dt;
                vz += pic[id].vz[ip]*dz/dt;
                vvar += pow(pic[id].vx[ip]*dx/dt,2) +
                  pow(pic[id].vy[ip]*dy/dt,2) +
                  pow(pic[id].vz[ip]*dz/dt,2);
              }
            }
            if (vpart > 0){
              printf("Subdomain %d, Dist %d Mean Vel: %f, Variance: %f, Num Part: %d\n", 
                     subdomain.id_number, id, 
                     sqrt(vx*vx+vy*vy+vz*vz)/vpart,
                     vvar/vpart, vpart);
            }else{
              printf("Subdomain %d, Dist %d, Num Part: %d (No paricltes)\n", 
                     subdomain.id_number, id,  vpart);
                     }
          }
      }
      for (int id=0; id < ndist; id++){
        printf("Subdomain %d, Dist %d, np %d, nabsent %d\n", 
               subdomain.id_number, id, pic[id].np, pic[id].nabsent);
      }*/

      // Check to see if an externally signaled error is set:
      if (end_asap) nt=it;
      // output (argv[1], pic, fspecie, Efield, rho, qden, divG, misc, it, it0);
      output (argv[1], pic, fspecie, Efield, rho, qden, 
	      INDICIES(Gx,Gy,Gz), misc, it, it0);
    }


  }//  End of main timestep loop 

  if (mpi_rank==0) {
    time(&timebin);
    strftime(timestr, 128, "%A %x %T",localtime(&timebin)); //"
    printf("EPPIC ending normally at: %s\n\n", timestr);
  }
  
  /* Note: PetscFinalize will call MPI_Finalize if PetscInitialize 
     was called before MPI_Init. In that case, terminate() needs
     to know to not call MPI_Finalize when PETSc is enabled.
     20Jun2013 (may)
  */

  //#ifdef HAVE_PETSC
#if HAVE_PETSC
  if (mpi_rank < petsc_np) {
    perr=MatDestroy(&A);CHKERRQ(perr);
    perr=KSPDestroy(&ksp);CHKERRQ(perr);
    perr=PetscFinalize();CHKERRQ(perr);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  terminate(end_asap,"\n");

  return 0;
} // End of Program 


  void finish_asap(int sig)
{
  if (mpi_rank==0) {
    time_t timebin;
    char timestr[128];
    time(&timebin);
    strftime(timestr, 128, "%A %x %T",localtime(&timebin));
    fprintf(stderr, "Caught signal %d: at %s\n  %s\n",sig,timestr,
	    "Will terminate after next output");
  }
  //    end_asap = TRUE;
}
