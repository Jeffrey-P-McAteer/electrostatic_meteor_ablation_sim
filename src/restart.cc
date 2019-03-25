/* Initialize those arrays necessary for a restart*/
#include "eppic.h"
#include "eppic-mpi.h"
#include <stdio.h> 
#include <string.h> 
#include <cmath> 

FILE* fopensafe(char* filename, char* mode);
FILE* fopensafe(char* filename, char* mode, unsigned long skip);

// intermediate array; needs to be external for restart purposes
extern fluid *f_old;

void restart(char* outdir, particle_dist *&pic, fluid *&fspecie, field &Efield, int &it0)
{
  /* Restart from restart.bin file */
  
  
  FILE *fdump;
  int id, len;
  char dir[128],name[128], *rest;
  
  /*Extract the directory */
  len=0;
  if (restart_nonlocal) {
    rest=strrchr(outdir,'\0');
    if (rest != NULL ) {
      len=(int) (rest-outdir+1);
      strncpy(dir,outdir,len);
      //strncpy(dir,outdir);
      strcat(dir,"restart/");
      len+=8;
    }
  } else {
    strncpy(dir,"restart/\0",len=8);
  }
  dir[len]='\0';

  // setup some fluid stuff
  void restart_fluid(fluid &fspecie, FILE *fdump, int id);
  void bcast_fluid (fluid &fspecie, int id);

  if (method.min() < 0) {
    if (iread == 0) {
      try {
	f_old = new fluid [ndist];
      } catch (const std::bad_alloc &e) {
	std::cout << "Allocation failed (fluid): " << e.what() << "\n";
	MPI_Abort(MPI_COMM_WORLD,-1);
      }
    }
  }

  // this loop insures that only nfiles_parallel processors 
  //inputs at a time. Something neccesary for the BU machines.
  int ndumps= static_cast<int>(ceil(FTYPE(mpi_np/nfiles_parallel)));
  if (ndumps==0) ndumps=1;

  for (int ipro=0; ipro<ndumps; ++ipro) {
    if (ipro == mpi_rank%ndumps) {
      sprintf(name,"%srestart%d.rst",dir,mpi_rank);
      fdump = fopensafe(name,"rb");
      fread(&it0,sizeof(it0),1,fdump);
      for (id=0; id<ndist; ++id) {
	
	if (method[id] >= 0) {
	  fread(&(pic[id].n0avg),sizeof(pic[id].n0avg),1,fdump);
	  fread(&(pic[id].np),sizeof(pic[id].np),1,fdump);
	  if (pic[id].np >= pic[id].x.size()) {
	    
	    terminate(-1,"Error:  Restart exceeds the number of particles\n");
	  }
	  fread(&(pic[id].x[0]),sizeof(pic[id].x[0]),pic[id].np,fdump);
	  if (ndim >= 2)
	    fread(&(pic[id].y[0]),sizeof(pic[id].y[0]),pic[id].np,fdump);
	  fread(&(pic[id].vx[0]),sizeof(pic[id].vx[0]),pic[id].np,fdump);
	  if (ndim >= 2)
	    fread(&(pic[id].vy[0]),sizeof(pic[id].vy[0]),pic[id].np,fdump);
	  if (vel_dim[id] == 3) 
	    fread(&(pic[id].vz[0]),sizeof(pic[id].vz[0]),pic[id].np,fdump);
	  if (ndim==3) // Read this last to make the program compatible with df2d
	    fread(&(pic[id].z[0]),sizeof(pic[id].z[0]),pic[id].np,fdump);
	  
	  // Fill the absent array.
	  pic[id].nabsent=0;
	  for (int i=0; i<pic[id].np; i++) {
	    if (pic[id].x(i) < -nx) {
	      // This particle is absent
	      if (pic[id].nabsent >= pic[id].absent.size()) 
		terminate(-1,"Error: Array absent in xpush_domain too small: increase part_pad");
	      pic[id].absent(pic[id].nabsent++)=i;
	    }
	  }
	// } else { // method[id] < 0, ie fluid
	} else if (method[id] < 0) {
	  if (subdomain.rank == subdomain.root){
	    if (method[id] != -4) {
	      restart_fluid(fspecie[id],fdump,id);
	      f_old[id].den = fspecie[id].den; 
	      f_old[id].vx = fspecie[id].vx;
	      f_old[id].vy = fspecie[id].vy;
	      f_old[id].vz = fspecie[id].vz;
	      restart_fluid(f_old[id],fdump,id);
	    }
	  }
	}
      }
      if (efield_algorithm == 2) {
	fread(&(Efield.phi_rho(INDICIES(-phix_guard_size[0],0,0))),
	      sizeof(Efield.phi_rho(INDICIES(-phix_guard_size[0],0,0))),
	      Efield.phi_rho.size(),fdump);      
      }
      fclose(fdump);
    }

    // The following avoids reading all files at once:
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  for (id=0; id<ndist; ++id) {
    // if (method[id] < 0) {	
    //   if (efield_algorithm!=2) {
    if (method[id] < 0 && method[id] != -4) {
	
      bcast_fluid(fspecie[id],id);
      if (subdomain.rank != subdomain.root){
	f_old[id].den = fspecie[id].den; 
	f_old[id].vx = fspecie[id].vx;
	f_old[id].vy = fspecie[id].vy;
	f_old[id].vz = fspecie[id].vz;
      }
	
      bcast_fluid(f_old[id],id);
    }
  }



  //#ifdef USE_MPI
  //MPI_Barrier(MPI_COMM_WORLD);
  //#endif

      


//   for (id=0; id<ndist; ++id) {
//     if (method[id] < 0) {

//       printf("Passing fluid (%d) info within subdomain. ME:%d\n",id,mpi_rank);
//       MPI_Barrier(MPI_COMM_WORLD);
//       bcast_fluid(fspecie[id],id);
      
//       f_old[id].den = fspecie[id].den; 
//       f_old[id].vx = fspecie[id].vx;
//       f_old[id].vy = fspecie[id].vy;
//       f_old[id].vz = fspecie[id].vz;


//       //#ifdef USE_MPI
//       // MPI_Barrier(MPI_COMM_WORLD);
//       //#endif
//       printf("Now passing f_old. ME:%d\n",mpi_rank);
//       MPI_Barrier(MPI_COMM_WORLD);
//       bcast_fluid(f_old[id],id);

//     }
//   }
  if (mpi_rank == 0) printf("Restarting at t=%6d \n\n",it0);
  /* Set the time step to start at the next iteration: */
  it0++;
    

  
#ifdef USE_DOMAINS
  MPI_Barrier(MPI_COMM_WORLD);
  void pass_fluid_guard(fluid&,int id);

	  // initialize guard cells
  for (id=0; id<ndist; ++id) {
    // if (method[id] < 0) {
    if (method[id] < 0 && method[id] != -4) {

//       if (subdomain.rank == 0) { // broadcast
	
//       } else { // receive
// 	    restart_fluid(fspecie[id],fdump,id);
// 	    f_old[id].den = fspecie[id].den; 
// 	    f_old[id].vx = fspecie[id].vx;
// 	    f_old[id].vy = fspecie[id].vy;
// 	    f_old[id].vz = fspecie[id].vz;
//       }
//       if (subdomain.rank == 0) { // broadcast

//       } else { // receive
// 	restart_fluid(f_old[id],fdump,id);
//       }
      
      pass_fluid_guard(fspecie[id],id);
    }
  }
#endif
  
  
}
