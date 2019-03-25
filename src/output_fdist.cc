// Output of phase-space distribution functions performed here

#include <stdio.h> 
#include <string.h> 

#include "eppic.h"
#include "eppic-mpi.h"

void output_fdist(particle_dist *pic, int it, char *dir)
{
  extern void gather_domain_scale(ArrayNd<OTYPE, 2*NDIM> &den, PTYPEAVec* a[], 
				  FTYPE *xmin, FTYPE *xmax, int *nmesh, 
				  int *wrap, int *wrapTo, 
				  FTYPE scale, int np, PTYPEAVec xmissing);

  static FILE *ffdist[MAXDIST];
  int ifd;
  int ngrid[] = {INDICIES(nx*nsubdomains,ny,nz)};
  static int first_entry=TRUE;
  if (first_entry) {
    first_entry=FALSE;


    // Open the files
    for (ifd=0; ifd < fndist; ifd++) {
      char name[256];
      sprintf(name,"%sfdist%d.bin",dir,ifd);
      unsigned long skip=((static_cast<unsigned long>(it)-1)/
			  (nout*fnout[ifd])+1)*
	sizeof(OTYPE)*fn(ifd,0)*fn(ifd,3);
      if (ndim >= 2) skip *= fn(ifd,1)*fn(ifd,4);
      if (ndim == 3) skip *= fn(ifd,2)*fn(ifd,5);
      
      char *openbtype;
      if (it == 0) openbtype="wb\0 ";
      else openbtype="ab\0";
      
      if (hdf_output_arrays == 0) {
	FILE* fopensafe(char* filename, char* mode, unsigned long skip);
	ffdist[ifd]  = fopensafe(name,openbtype,skip);
      }
    }

    // Check than fdist requests do not go out of range 
    for (ifd=0; ifd < fndist; ifd++) {
      for (int in=0; in < fnof[ifd]; in++) {
	int id=fof[ifd][in];
	if ( id >= ndist || id < 0 )
	  terminate (0,"Error: f dist requests nonexistant distribution\n");
      }

//       // catch max/min errors:
//       int nperdim[] = {nx,ny,nz};
//       for (int idim=0;idim<3;idim++) {
// 	if (fnmin(ifd,idim) < 0) {
// 	  if(mpi_rank==0) 
// 	    cout << endl
// 		 << ";; Warning: fnmin for fdist("
// 		 << ifd
// 		 << ") and dimension "
// 		 << idim
// 		 << " is being reset to 0!" 
// 		 << endl
// 		 << ";           ";
// 	  fnmin(ifd,idim) = 0;
// 	}
// 	if (fnmax(ifd,idim) > dx*nperdim[idim]*nsubdomains) {
// 	  fnmax(ifd,idim) = dx*nperdim[idim]*nsubdomains;
// 	  if(mpi_rank==0) 
// 	    cout << endl
// 		 << ";; Warning: fnmax for fdist("
// 		 << ifd
// 		 << ") and dimension "
// 		 << idim
// 		 << " is being reset to "
// 		 << nx*nsubdomains*dx
// 		 << "!" 
// 		 << endl
// 		 << ";           ";
// 	}

//       }

    }

    
  } // end if (first_entry) 
  
  for (ifd=0; ifd < fndist; ifd++) 
    if (it%(fnout[ifd]*nout)==0) {
      
      // Define the array to sum onto
      ArrayNd<OTYPE,2*NDIM> f(INDICIES( fn(ifd,0), fn(ifd,1), fn(ifd,2) ),
			  INDICIES( fn(ifd,3), fn(ifd,4), fn(ifd,5) ));

      // Define array limits

      int gn[]={INDICIES( fn(ifd,0), fn(ifd,1), fn(ifd,2) ),
		INDICIES( fn(ifd,3), fn(ifd,4), fn(ifd,5) )};
      
      // the -subdomain.id_number*nx shifts the min,max to the frame of the
      // current domain.
      FTYPE gmin[]={INDICIES( fnmin(ifd,0)/dx-subdomain.id_number*nx, fnmin(ifd,1)/dy, fnmin(ifd,2)/dz ),
		    INDICIES( fnmin(ifd,3)*dt/dx, fnmin(ifd,4)*dt/dy, 
			      fnmin(ifd,5)*dt/dz )};
      FTYPE gmax[]={INDICIES( fnmax(ifd,0)/dx-subdomain.id_number*nx, fnmax(ifd,1)/dy, fnmax(ifd,2)/dz ),
		    INDICIES( fnmax(ifd,3)*dt/dx, fnmax(ifd,4)*dt/dy, 
			      fnmax(ifd,5)*dt/dz )};
      
      // Loop over the distributions to be summed 
      for (int in=0; in < fnof[ifd]; in++) {
	int id=fof[ifd][in];
	
	// Generate a set of pointers to vectors of particles
	PTYPEAVec *marker[NDIM*2]={INDICIES(&pic[id].x, &pic[id].y, &pic[id].z),
				   INDICIES(&pic[id].vx,&pic[id].vy,&pic[id].vz)};
	
	// the wrap in x occurs when the particle position is nx*nsubdomains 
	// regardless of domain decomposition.
	// this is because gather_scale shifts the position relative to xmin
	// this works because xmin and xmax are limited to the size of the box
	// there would be issues if the were allowed to have any value
	// there will be issues if the wrapping should be turned off, 
	// i.e. non-periodic in x runs. 
// 	int wrap[]={
// 	  INDICIES(static_cast<int>(nx*nsubdomains*fn(ifd,0)/(gmax[0]-gmin[0])),fn(ifd,1),fn(ifd,2)),
// 	  INDICIES(fn(ifd,3),fn(ifd,4),fn(ifd,5))};
	int wrap[]={
	  INDICIES(static_cast<int>((nx*nsubdomains-fnmin(ifd,0)/dx)*fn(ifd,0)/(gmax[0]-gmin[0])),static_cast<int>((ny-fnmin(ifd,1)/dy)*fn(ifd,1)/(gmax[1]-gmin[1])),static_cast<int>((-fnmin(ifd,2)/dz)*fn(ifd,2)/(gmax[2]-gmin[2]))),

	  INDICIES(fn(ifd,3),fn(ifd,4),fn(ifd,5))};
	int wrapTo[]={
	  INDICIES(static_cast<int>((-fnmin(ifd,0))/dx*fn(ifd,0)/(gmax[0]-gmin[0])),static_cast<int>((-fnmin(ifd,1))/dy*fn(ifd,1)/(gmax[1]-gmin[1])),static_cast<int>((-fnmin(ifd,2))/dz*fn(ifd,2)/(gmax[2]-gmin[2]))),
	  INDICIES(0,0,0)};

	
	FTYPE fscale=pic[id].n0avg/npd[id];
	for (int i=0; i<ndim; i++) fscale *= fn(ifd,i)*ngrid[i]/(gmax[i]-gmin[i]);
	for (int i=ndim; i<ndim*2; i++) fscale *= fn(ifd,i);
	// Gather the total distribution in x,y,z,vx,vy,vz of (marker) particles
	gather_domain_scale(f, marker, gmin, gmax, gn, wrap, wrapTo, fscale, pic[id].np, pic[id].x);
	
      }
      
#ifdef USE_MPI
      
      { // Sum the f arrays accross all processors 
	
	ArrayNd<OTYPE,2*NDIM> f_sum(INDICIES( fn(ifd,0), fn(ifd,1), fn(ifd,2) ),
				    INDICIES( fn(ifd,3), fn(ifd,4), fn(ifd,5) ));
	
	int mpi_err=MPI_Reduce( (void*) &f(0), (void*) &f_sum(0), f.length(), 
				MPI_OTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
	if ( mpi_err != MPI_SUCCESS) 
	  mpi_error(mpi_rank, mpi_err,"Reduce f call in output_fdist failed");
	f = f_sum;
	f /= mpi_np;
      }
#endif    
      
      //Print everything:
      if (mpi_rank == 0) {
	if (hdf_output_arrays == 0) fwrite(&f(0), sizeof(OTYPE), f.length(), ffdist[ifd]);
      }

    }

  return;
}

