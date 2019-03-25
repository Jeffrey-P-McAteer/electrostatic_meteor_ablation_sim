// Initialize the paticle quantities 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-system.h"
#include "eppic-misc.h"

PTYPE gauss_scale,gauss_width,den_not;
int curr_id;


// include "den_func"'s stored in separate files
#include "init_func_hole.h"
#include "init_func_gaussian.h"
#include "init_func_sine.h"
#include "init_func_graddrift.h"
#include "init_func_graddrift_inject.h"
#include "init_func_linear.h"
#include "init_func_dust_test.h"

/** A routine to initialize particle distributions. 
    \todo Add optional calculation of standard deviation from ideal distributions, use this figure as a metric for particle noise. 
*/

void init_particles(particle_dist *&pic, particle_misc *&misc)
{

    //Local variables: 
    int id, i, npad;

    // pad size of arrays in pic, since they can grow due to
    // particle migration for multiple sub domain cases


    // Change the random number seed for each processor 
    static int iran = -33 -cc_rand - mpi_rank;
    FTYPE gasdev(int &);
    FTYPE ran3(int *);
    int nxpad=0, nxshift=0;
    int nypad=0, nyshift=0;
    int nzpad=0, nzshift=0;

    int prime(int);
    FTYPE gasdev(int &, int );
    FTYPE reverse(long long int,int);

    // use enum to give different initializations memorable names
    enum {
      injection = -1,
      rand_flat = 0,
      prime_rev = 1,
      prime_rev_shell = 2,
      gaussian = 3,
      spacial_hole = 4,
      gaussian_three = 6,
      sine = 7,
      meteor = 8,
      graddrift = 9,
      external_den = 10,
      ablating_meteor = 11,
      linear_gradient = 12,
      dust_test = 13
    };
  
    extern FTYPE gdist(FTYPE , FTYPE , int );
    extern FTYPE gxdist_normal(FTYPE , int );
    extern FTYPE gydist_normal(FTYPE , int );
    extern FTYPE gzdist_normal(FTYPE , int );
    extern FTYPE gdist_normal(FTYPE, FTYPE, int );
    extern FTYPE randist(FTYPE (*func)(FTYPE, int), int id, int &);
    extern void randist(FTYPE (*func)(FTYPE, FTYPE, int),  
			FTYPE &, FTYPE &, int, int &);
    extern FTYPE f0dist(FTYPE , FTYPE , int );
    void rejectND(particle_dist &, int, int &, 
		  PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
    
    if (method.max() >= 0) {
	// Create a set of particles, one for each distribution 
	// We have ndist distributions of charged particles 
      try {
	pic=new particle_dist [ndist];
      } catch (const std::bad_alloc &e) {
	std::cout << "Allocation failed (pic): " << e.what() << "\n";
	MPI_Abort(MPI_COMM_WORLD,-1);
      }
      try {	
	misc=new particle_misc [ndist];
      } catch (const std::bad_alloc &e) {
	std::cout << "Allocation failed (misc): " << e.what() << "\n";
	MPI_Abort(MPI_COMM_WORLD,-1);
      }
    }

    // For each particle type load particle characteristics: 
    for (id=0; id<ndist; ++id) {      
	// padded number of particles
  	npad = (int)(npd[id]*part_pad[id]);
        if (npd[id] == 0){
          npad = (int) part_pad[id];
        }

	if (method[id] >= 0) {
	  // Load particle characteristics 
          pic[id].m=md[id];
          pic[id].q=qd[id];
          pic[id].n0avg=n0d[id];
          pic[id].np=npd[id];
          pic[id].np_all = npd[id]*nsubdomains;
          
          // Injection: If this is the leftmost domain then we must 
          // include enough particles for the extra cell on the left:
          
          if (boundary_type[0] !=PERIODIC && field_boundary_type[0][0] != periodic
              && subdomain.id_number==0) {
            FTYPE scale = (nx+1);
            scale /= nx;
            pic[id].np = static_cast<int>(npd[id]*scale + 0.5);
            nxpad = 1;
            nxshift = -1;
          }

          /*
          if (boundary_type[1] != PERIODIC){
            FTYPE scale = (ny+1);
            scale /= ny;
            pic[id].np = static_cast<int>(npd[id]*scale + 0.5);
            nypad = 1;
            nyshift = -1;
          }
          
          if (boundary_type[2] != PERIODIC){
            FTYPE scale = (nz+1);
            scale /= nz;
            pic[id].np = static_cast<int>(npd[id]*scale + 0.5);
            nzpad = 1;
            nzshift = -1;
          } 
          */

          // Initiailize particle position arrays:
          pic[id].x=PTYPEAVec(npad) = 0.;
          if (ndim >= 2) pic[id].y=PTYPEAVec(npad) = 0.; 
          if (ndim == 3) pic[id].z=PTYPEAVec(npad) = 0.;
          
          // Allocate the absent array to keep track of particles
          // that have migrated out of the domain
          if (nsubdomains > 1 || boundary_type[0] == INJECT) {
            // Deal with situations when np=0
            if (pic[id].np == 0 ){
              pic[id].absent=intAVec((int) part_pad[id])=-1;
            }
            else{
              pic[id].absent=intAVec((int) (pic[id].np*(part_pad[id]-1.0)))=-1;
            }
            pic[id].nabsent=0;
          }
          // Initiailize particle velocity arrays:
          pic[id].vx=PTYPEAVec(npad) = 0.;
          if (ndim >= 2) pic[id].vy=PTYPEAVec(npad) = 0.;
          if (ndim == 3 || vel_dim[id] == 3) pic[id].vz=PTYPEAVec(npad) = 0.;
          
          if (iread == 0) {	      
            // Set the number of particles:
            long long int np_bulk=pic[id].np;
            // Adjust for existence of another subpopulation
            if (n0_shell[id] > 0.0) {
              np_bulk=static_cast<int>(pic[id].np*(pic[id].n0avg-n0_shell[id])/pic[id].n0avg);
              //np_shell=pic[id].np-np_bulk;
            }
            // random distribution of particles positions: 
            
            // Load particle positions - normalized to dx 
            
            // The various initial loading distributions are represented 
            // by init_dist
            if (init_dist[id] == rand_flat) {
              //if (boundary_type[0] != PERIODIC){
              for (i = 0; i < np_bulk ; ++i) pic[id].x[i] = ran3(&iran) * (nx+nxpad)+nxshift;
              //              else{
              //for (i = 0; i < np_bulk ; ++i) pic[id].x[i] = ran3(&iran) * nx;
              //}
              // changed to if else above to deal with open boundaries
              //if (boundary_type[0] == PERIODIC) for (i = 0; i < np_bulk ; ++i) pic[id].x[i] = ran3(&iran) * nx;
              if (ndim >= 2){
                for (i = 0; i < np_bulk ; ++i)  pic[id].y[i] = ran3(&iran) * (ny+nypad)+nyshift;
              }
              if (ndim == 3) {
                //pic[id].z = fabsent;
                for (i = 0; i < np_bulk ; ++i)  pic[id].z[i] = ran3(&iran) * (nz+nzpad)+nzshift;
              }
              
            } else if (init_dist[id] == prime_rev || 
                       init_dist[id] == prime_rev_shell) {
              
              PTYPE smallest_float_below(PTYPE x);
              static const PTYPE fnx_eps=1.0+smallest_float_below(-1.0);
              
              // This mimics the particle placement for systems with identicle size and 
              // number of particles but no subdomains
              void gasuniform(PTYPEAVec &, PTYPEAVec &, int, int, int);
              long long int ic=0;
              long long int ipart;
              
              // Bit reversed distribution of particles positions: 
              // We need prime numbers seeds - one for each 
              int prime_seed=prime(id*6+2);
              int prime_seed3=prime(id*6+3);
              int prime_seed4=prime(id*6+4);
              int prime_seed5=prime(id*6+5);
              int prime_seed6= prime(id*6+6);
              
              int prime_seed7=prime(id*6+7);
              int prime_seed8= prime(id*6+8);
              
              long long int np_all = npd[id]*nsubdomains;
              
              
              /*
                if ((boundary_type[0] == INJECT)) {
                if (subdomain.id_number == 0) {
                np_all += pic[id].np - npd[id];
                } else
                np_all += 
                static_cast<int>(npd[id]*(nx*nsubdomains+nsubdomains)/(nx*nsubdomains))-npd[id];
                nxshift = -1;
                nxpad = 1;
                }
              */
              
              pic[id].np_all = np_all;
              np_bulk=static_cast<int>(pic[id].np_all*pic[id].n0avg/(pic[id].n0avg+n0_shell[id]));
              for (long long int i = 0; i < np_all ; ++i) {
                ipart=i+np_all*subdomain.rank+1;
                PTYPE xtmp= reverse(ipart, prime_seed) *(nx+nxpad+nx*(nsubdomains-1))+nxshift+fnx_eps;
                if ((xtmp >= nx*subdomain.id_number  || 
                     (boundary_type[0]==INJECT && subdomain.id_number==0) )
                    && xtmp < nx*(subdomain.id_number+1) ){
                  xtmp -= nx*subdomain.id_number;
                  pic[id].x[ic] = xtmp;
                  
                  if (ndim >= 2) {
                    pic[id].y[ic] = reverse(ipart, prime_seed3) * ny;
                  }
                  
                  if (ndim == 3) {
                    pic[id].z[ic] = reverse(ipart, prime_seed4) * nz;
                  }
                  
                  double vmag, theta, fac;
                  vmag=reverse(ipart,prime_seed5);
                  if (vmag==0) {
                    cout << "Error: " << vmag << " " << ipart << " " << prime_seed5 << "\n";
                    terminate(1," Error in gasuniform 1: infinite velocity");
                  }
                  if (i <  np_bulk) {// This part handles the bulk distribution:
                    theta = (2.*M_PI) * reverse(ipart,prime_seed6);
                    fac = sqrt ((-2.)*log(vmag));
                    if (init_dist[id]==prime_rev) { // Standard Gaussian Uniform distribution
                      pic[id].vx[ic] = (fac * cos(theta))*vxthd[id]*dt/dx;
                      if (vel_dim[id] >= 2) pic[id].vy[ic] = (fac * sin(theta))*vythd[id]*dt/dy;
                      
                      pic[id].vx[ic] += vx0d[id] * (dt / dx);
                      if (vel_dim[id] >= 2) pic[id].vy[ic] += vy0d[id] * (dt / dy);
                      
                      if (vel_dim[id] == 3) {
                        vmag=reverse(ipart,prime_seed7);
                        if (vmag==0) terminate(1," Error in gasuniform 2: infinite velocity");
                        theta = (2.*M_PI) * reverse(ipart,prime_seed8);
                        fac = sqrt ((-2.)*log(vmag));
                        pic[id].vz[ic] = (fac * sin(theta))*vzthd[id]*dt/dz;
                      }
                    } else { 
                      // init_dist 2 generates a shell in phase space where 
                      // vth is the thermal width but 
                      // vx0d is the avg. speed in phase space
                      PTYPE vr= ( (fac * cos(theta))*vxthd[id] + vx0d[id]) *dt/dx;
                      // We need 2 more angles to rotate the distribution around in 3D
                      do { // Begin rejection loop until we find a place for this particle
                        PTYPE theta2 = (M_PI) * ran3(&iran);
                        PTYPE phi = (2.*M_PI) * ran3(&iran);
                        pic[id].vx[ic] = vr * sin(theta2)*cos(phi);
                        pic[id].vy[ic] = vr * sin(theta2)*sin(phi)*dx/dy;
                        pic[id].vz[ic] = vr * cos(theta2)*dx/dz;
                        // This creates a gaussian ring around (vx,vy) but then 
                        // rotates it around the vz axis, causing a greatly enhanzed 
                        // density on this axis.  We can repair this by rejecting 
                        // particles in proportion to (vr^2-vz^2)/vr^2
                      } while (ran3(&iran) > sqrt((Sqr(vr)-Sqr(pic[id].vz[ic]))/Sqr(vr)) );
                      
                    }
                    
                  } else { // This part adds shell particles:
                    void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz);
                    create_shell_particle(v0_shell[id]*dt/dx,vth_shell[id]*dt/dx,
                                          pic[id].vx(ic),pic[id].vy(ic),pic[id].vz(ic));
                    pic[id].vy(ic) *= dx/dy;
                    pic[id].vz(ic) *= dx/dz;
                  }
                  
                  /*
                    if (id==0) 
                    printf("%2d %4d %4d, %16.8f, %16.8f, %16.8f, %16.8f\n",
                    
                    mpi_rank,ic,ipart,
                    pic[id].x[ic],pic[id].y[ic],
                    pic[id].vx[ic],pic[id].vy[ic]);
                  */
                  ic++;
                  
                } else {
                  // cout << "Warning: rank " << mpi_rank << " particle " << ic << 
                  // " outside: x=" << xtmp << endl;		
                }
              }
              
              pic[id].np=ic;
              
              //#endif
            } else if (init_dist[id] == gaussian) {
              // // Place particles according to the subroutine position 
              // // This is for nonuniform distribution in space 
              void position(particle_dist &, int, int & );
              position(pic[id], id, iran);
            } else if (init_dist[id] == spacial_hole) {
              // Place particles according to the subroutine rejectND
              // This is for nonuniform distribution in space 
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              // YANN EDIT:
              if (mpi_rank==0)
                cout << ";; Warning: Parameter definitions have changed" 
                     << " as of 1/25/11. If this input was generated"
                     << " before this date, scale each value by"
                     << " 1/(2*sqrt(2))."
                     << endl;

#if NDIM == 2
              if (param4[id]==0) param4[id]=nx*nsubdomains/4;
              if (param4[id]>nx*nsubdomains/2) param4[id]=nx*nsubdomains/2;
#else 
              // In 3D make a cylindrical hole aligned along x
              if (param4[id]==0) param4[id]=nz/4.; 
              if (param4[id]>nz/2) param4[id]=nz/2;
#endif
              if (param4[id]>ny) param4[id]=ny;
              if (param5[id]==0) param5[id]=param4[id]/8.;
              if (param5[id]>=param4[id]) param5[id]=param4[id]/8.;
              
              if (param3[id]==0) param3[id]=(param4[id]-param5[id])/2.;
              
              if (mpi_rank==0) {
                cout << "Hole parameters for id= "<< id 
                     << " param3 = " << param3[id] 
                     <<" param4 = " << param4[id] 
                     << " param5 = " << param5[id] <<" den_min = "
                     << 1./exp(1./param3[id]*(param4[id]-param5[id])) <<endl;
              } 
              
              rejectND(pic[id], id, iran, hole);
              // Place positions randomly:
              for (i = 0; i < pic[id].np; ++i) 
                pic[id].vx[i] = (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
              if (ndim >= 2) for (i = 0; i < pic[id].np; ++i) 
                               pic[id].vy[i] = (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
              if (vel_dim[id] == 3) 
                for (i = 0; i < pic[id].np; ++i) 
                  pic[id].vz[i] = (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
            } else if (init_dist[id] == gaussian_three) {
              // Place particles according to the subroutine rejectND
              // This is for nonuniform distribution in space 
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              param7[id]=Sqr(dx/param7[id]);
              //		    gauss_scale=1.0;
              //		    gauss_width=0.05;
              rejectND(pic[id], id, iran, gaussian_axis);
              
              // Place positions randomly:
              for (i = 0; i < pic[id].np; ++i) 
                pic[id].vx[i] = (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
              if (ndim >= 2) for (i = 0; i < pic[id].np; ++i) 
                               pic[id].vy[i] = (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
              if (vel_dim[id] == 3) 
                for (i = 0; i < pic[id].np; ++i) 
                  pic[id].vz[i] = (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
              
            } else if (init_dist[id] == meteor) {
              // Place particles according to the subroutine rejectND
              // This is for nonuniform distribution in space 
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              param7[id]=Sqr(dx/param7[id]);
              param8[id]=1.0;
              rejectND(pic[id], id, iran, gaussian_axis);
              
              for (int ip= 0; ip < pic[id].np; ++ip) {
                if (ran3(&iran)*
                    gaussian_axis(id,
                                  INDICIES(pic[id].x[ip],pic[id].y[ip],pic[id].z[ip]))
                    <= 1) {
                  pic[id].vx[ip] = 
                    (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
                  if (ndim >= 2) 
                    pic[id].vy[ip] = 
                      (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
                  if (ndim == 3) 
			  pic[id].vz[ip] = 
			    (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
                  
                } else {
                  pic[id].vx[ip] = 
                    (param9[id]*vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
                  if (ndim >= 2) 
                    pic[id].vy[ip] = 
                      (param9[id]*vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
                  if (ndim == 3) 
                    pic[id].vz[ip] = 
                      (param9[id]*vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
                }
              }
              
            } else if (init_dist[id] == graddrift) {
              
              // 		   this definition of system params is done here
              // 		     because the whole code has not been updated in 
              // 		     such a way
              
              eppic_system sys_params;
              sys_params.ndim=ndim;
              sys_params.nx=nx;
              sys_params.ny=ny;
              sys_params.nz=nz;
              sys_params.dx=dx;
              sys_params.dy=dy;
              sys_params.dz=dz;
              sys_params.eps=eps;
              sys_params.boundary_type[0]=boundary_type[0];
              sys_params.boundary_type[1]=boundary_type[1];
              sys_params.boundary_type[2]=boundary_type[2];
              
              
              // Same parameters as hole function
              // param3 -- gradient length scale in grid units
              FTYPE gradNn = 1./param3[id];
              // param4 -- outside radius in grid units
              PTYPE routside_hole=param4[id];
              // param5 -- inside radius in grid units
              PTYPE rinside_hole=param5[id];
              
              // this get's used by the smooth laying out of particles
              // routine
              int npxdirection = 
                static_cast<int>(sqrt(pic[id].np*
                                      mpi_np/
                                      nsubdomains/
                                      static_cast<FTYPE>(nx*ny))
                                 *nx);
              if (npxdirection==0) 
                npxdirection = static_cast<int>(sqrt(pic[id].np));
              
              int xstart=0,xend=nx;
              if ((boundary_type[0]==INJECT)) {
                if (subdomain.id_number==0)
                  xstart=0;
              }
              ArrayNd<FTYPE,1> prob_distx= 
                ArrayNd<FTYPE,1>(xend-xstart+1);
              ArrayNd<FTYPE,1> prob_disty= 
                ArrayNd<FTYPE,1>(sys_params.ny+1);
              
              int xshift = subdomain.id_number*nx+xstart;
              FTYPE total_nx=nx*nsubdomains;
              
              
              
              if (mpi_rank==0)
                cout << ";; Warning: Parameter definitions have changed" 
                     << " as of 1/25/11. If this input was generated"
                     << " before this date, scale each value by"
                     << " 1/(2*sqrt(2))."
                     << endl;
              if (param6[id]==0) {
                // then gradient along x direction
                if (boundary_type[0] == INJECT) {
                  if (routside_hole==0) routside_hole = total_nx;
                  if (routside_hole>total_nx) routside_hole = total_nx;
                  if (rinside_hole<0) rinside_hole = 0;
                  if (rinside_hole>routside_hole) rinside_hole = routside_hole/4;
                  FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
                  PTYPE routside_hole_sqr=Sqr(param4[id]);
                  PTYPE rinside_hole_sqr=Sqr(param5[id]);
                  for (int ix=xshift;ix<xshift+xend+1;ix++) {
                    FTYPE r = (ix);
                    if (r>routside_hole) 
                      prob_distx(ix-xshift-xstart)=1.0;
                    else 
                      if (r<=routside_hole && 
                          r>rinside_hole) {
                        prob_distx(ix-xshift-xstart)=
                          den_min*exp(gradNn*(r-rinside_hole));
                      }
                      else prob_distx(ix-xshift-xstart)=den_min;
                  }
                } else {
                  if (routside_hole==0) routside_hole = total_nx/2;
                  if (routside_hole>ny) routside_hole = total_nx/2;
                  if (rinside_hole==0) rinside_hole = routside_hole/4;
                  if (rinside_hole>routside_hole) rinside_hole = routside_hole/4;
                  FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
                  PTYPE routside_hole_sqr=Sqr(param4[id]);
                  PTYPE rinside_hole_sqr=Sqr(param5[id]);
                  for (int ix=xshift;ix<xshift+xend+1;ix++) {
                    PTYPE rsqr=Sqr(ix-nx*nsubdomains/2.);
                    if (rsqr>routside_hole_sqr) 
                      prob_distx(ix-xshift-xstart)=1.0;
                    else 
                      if (rsqr<=routside_hole_sqr && 
                          rsqr>rinside_hole_sqr) {
                        prob_distx(ix-xshift-xstart)=
                          den_min*exp(gradNn*(sqrt(rsqr)-rinside_hole));
                      }
                      else prob_distx(ix-xshift-xstart)=den_min;
                  }
                }
                for (int iy=0;iy<sys_params.ny+1;iy++) {
                  prob_disty(iy) = 
                    1.0;
                }
              } else {
                // gradient along y
                if (routside_hole==0) routside_hole = ny/2;
                if (routside_hole>ny) routside_hole = ny/2;
                if (rinside_hole==0) rinside_hole = routside_hole/4;
                if (rinside_hole>routside_hole) rinside_hole = routside_hole/4;
                FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
                PTYPE routside_hole_sqr=Sqr(param4[id]);
                PTYPE rinside_hole_sqr=Sqr(param5[id]);
                for (int iy=0;iy<ny+1;iy++) {
                  PTYPE rsqr=Sqr(iy-ny/2.);
                  if (rsqr>routside_hole_sqr) 
                    prob_disty(iy)=1.0;
                  else if (rsqr<=routside_hole_sqr 
                           && rsqr>rinside_hole_sqr) {
                    prob_disty(iy)=
                      den_min*exp(gradNn*(sqrt(rsqr)-rinside_hole));
                  }
                  else prob_disty(iy)=den_min;
                }
                
		
                for (int ix=xshift;ix<xshift+xend+1;ix++) {
                  prob_distx(ix-xshift-xstart) =  1.0;
                }
              }
              // get max and avg values for scaling np and n0avg
              FTYPE global_max=0,global_avg_max=0,right_area,left_area;
              prob_dist_1D_ddecomposed(prob_distx,global_max,
                                       global_avg_max,
                                       right_area,left_area);
              
              // np is domain specific
              pic[id].np = int(pic[id].np/(prob_distx.length())*
                               FTYPE((total_nx-xstart))/FTYPE(total_nx)*
                               prob_distx.sum()/global_avg_max);

              // n0avg has a domain specific part and a common part
              // related to prob_disty
              pic[id].n0avg*=global_avg_max/global_max;
              pic[id].n0avg*=
                prob_disty.sum()/prob_disty.max()/prob_disty.size(0);
              long long int total_np=0,local_np = pic[id].np;
              MPI_Allreduce(&local_np,&total_np,1,
                            MPI_LONG_LONG_INT,MPI_SUM,
                            subdomain.neighbor_comm);

              
              if (mpi_rank==0) {
                printf("; WARNING: N0avg for PIC adjusted!\n");
                printf("\tn0avgd%1d = %g\n",id,pic[id].n0avg);
              }
              
              
              /* Change boudary conditions n0 values to match specific case
                 Everything should be in units of global_max, so use that
                 to scale prob_dist(0), ie -1, and prob_dist(nx) to 
              */
              
              
              if (subdomain.id_number==0) {
                n0lhsd[id] = prob_distx(1)/global_avg_max*pic[id].n0avg;
                if (subdomain.rank == 0) 
                  cout << "; For dist "
                       << id 
                       << " n0lhs is changed."
                       << endl 
                       << " n0lhsd" << id << " = " << n0lhsd[id] << endl;
              }
              
              if (subdomain.id_number==nsubdomains-1) {
                n0rhsd[id] = prob_distx(nx-1)/global_avg_max*pic[id].n0avg;
                if (subdomain.rank == 0) 
                  cout << "; For dist "
                       << id 
                       << " n0rhs is changed"
                       << endl 
                       << " n0rhsd" << id << " = " << n0rhsd[id] << endl;
              }
              
              
              non_uniform_sampling2x1D(pic[id],
                                       prob_distx,
                                       prob_disty,
                                       npxdirection,
                                       pic[id].np/npxdirection,
                                       right_area,
                                       left_area,
                                       total_np,
                                       sys_params,
                                       xstart,
                                       id*3);
              
              
              if (boundary_type[0] == INJECT) {
                /* force local sine wave in y */
		
                int sin_int = int(param7[id]);
                if (sin_int != 0) {
                  if (sin_int%2 ==0) sin_int+=1;
                  //FTYPE sin_exponent = FTYPE(sin_int);
                  FTYPE old_yval=0;
                  for (int ipart=0;ipart<pic[id].np;ipart++) {
                    if ((pic[id].x(ipart)+xshift>rinside_hole) &
                        (pic[id].x(ipart)+xshift<routside_hole)) {
                      old_yval=pic[id].y[ipart];
                      old_yval+=
                        //  			    exp(-Sqr(2*(total_nx/2-pic[id].x[ipart]-xshift)/
                        // 				     (2.0*total_nx))*param8[id])*
                        //			    param4[id]*ny/(param5[id]*PI)*
                        param8[id]/(2*PI*sin_int)*ny*
                        sin(sin_int*2*PI*pic[id].y[ipart]/ny);
                      pic[id].y[ipart]=old_yval;
                      if (pic[id].y[ipart]>=ny)
                        pic[id].y[ipart]-=ny*ceil(pic[id].y[ipart]/ny);
                      if (pic[id].y[ipart]<0)
                        pic[id].y[ipart]+=ny*ceil(-pic[id].y[ipart]/ny);
                    }// if inside radius
                  }// for each particle
                  
                } // if sin_int, ie param7 is > 0
              } // if injection
              
              
              if (subdomain.rank==0) {
                cout << "domain = " << subdomain.id_number
                     << "   np  = " << pic[id].np << endl;
              }
              
              /* Load velocities according to thermal vel. */
              maxwel_vel_dist_loading(pic[id],
                                      vel_dim[id],
                                      vxthd[id]*dt/dx,
                                      vythd[id]*dt/dy,
                                      vzthd[id]*dt/dz,
                                      id+3);
              
              /* Add velocity drifts to particle velocities */
              if (vel_dim[id]==1) {
                for (long long int ipart=0;ipart<pic[id].np;ipart++){
                  pic[id].vx[ipart]+=vx0d[id]*dt/dx;
                }
              }
              else if (vel_dim[id]==2) {
                for (long long int ipart=0;ipart<pic[id].np;ipart++){
                  pic[id].vx[ipart]+=vx0d[id]*dt/dx;
                  pic[id].vy[ipart]+=vy0d[id]*dt/dy;
                }
              }
              else {
                
                for (long long int ipart=0;ipart<pic[id].np;ipart++){
                  pic[id].vx[ipart]+=vx0d[id]*dt/dx;
                  pic[id].vy[ipart]+=vy0d[id]*dt/dy;
                  pic[id].vz[ipart]+=vz0d[id]*dt/dz;
                }
              }
              
            } else if (init_dist[id] == external_den) {
              // Place particles according to the subroutine rejectND
              // This is for nonuniform distribution in space 
              
              /* define array from reading it in */
              
              
              void gridArrayNd(particle_dist &dist, int id, int &iran, 
                               FArrayND den_arraynd);
              
              void read_domains(FArrayND_ranged &array,
                                MPI_Comm across_comm,
                                char *path_name,
                                int nghost_pts[2]);
              
              if (nsubdomains==1) {
// 		    FArrayND tmp_den=FArrayND(INDICIES(nx,ny,nz));
// 		    FILE* fopensafe(char* filename, char* mode);
// 		    char filename[256];
// 		    FILE *fp = fopensafe(filename,"rb");
// 		    tmp_den.bin_input(fp);
// 		    fclose(fp);
// 		    gridArrayNd(pic[id],id,iran,tmp_den);
                terminate(1,
                          "external_den particle loading not available for 1 domain runs!");
              } else {
                int nghost[2]={0,1};
                int start[]={INDICIES(0,0,0)};
                int range[]={INDICIES(nx+1,ny,nz)};
                FArrayND_ranged tmp_den=
                  FArrayND_ranged(start,range);
                void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, 
                                 const int nx_guard[]);
                char filename[256];
                strcpy(filename,den_file[id].c_str());
                read_domains(tmp_den,subdomain.neighbor_comm,
                             filename,nghost);
                pass_guards(tmp_den,nghost);
                gridArrayNd(pic[id],id,iran,tmp_den);
              }
              
              for (int ip=0; ip<pic[id].np; ++ip) {
                pic[id].vx[ip] = 
                  (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
                if (ndim >= 2) 
                  pic[id].vy[ip] = 
                    (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
                if (ndim == 3) 
                  pic[id].vz[ip] = 
                    (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
              }
            } else if (init_dist[id] == injection) {
              // This leaves the distribution empty to be filled in by injection
              pic[id].np=0;
              
            } else if (init_dist[id] == sine) {
              // Place particles according to the subroutine rejectND
              // This is for nonuniform distribution in space 
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              curr_id=id;
              // if (mpi_rank == 0) printf("%s: Using sine initial distribution \n" \ 
              // 			    "param3 = %f\nparam4 = %f\nparam5 = %f\nparam6 = %f\n",
              // 			    __func__,param3[id],param4[id],param5[id],param6[id]);
              param4[id] /= nx*nsubdomains;
              param6[id] /= ny;
              if (mpi_rank == 0) {
                cout << "WARNING:" << endl;
                cout << "init_particles redefined the following parameters." << endl;
                cout << "  param4[" << id << "] = " << param4[id] << endl;
                cout << "  param6[" << id << "] = " << param6[id] << endl;
              }
              rejectND(pic[id], id, iran,sine_init);
              
              // Place positions randomly:
              for (i = 0; i < pic[id].np; ++i) 
                pic[id].vx[i] = (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
              if (ndim >= 2) for (i = 0; i < pic[id].np; ++i) 
                               pic[id].vy[i] = (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
              if (vel_dim[id] == 3) 
                for (i = 0; i < pic[id].np; ++i) 
                  pic[id].vz[i] = (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
            } else if (init_dist[id] == ablating_meteor) {
              // This section allows for a randomly emitted set of particles.  
              // Initially all the np initialized particles are absent
              pic[id].nabsent = pic[id].np;
              pic[id].np_all = 0;
              for (i=0; i < pic[id].np; ++i){
                pic[id].x[i] = fabsent;
                pic[id].absent[i] = pic[id].np-i-1;
              }
              
              /*
              // Set all particles to absent
              pic[id].x = fabsent; //This specifies that all particles are absent initially.
              pic[id].nabsent=pic[id].absent.size();
              for (i = 0; i < pic[id].absent.size(); ++i) {
                pic[id].absent[i] = pic[id].absent.size()-i-1;
              }
              pic[id].np=0; // This says that initially no particles are present.
              pic[id].np_all=0;
              */
            } else if (init_dist[id] == linear_gradient) {
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              rejectND(pic[id],id,iran,linear_gradient_init);
            } else if (init_dist[id] == dust_test) {
              void rejectND(particle_dist &, int, int &, 
                            PTYPE (den_func)(int id,INDICIES(PTYPE x, PTYPE y, PTYPE z)));
              rejectND(pic[id],id,iran,dust_test_init);
            } else {		  
              cerr << "Unkown Particle initialization method, "<<init_dist[id]<< 
                " in distribution"<< id<<endl;
              terminate(-1," Unknown specified particle initialization");
            }
            
            // Attempt to put particles slightly outside the system back in: 
            if(mpi_rank==0) cout << "; --- Resetting external particles " << endl;
            for (i = 0; i < pic[id].np; ++i) {
              // Make sure particle is not absent (added by Glenn 20180418 becuase of ablating meteor)
              if (pic[id].x[i] > fabsent){
                if (pic[id].x[i] < 0. && boundary_type[0] == PERIODIC) {
                  cout << "Warning: rank=" << mpi_rank << " particle " << i 
                       << " outside: x=" << pic[id].x[i] << " being reentered\n";
                  pic[id].x[i] += nx;
                }
                if (pic[id].x[i] >= nx) {
                  cout << "Warning: rank=" << mpi_rank << " particle " << i 
                       << " outside: x=" << pic[id].x[i] << " being reentered\n";
                  pic[id].x[i] += -nx;
                }
                if (ndim >= 2) {
                  if (pic[id].y[i] < 0.) {
                    cout << "Warning: rank="  << mpi_rank << " particle " << i 
                         << " outside: y=" << pic[id].y[i] << " being reentered\n";
                    pic[id].y[i] += ny; 
                  }
                  if (pic[id].y[i] >= ny) {
                    cout<< "Warning: rank="  << mpi_rank << " particle " << i 
                        << " outside: y=" << pic[id].y[i] << " being reentered\n";
                    pic[id].y[i] += -ny;
                  }
                }
                if (ndim == 3) {
                  if (pic[id].z[i] < 0.) pic[id].z[i] += nz;
                  if (pic[id].z[i] >= nz) pic[id].z[i] += -nz;
                }
              }
            }
            
            // Load particle velocities - normalized to dx/dt: 
            
            if (method[id] < 10) {
              // The particles should be distributed according to the 
              //   actual distribution at t=0, f0 
              
              if (init_dist[id] == 0 || init_dist[id] == 3) {
                if(mpi_rank==0) cout << "; --- Setting up particle velocities " << endl;
                // Use a random gassian distribution 
                for (i = 0; i < pic[id].np; ++i) 
                  pic[id].vx[i] = (vxthd[id]*gasdev(iran) + vx0d[id]) * (dt / dx);
                if (ndim >= 2) for (i = 0; i < pic[id].np; ++i) 
                                 pic[id].vy[i] = (vythd[id]*gasdev(iran) + vy0d[id]) * (dt / dy);
                if (vel_dim[id] == 3) 
                  for (i = 0; i < pic[id].np; ++i) 
                    pic[id].vz[i] = (vzthd[id]*gasdev(iran) + vz0d[id]) * (dt / dz);
              }
              /*
                else if (init_dist[id] == 1) {
                // Use a uniform gaussian distribution 
                void gasuniform(PTYPEAVec &, PTYPEAVec &, int, int, int);
		
                if (ndim == 1)
                gasuniform(pic[id].vx, pic[id].vx, pic[id].np*seed_mod+subdomain.id_number,
                prime(id*6+5), prime(id*6+6) );
                else
                gasuniform(pic[id].vx, pic[id].vy, pic[id].np*seed_mod+subdomain.id_number,
                prime(id*6+5), prime(id*6+6) );
		
                pic[id].vx *= vxthd[id] * (dt / dx);
                if (vx0d[id] !=0.0) pic[id].vx += vx0d[id] * (dt / dx);
		
                if (ndim >=2) {
                pic[id].vy *= vythd[id] * (dt / dy);
                if (vy0d[id] !=0.0) pic[id].vy += vy0d[id]  * (dt / dy);
                }
		
                if (vel_dim[id] == 3) {
                gasuniform(pic[id].vz, pic[id].vz, pic[id].np*seed_mod+subdomain.id_number,
                prime(id*6+7), prime(id*6+8) );
                pic[id].vz *= vzthd[id] * (dt / dz);
                if (vz0d[id] !=0.0) pic[id].vz += vz0d[id] * (dt / dz);
			}
                        
			// Shell distribution particle velocities should be placed as a random gaussian:
			void create_shell_particle(PTYPE vradius, PTYPE vth, PTYPE &vx, PTYPE &vy, PTYPE &vz);
			for (i =  np_bulk; i < np_bulk+np_shell ; ++i) {
                        create_shell_particle(v0_shell[id]*dt/dx, vth_shell[id]*dt/dx, 
                        pic[id].vx(i), pic[id].vy(i), pic[id].vz(i) );
			  pic[id].vy(i) *= dx/dy;
			  pic[id].vz(i) *= dx/dz;
                          }
                          
                          
                          }
              */
            }
          } // iread == 0
          else if (init_dist[id] == graddrift) {
            
            // 		   this definition of system params is done here
            // 		     because the whole code has not been updated in 
            // 		     such a way
            
            eppic_system sys_params;
            sys_params.ndim=ndim;
            sys_params.nx=nx;
            sys_params.ny=ny;
            sys_params.nz=nz;
            sys_params.dx=dx;
            sys_params.dy=dy;
            sys_params.dz=dz;
            sys_params.eps=eps;
            sys_params.boundary_type[0]=boundary_type[0];
            sys_params.boundary_type[1]=boundary_type[1];
            sys_params.boundary_type[2]=boundary_type[2];
            
            
            // Same parameters as hole function
            // param3 -- gradient length scale in grid units
            FTYPE gradNn = 1./param3[id];
            // param4 -- outside radius in grid units
            PTYPE routside_hole=param4[id];
            // param5 -- inside radius in grid units
            PTYPE rinside_hole=param5[id];
            
            // this get's used by the smooth laying out of particles
            // routine
            int npxdirection = 
              static_cast<int>(sqrt(pic[id].np*
                                    mpi_np/
                                    nsubdomains/
                                    static_cast<FTYPE>(nx*ny))
                               *nx);
            if (npxdirection==0) 
              npxdirection = static_cast<int>(sqrt(pic[id].np));
            
            int xstart=0,xend=nx;
            if ((boundary_type[0]==INJECT)) {
              if (subdomain.id_number==0)
                xstart=0;
            }
            ArrayNd<FTYPE,1> prob_distx= 
              ArrayNd<FTYPE,1>(xend-xstart+1);
            ArrayNd<FTYPE,1> prob_disty= 
              ArrayNd<FTYPE,1>(sys_params.ny+1);
            
            int xshift = subdomain.id_number*nx+xstart;
            FTYPE total_nx=nx*nsubdomains;
            
            
            
            if (mpi_rank==0)
              cout << ";; Warning: Parameter definitions have changed" 
                   << " as of 1/25/11. If this input was generated"
                   << " before this date, scale each value by"
                   << " 1/(2*sqrt(2))."
                   << endl;
            if (param6[id]==0) {
              // then gradient along x direction
              if (boundary_type[0] == INJECT) {
                if (routside_hole==0) routside_hole = total_nx;
                if (routside_hole>total_nx) routside_hole = total_nx;
                if (rinside_hole<0) rinside_hole = 0;
                if (rinside_hole>routside_hole) 
                  rinside_hole = routside_hole/4;
                FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
                PTYPE routside_hole_sqr=Sqr(param4[id]);
                PTYPE rinside_hole_sqr=Sqr(param5[id]);
                for (int ix=xshift;ix<xshift+xend+1;ix++) {
                  FTYPE r = (ix);
                  if (r>routside_hole) 
                    prob_distx(ix-xshift-xstart)=1.0;
                  else 
                    if (r<=routside_hole && r>rinside_hole) {
                      prob_distx(ix-xshift-xstart) =
                        den_min*exp(gradNn*(r-rinside_hole));
                    }
                    else prob_distx(ix-xshift-xstart)=den_min;
                }
              } else {
                if (routside_hole==0) routside_hole = total_nx/2;
                if (routside_hole>ny) routside_hole = total_nx/2;
                if (rinside_hole==0) rinside_hole = routside_hole/4;
                if (rinside_hole>routside_hole) 
                  rinside_hole = routside_hole/4;
                FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
                PTYPE routside_hole_sqr=Sqr(param4[id]);
                PTYPE rinside_hole_sqr=Sqr(param5[id]);
                for (int ix=xshift;ix<xshift+xend+1;ix++) {
                  PTYPE rsqr=Sqr(ix-nx*nsubdomains/2.);
                  if (rsqr>routside_hole_sqr) 
                    prob_distx(ix-xshift-xstart)=1.0;
                  else 
                    if (rsqr<=routside_hole_sqr && 
                        rsqr>rinside_hole_sqr) {
                      prob_distx(ix-xshift-xstart)=
                        den_min*exp(gradNn*(sqrt(rsqr)-rinside_hole));
                    }
                    else prob_distx(ix-xshift-xstart)=den_min;
                }
              }
              for (int iy=0;iy<sys_params.ny+1;iy++) {
                prob_disty(iy) = 1.0;
              }
            } else {
              // gradient along y
              if (routside_hole==0) routside_hole = ny/2;
              if (routside_hole>ny) routside_hole = ny/2;
              if (rinside_hole==0) rinside_hole = routside_hole/4;
              if (rinside_hole>routside_hole) 
                rinside_hole = routside_hole/4;
              FTYPE den_min=1./exp(gradNn*(routside_hole-rinside_hole));
              PTYPE routside_hole_sqr=Sqr(param4[id]);
              PTYPE rinside_hole_sqr=Sqr(param5[id]);
              for (int iy=0;iy<ny+1;iy++) {
                PTYPE rsqr=Sqr(iy-ny/2.);
                if (rsqr>routside_hole_sqr) 
                  prob_disty(iy)=1.0;
                else if (rsqr<=routside_hole_sqr 
                         && rsqr>rinside_hole_sqr) {
                  prob_disty(iy)=
                    den_min*exp(gradNn*(sqrt(rsqr)-rinside_hole));
                }
                else prob_disty(iy)=den_min;
              }
              
              
              for (int ix=xshift;ix<xshift+xend+1;ix++) {
                prob_distx(ix-xshift-xstart) =  1.0;
              }
            }
            // get max and avg values for scaling np and n0avg
            FTYPE global_max=0,global_avg_max=0,right_area,left_area;
            prob_dist_1D_ddecomposed(prob_distx,global_max,
                                     global_avg_max,
                                     right_area,left_area);
            
            // np is domain specific
            pic[id].np = int(pic[id].np/(prob_distx.length())*
                             FTYPE((total_nx-xstart))/FTYPE(total_nx)*
                             prob_distx.sum()/global_avg_max);
            
            // n0avg has a domain specific part and a common part
            // related to prob_disty
            pic[id].n0avg*=global_avg_max/global_max;
            pic[id].n0avg*=
              prob_disty.sum()/prob_disty.max()/prob_disty.size(0);
            long long int total_np=0,local_np = pic[id].np;
            MPI_Allreduce(&local_np,&total_np,1,
                          MPI_LONG_LONG_INT,MPI_SUM,
                          subdomain.neighbor_comm);
            
            
            if (mpi_rank==0) {
              printf("; WARNING: N0avg for PIC adjusted!\n");
              printf("\tn0avgd%1d = %g\n",id,pic[id].n0avg);
            }
            
            
            /* Change boudary conditions n0 values to match specific case
               Everything should be in units of global_max, so use that
               to scale prob_dist(0), ie -1, and prob_dist(nx) to 
            */
            
            
            if (subdomain.id_number==0) {
              n0lhsd[id] = prob_distx(1)/global_avg_max*pic[id].n0avg;
              if (subdomain.rank == 0) 
                cout << "; For dist "
                     << id 
                     << " n0lhs is changed."
                     << endl 
                     << " n0lhsd" << id << " = " << n0lhsd[id] << endl;
            }
            
            if (subdomain.id_number==nsubdomains-1) {
              n0rhsd[id] = prob_distx(nx-1)/global_avg_max*pic[id].n0avg;
              if (subdomain.rank == 0) 
                cout << "; For dist "
                     << id 
                     << " n0rhs is changed"
                     << endl 
                     << " n0rhsd" << id << " = " << n0rhsd[id] << endl;
            }
            
            
            
            if (subdomain.rank==0) {
              cout << "domain = " << subdomain.id_number
                   << "   np  = " << pic[id].np << endl;
            }
            
          } // init_dist == graddrift
          else {}
          
          // If iread !=0 then chargeon must = 0
          if(mpi_rank==0) cout << "; --- Setting up chargeon array " << endl;
          if (iread !=0) chargeon[id]=0;
          // If chargeon == 1, each particle needs its own distinct charge 
          if (chargeon[id] == 1)  {
            misc[id].charge=PTYPEAVec(npad) = 0.;
            if (method[id] != 0) {
              terminate(0,"chargeon only implemented for method[id]=0\n");
            }
          }
          
          
	} // END if (method[id] >= 0) 
	
    } // END for (id=0; id<ndist; ++id) 
    
    
    // Setup the velocity damping 
    if(mpi_rank==0) cout << "; --- Setting up velocity damping " << id << endl;
    if (damp_nx < 0) 
      for (id=0; id < ndist; ++id) if (pdamp_nu[id]>0.) {
          printf("\nWarning: pdamp_nu[%d] set while damp_nx not set ... ignored\n",id);
          pdamp_nu[id]=-1;
        }
    
    if (damp_nx >= 0 && damp_nx < nx) {
      for (id=0; id < ndist; ++id) {
	if (pdamp_nu[id]<0.) {
          pdamp_nu[id]=damp_nu;
        }
	misc[id].vx0 = pic[id].vx;
	if (ndim >= 2) {
	  misc[id].vy0 = pic[id].vy;
        }
	if (vel_dim[id]==3 && ndim < 3) {
	  for (i = 0; i < pic[id].np; ++i) 
	    misc[id].vz0[i] = sqrt(Sqr(pic[id].vy[i])+Sqr(pic[id].vz[i]));
	}
	if (ndim >= 3) {
	  misc[id].vz0 = pic[id].vz;
        }
      }
    }
    
    // Velocity modification 
    if(mpi_rank==0) cout << "; --- Setting up velocity modification " << endl;
    void velocity_mod(particle_dist *&pic);
    if (iread == 0) velocity_mod(pic);
    
    // Velocity distribution gather - setup 
    // If user did not set limits - do it 
    FTYPE vmax,vmin;
    
    for (id=0; id < ndist; ++id) {
      if (method[id] > 0) {
	if(mpi_rank==0) cout << "; --- Setting up velocity gather for dist " << id << endl;
	if (nx == 1) pnvx[id]=1;
	if (pvxmax[id] <= pvxmin[id]) { 
	  vmax = pic[id].vx.max()*dx/dt;
	  vmin = pic[id].vx.min()*dx/dt;
	  pvxmin[id] = vmin - (vmax-vmin);
	  pvxmax[id] = vmax + (vmax-vmin);
	}
        
        
	if (ny == 1) pnvy[id]=1;
	if (vel_dim[id] >= 2) {
	  if (pvymax[id] <= pvymin[id]) { 
	    vmax = pic[id].vy.max()*dy/dt;
	    vmin = pic[id].vy.min()*dy/dt;
	    pvymax[id] = vmax + (vmax-vmin);
	    pvymin[id] = vmin - (vmax-vmin);
	  }
	}
	
	// if (nz == 1) pnvz[id]=1;
	if (vel_dim[id] >= 3) {
	  if (pvzmax[id] <= pvzmin[id]) { 
	    vmax = pic[id].vz.max()*dz/dt;
	    vmin = pic[id].vz.min()*dz/dt;
	    pvzmax[id] = vmax + (vmax-vmin);
	    pvzmin[id] = vmin - (vmax-vmin);
	  }
	}
      } // end if (method > 0)
    }
    if(mpi_rank==0) cout << "; --- Particle initialization complete" << endl;
}
