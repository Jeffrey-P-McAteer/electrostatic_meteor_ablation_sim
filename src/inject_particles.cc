
//----------------------------------------------------

// inject particles and remove particles on the boundaries
// by placing new particles to be injected into the absent arrays.

// Meers Oppenheim June 2008
// 
// Yann Tambouret, Updated May 2011
// Glenn Sugar, Updated April 2018
//----------------------------------------------------

#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-math.h"
#include "eppic-velsample.h"

void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                         bool add_absent=false);
void inject_particles(particle_dist *adv, int id, int dim, 
                      FTYPE nwest_inject_float, FTYPE neast_inject_float,
                      velsample velwest, velsample veleast);

// Inject particles into the particle distribution array
void inject_particles_x(particle_dist *adv, int id, FArrayND &rho){
  FTYPE gasdev(int &idum);  
  static bool first_entry = true;
  static PTYPEAVec nwest_inject_float(ndist);
  static PTYPEAVec neast_inject_float(ndist);
  static velsample *velwest;
  static velsample *veleast;
  static FTYPE delta_rho = 0;

  if (first_entry) {
    first_entry = false;
    if (subdomain.id_number == 0) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++) {
	sprintf(buffer,"%s/vxinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    } 
    if (subdomain.id_number == nsubdomains-1) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++){
	sprintf(buffer,"%s/vxinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    }
    // Deal with velocity sampling
    veleast = new velsample[ndist];
    velwest = new velsample[ndist];

    FTYPE area = ny*nz;
    FTYPE subdomain_volume = nx*ny*nz;
    // Loop through each id distribution to assign static variables
    for (int id2=0; id2<ndist; id2++){
      // Calculate the number of particles to inject each time step
      // Input n0lhs * Area_plane * dt * subcycle * mean_velocity (0 through infinity)
      FTYPE mean_vel_lhs;
      FTYPE mean_vel_rhs;

      // Calculate the area of the plane and the mean velocity
      // They are dependent on the dimension/axes we are injecting from
      if (subdomain.id_number==0){
        mean_vel_lhs = vxthd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vx0lhsd[id2]/vxthd[id2])/2.0)
          +vx0lhsd[id2]/2.0*(1+erf(vx0lhsd[id2]/vxthd[id2]/sqrt(2.0)));
        init_velsample(velwest[id2],vxthd[id2],vx0lhsd[id2],6,1,55+mpi_rank);
      }
      else{
        mean_vel_lhs = 0;
        init_velsample(velwest[id2],0,0,1,1,55+mpi_rank);
      }
      if (subdomain.id_number==nsubdomains-1){
        mean_vel_rhs = vxthd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vx0rhsd[id2]/vxthd[id2])/2.0)
          -vx0rhsd[id2]/2.0*(1-erf(vx0rhsd[id2]/vxthd[id2]/sqrt(2.0)));
        init_velsample(veleast[id2],vxthd[id2],vx0rhsd[id2],6,-1,56+mpi_rank);
      }
      else{
        mean_vel_rhs = 0;
        init_velsample(veleast[id2],0,0,1,-1,56+mpi_rank);
      }

      // Determine the number of particles to inject
      //nwest_inject_float(id2) = n0lhsd[id2]*area*dt*subcycle[id2]*mean_vel_lhs/subdomain_volume;
      //neast_inject_float(id2) = n0rhsd[id2]*area*dt*subcycle[id2]*mean_vel_rhs/subdomain_volume;
      // Use cell units not length units (1 length unit = 1 cell) area=nx*ny, volume=nx*ny*nz
      nwest_inject_float(id2) = n0lhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_lhs*dt/dx;
      neast_inject_float(id2) = n0rhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_rhs*dt/dx;

      // deal with scaled densities
      if (!unscale_density){
        nwest_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
        neast_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
      }
      
      //      if (nwest_inject_float(id2) > 0 || neast_inject_float(id2) > 0){
      //        printf("proc %d, x inject: nwest_inject_float(%d)= %f, neast_inject_float=%f\n", 
      //               mpi_rank, id2, nwest_inject_float(id2), neast_inject_float(id2));
      //      }
      delta_rho += nwest_inject_float[id2]*adv[id2].q + neast_inject_float[id2]*adv[id2].q;
      
    } // End for id2
  } // End first_entry

  // Adjust the number of particles to inject if it is proporional to rho_dc
  if (inject_prop_dc){
    // This will inject a number proportional to the net charge imbalance in the system
    // Calculate the net charge inside the subdomain
    FTYPE rho_dc=delta_rho;
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<ny;iy++) 
        for(int iz=0;iz<nz;iz++) 
          rho_dc+=rho(INDICIES(ix,iy,iz));
    
    // See if the charge is the same sign as rho_dc
    //if ((rho_dc > 0.5 && adv[id].q > 0) || (rho_dc < 0.5 && adv[id].q < 0)){
    if ((rho_dc > 1.5*adv[id].q && adv[id].q > 0) || (rho_dc < 1.5*adv[id].q && adv[id].q < 0)){
      // Determine the fractional split
      FTYPE ninject_east = rho_dc*neast_inject_float[id]/(nwest_inject_float[id]+
                                                            neast_inject_float[id]);
      FTYPE ninject_west = rho_dc - ninject_east;
      inject_particles(adv, id, 0, nwest_inject_float[id]+ninject_west,
                       nwest_inject_float[id]+ninject_west,
                       velwest[id], veleast[id]);
    }
    else{
      // No modified injection neccessary.  Use default
      inject_particles(adv, id, 0, nwest_inject_float[id], neast_inject_float[id],
                       velwest[id], veleast[id]);
    }
  }
  else{
    inject_particles(adv, id, 0, nwest_inject_float[id], neast_inject_float[id],
                     velwest[id], veleast[id]);
  }
}  // END INJECT_PARTICLES_X

void inject_particles_y(particle_dist *adv, int id, FArrayND &rho){
  static int iran = -55-mpi_rank;
  FTYPE gasdev(int &idum);  
  static bool first_entry = true;
  static PTYPEAVec nwest_inject_float(ndist);
  static PTYPEAVec neast_inject_float(ndist);
  static velsample *velwest;
  static velsample *veleast;
  static FTYPE delta_rho = 0;

  if (first_entry) {
    first_entry = false;
    if (subdomain.id_number == 0) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++) {
	sprintf(buffer,"%s/vyinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    } 
    if (subdomain.id_number == nsubdomains-1) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++){
	sprintf(buffer,"%s/vyinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    }
    // Deal with velocity sampling
    veleast = new velsample[ndist];
    velwest = new velsample[ndist];

    //FTYPE area = nx*dx*nz*dz;
    FTYPE area = nx*nz;
    //FTYPE subdomain_volume = nx*dx*ny*dy*nz*dz;    
    FTYPE subdomain_volume = nx*ny*nz;
    // Loop through each id distribution to assign static variables
    for (int id2=0; id2<ndist; id2++){
      // Calculate the number of particles to inject each time step
      // Input n0lhs * Area_plane * dt * subcycle * mean_velocity (0 through infinity)
      FTYPE mean_vel_lhs;
      FTYPE mean_vel_rhs;

      // Calculate the area of the plane and the mean velocity
      // They are dependent on the dimension/axes we are injecting from
      mean_vel_lhs = vythd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vy0lhsd[id2]/vythd[id2])/2.0)
        +vy0lhsd[id2]/2.0*(1+erf(vy0lhsd[id2]/vythd[id2]/sqrt(2.0)));
      mean_vel_rhs = vythd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vy0rhsd[id2]/vythd[id2])/2.0)
        -vy0rhsd[id2]/2.0*(1-erf(vy0rhsd[id2]/vythd[id2]/sqrt(2.0)));
      init_velsample(velwest[id2],vythd[id2],vy0lhsd[id2],6,1,55+mpi_rank);
      init_velsample(veleast[id2],vythd[id2],vy0rhsd[id2],6,-1,56+mpi_rank);

      // compute the number of particles to inject on each side of domain
      //nwest_inject_float(id2) = n0lhsd[id2]*area*dt*subcycle[id2]*mean_vel_lhs/subdomain_volume;
      //neast_inject_float(id2) = n0rhsd[id2]*area*dt*subcycle[id2]*mean_vel_rhs/subdomain_volume;
      // Use cell units not length units (1 length unit = 1 cell) area=nx*ny, volume=nx*ny*nz
      nwest_inject_float(id2) = n0lhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_lhs*dt/dy;
      neast_inject_float(id2) = n0rhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_rhs*dt/dy;
      
      // deal with scaled densities
      if (!unscale_density){
        nwest_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
        neast_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
      }

      if (nwest_inject_float(id2) > 0 || neast_inject_float(id2) > 0){
        printf("proc %d, y inject: nwest_inject_float(%d)= %f, neast_inject_float=%f\n", 
               mpi_rank, id2, nwest_inject_float(id2), neast_inject_float(id2));
      }

      delta_rho += nwest_inject_float[id2]*adv[id2].q + neast_inject_float[id2]*adv[id2].q;
    } // End for id2
  } // End first_entry

  if (inject_prop_dc){
    // This will inject a number proportional to the net charge imbalance in the system
    // Calculate the net charge inside the subdomain
    FTYPE rho_dc=delta_rho;
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<ny;iy++) 
        for(int iz=0;iz<nz;iz++) 
          rho_dc+=rho(INDICIES(ix,iy,iz));
    
    // See if the charge is the same sign as rho_dc (so we should inject more)
//    if ((rho_dc > 0.5 && adv[id].q > 0) || (rho_dc < 0.5 && adv[id].q < 0)){
    if ((rho_dc > 1.5*adv[id].q && adv[id].q > 0) || (rho_dc < 1.5*adv[id].q && adv[id].q < 0)){
      // Determine the fractional split
      FTYPE ninject_east = rho_dc*neast_inject_float[id]/(nwest_inject_float[id]+
                                                            neast_inject_float[id]);
      FTYPE ninject_west = rho_dc - ninject_east;
      inject_particles(adv, id, 1, nwest_inject_float[id]+ninject_west,
                       nwest_inject_float[id]+ninject_west,
                       velwest[id], veleast[id]);
    }
    else{
      // No modified injection neccessary.  Use default
      inject_particles(adv, id, 1, nwest_inject_float[id], neast_inject_float[id],
                       velwest[id], veleast[id]);
    }
  }
  else{
    inject_particles(adv, id, 1, nwest_inject_float[id], neast_inject_float[id],
                     velwest[id], veleast[id]);
  }

}  // END INJECT_PARTICLES_Y

void inject_particles_z(particle_dist *adv, int id, FArrayND &rho){
  static int iran = -55-mpi_rank;
  FTYPE gasdev(int &idum);  
  static bool first_entry = true;
  static PTYPEAVec nwest_inject_float(ndist);
  static PTYPEAVec neast_inject_float(ndist);
  static velsample *velwest;
  static velsample *veleast;
  static FTYPE delta_rho = 0;

  if (first_entry) {
    first_entry = false;
    if (subdomain.id_number == 0) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++) {
	sprintf(buffer,"%s/vzinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    } 
    if (subdomain.id_number == nsubdomains-1) {
      char buffer[256];
      for (int id2=0;id2<ndist;id2++){
	sprintf(buffer,"%s/vzinject_id%d_rank%d.dat\0",outdir,id2,mpi_rank);
      }
    }
    // Deal with velocity sampling
    veleast = new velsample[ndist];
    velwest = new velsample[ndist];
    
    // Loop through each id distribution to assign static variables
    //FTYPE area = nx*dx*ny*dy;
    FTYPE area = nx*ny;
    //FTYPE subdomain_volume = nx*dx*ny*dy*nz*dz;
    FTYPE subdomain_volume = nx*ny*nz;
    for (int id2=0; id2<ndist; id2++){
      // Calculate the number of particles to inject each time step
      // Input n0lhs * Area_plane * dt * subcycle * mean_velocity (0 through infinity)
      FTYPE mean_vel_lhs;
      FTYPE mean_vel_rhs;

      // Calculate the area of the plane and the mean velocity
      // They are dependent on the dimension/axes we are injecting from
      // from PICARD: (vth/(sqrt(2*pi))*exp(-v0*v0/(2*vth*vth)) + (v0/2)*(direction + erf(v0/(vth*sqrt(2))))
      mean_vel_lhs = vzthd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vz0lhsd[id2]/vzthd[id2])/2.0)
        +vz0lhsd[id2]/2.0*(1+erf(vz0lhsd[id2]/vzthd[id2]/sqrt(2.0)));
      mean_vel_rhs = vzthd[id2]/sqrt(2.0*M_PI)*exp(-Sqr(vz0rhsd[id2]/vzthd[id2])/2.0)
        -vz0rhsd[id2]/2.0*(1-erf(vz0rhsd[id2]/vzthd[id2]/sqrt(2.0)));
      init_velsample(velwest[id2],vzthd[id2],vz0lhsd[id2],6,1,55+mpi_rank);
      init_velsample(veleast[id2],vzthd[id2],vz0rhsd[id2],6,-1,56+mpi_rank);

      //nwest_inject_float(id2) = n0lhsd[id2]*area*dt*subcycle[id2]*mean_vel_lhs/subdomain_volume;
      //neast_inject_float(id2) = n0rhsd[id2]*area*dt*subcycle[id2]*mean_vel_rhs/subdomain_volume;
      // Use cell units not length units (1 length unit = 1 cell) area=nx*ny, volume=nx*ny*nz
      nwest_inject_float(id2) = n0lhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_lhs*dt/dz;
      neast_inject_float(id2) = n0rhsd[id2]/subdomain_volume*area*subcycle[id2]*mean_vel_rhs*dt/dz;
      
      // deal with scaled densities
      if (!unscale_density){
        nwest_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
        neast_inject_float(id2) *= npd[id2]/adv[id2].n0avg;
      }

      //      if (nwest_inject_float(id2) > 0 || neast_inject_float(id2) > 0){
      //printf("proc %d, z inject: nwest_inject_float(%d)= %f, neast_inject_float=%f\n", 
      //       mpi_rank, id2, nwest_inject_float(id2), neast_inject_float(id2));
      //}

      delta_rho += nwest_inject_float[id2]*adv[id2].q + neast_inject_float[id2]*adv[id2].q;
    } // End for id2
  } // End first_entry

  // Adjust the number of particles to inject if it is proporional to rho_dc

  if (inject_prop_dc){
    // This will inject a number proportional to the net charge imbalance in the system
    // Calculate the net charge inside the subdomain
    FTYPE rho_dc=delta_rho;
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<ny;iy++) 
        for(int iz=0;iz<nz;iz++) 
          rho_dc+=rho(INDICIES(ix,iy,iz));
    
    // See if the charge is the same sign as rho_dc
    //if ((rho_dc > 0.5 && adv[id].q > 0) || (rho_dc < 0.5 && adv[id].q < 0)){
    if ((rho_dc > 1.5*adv[id].q && adv[id].q > 0) || (rho_dc < 1.5*adv[id].q && adv[id].q < 0)){
      // Determine the fractional split
      FTYPE ninject_east = rho_dc*neast_inject_float[id]/(nwest_inject_float[id]+
                                                            neast_inject_float[id]);
      FTYPE ninject_west = rho_dc - ninject_east;
      inject_particles(adv, id, 2, nwest_inject_float[id]+ninject_west,
                       nwest_inject_float[id]+ninject_west,
                       velwest[id], veleast[id]);
    }
    else{
      // No modified injection neccessary.  Use default
      inject_particles(adv, id, 2, nwest_inject_float[id], neast_inject_float[id],
                       velwest[id], veleast[id]);
    }
  }
  else{
    inject_particles(adv, id, 2, nwest_inject_float[id], neast_inject_float[id],
                     velwest[id], veleast[id]);
  }
}  // END INJECT_PARTICLES_Z


void inject_particles(particle_dist *adv, int id, int dim, 
                      FTYPE nwest_inject_float, FTYPE neast_inject_float,
                      velsample velwest, velsample veleast){
  static int iran = -55-mpi_rank;
  //static PTYPE fnx = (PTYPE) nx;
  PTYPE smallest_float_below(PTYPE x);
  //static PTYPE fnx_minus_epsilon=smallest_float_below(fnx);
  static PTYPE fnx_minus_epsilon=smallest_float_below(static_cast<PTYPE>(nx));
  static PTYPE fny_minus_epsilon=smallest_float_below(static_cast<PTYPE>(ny));
  static PTYPE fnz_minus_epsilon=smallest_float_below(static_cast<PTYPE>(nz));
  static PTYPE epsilon=smallest_float_above(PTYPE(0));

  // See if we are in the 0th subdomain or we are injecting in the y,z dimension:
  if (subdomain.id_number==0 || dim!=0 ) {
    // We should inject particles from the left (west) with positive velocities
    int n_i_west=static_cast<int>(nwest_inject_float);
    // Deal with decimal valued n_inject
    if ( (nwest_inject_float-n_i_west) > ran3(&iran) ) n_i_west += 1;

    if (n_i_west > 0){
      int nwest=0;
      PTYPEAVec recv_west((ndim+vel_dim[id])*n_i_west);
      for (int i=0;i<n_i_west;i++) {
        //FTYPE vx,x;
        FTYPE v_inject,p_inject;
        //vx = sampleVel(velwest[id])*dt/dx;
        
        v_inject = sampleVel(velwest);
	
        if (dim == 0){
          v_inject *= dt/dx;
          p_inject=v_inject*(1-ran3(&iran))*subcycle[id]+epsilon;
          recv_west(nwest++)=p_inject;
          if(ndim > 1) recv_west(nwest++)=ran3(&iran)*fny_minus_epsilon;
          if(ndim > 2) recv_west(nwest++)=ran3(&iran)*fnz_minus_epsilon;
          recv_west(nwest++)=v_inject;
          if (vel_dim[id]>=2) recv_west(nwest++)=(vythd[id]*gasdev(iran)+vy0d[id])*(dt/dy);
          if (vel_dim[id]==3) recv_west(nwest++)=(vzthd[id]*gasdev(iran)+vz0d[id])*(dt/dz);
        }
        else if (dim == 1){
          v_inject *= dt/dy;
          p_inject=v_inject*(1-ran3(&iran))*subcycle[id]+epsilon;
          recv_west(nwest++)=ran3(&iran)*fnx_minus_epsilon;
          if (ndim > 1) recv_west(nwest++)=p_inject;
          if (ndim > 2) recv_west(nwest++)=ran3(&iran)*fnz_minus_epsilon;
          recv_west(nwest++)=(vxthd[id]*gasdev(iran)+vx0d[id])*(dt/dx);
          if (vel_dim[id]>=2) recv_west(nwest++)=v_inject;
          if (vel_dim[id]==3) recv_west(nwest++)=(vzthd[id]*gasdev(iran)+vz0d[id])*(dt/dz);
        }
        else if (dim == 2){
          v_inject *= dt/dz;
          p_inject=v_inject*(1-ran3(&iran))*subcycle[id]+epsilon;
          recv_west(nwest++)=ran3(&iran)*fnx_minus_epsilon;
          if (ndim > 1) recv_west(nwest++)=ran3(&iran)*fny_minus_epsilon;;
          if (ndim > 2) recv_west(nwest++)=p_inject;
          recv_west(nwest++)=(vxthd[id]*gasdev(iran)+vx0d[id])*(dt/dx);
          if (vel_dim[id]>=2) recv_west(nwest++)=(vythd[id]*gasdev(iran)+vy0d[id])*(dt/dy);
          if (vel_dim[id]==3) recv_west(nwest++)=v_inject;
        }
      } // end for n_i_west
      
      // Put the particles into the particle structure
      fill_particle_holes(recv_west, nwest, adv[id], vel_dim[id]);
    } // end if n_i_west>0
  } // end if subdomain.id==0

  // See if we are in the last subdomain:
  if (subdomain.id_number==nsubdomains-1 || dim!=0) {
    // We should inject particles from the right (east) with negative velocities
    int n_i_east=static_cast<int>(neast_inject_float);
    // Deal with decimal valued n_inject
    if ( (neast_inject_float-n_i_east) > ran3(&iran) ) n_i_east += 1;

    if (n_i_east>0){
      int neast=0;
      PTYPEAVec recv_east((ndim+vel_dim[id])*n_i_east);
      for (int i=0;i<n_i_east;i++) {
        //FTYPE vx,x;
        FTYPE v_inject,p_inject;
        v_inject = sampleVel(veleast);
	
        if (dim == 0){
          v_inject *= dt/dx;
          p_inject = fnx_minus_epsilon+v_inject*(1-ran3(&iran))*subcycle[id];
          recv_east(neast++)=p_inject;
          if(ndim > 1) recv_east(neast++)=ran3(&iran)*fny_minus_epsilon;
          if(ndim > 2) recv_east(neast++)=ran3(&iran)*fnz_minus_epsilon;
          recv_east(neast++)=v_inject;
          if (vel_dim[id]>=2) recv_east(neast++)=(vythd[id]*gasdev(iran)+vy0d[id])*(dt/dy);
          if (vel_dim[id]==3) recv_east(neast++)=(vzthd[id]*gasdev(iran)+vz0d[id])*(dt/dz);
        }
        else if (dim == 1){
          v_inject *= dt/dy;
          p_inject = fny_minus_epsilon+v_inject*(1-ran3(&iran))*subcycle[id];
          recv_east(neast++)=ran3(&iran)*fnx_minus_epsilon;
          if (ndim > 1) recv_east(neast++)=p_inject;
          if (ndim > 2) recv_east(neast++)=ran3(&iran)*fnz_minus_epsilon;
          recv_east(neast++) = (vxthd[id]*gasdev(iran)+vx0d[id])*(dt/dx);
          if (vel_dim[id]>=2) recv_east(neast++)=v_inject;
          if (vel_dim[id]==3) recv_east(neast++)=(vzthd[id]*gasdev(iran)+vz0d[id])*(dt/dz);
        }
        else if (dim == 2){
          v_inject *= dt/dz;
          p_inject = fnz_minus_epsilon+v_inject*(1-ran3(&iran))*subcycle[id];
          recv_east(neast++)=ran3(&iran)*fnx_minus_epsilon;
          if (ndim > 1) recv_east(neast++)=ran3(&iran)*fny_minus_epsilon;;
          if (ndim > 2) recv_east(neast++)=p_inject;
          recv_east(neast++)=(vxthd[id]*gasdev(iran)+vx0d[id])*(dt/dx);
          if (vel_dim[id]>=2) recv_east(neast++)=(vythd[id]*gasdev(iran)+vy0d[id])*(dt/dy);
          if (vel_dim[id]==3) recv_east(neast++)=v_inject;
        }
      } // end for n_i_east
      
      // Put the particles into the particle structure
      fill_particle_holes(recv_east, neast, adv[id], vel_dim[id]);
    } // end if n_i_east > 0
  } // end if subdomain.id==nsubdomains-1
} // END INJECT_PARTICLES
