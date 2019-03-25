//----------------------------------------------------

// create particles particles in an arbitrarilly defined distribution
// by placing new particles into the absent arrays.
// This distribution mimics what happens when an ablating neutral population is expanding from a hypersonic particle.

// Meers Oppenheim July 2016
// 

#include "eppic.h"
#include "eppic-mpi.h"
#include "math.h"
#include "Random.h"
#include "eppic-velsample.h"
#include <cmath>

void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                         bool add_absent=false);

// Create a distribution of ablating particles, in distribution id1
void create_particles_ablating(int id1, particle_dist *adv)
{
  static int first_entry=TRUE;
  static PTYPEAVec n_inject_Vpos;
  static velsample *vel_in;
  //int nabsent_start = adv[id1].nabsent;

  static FTYPEAVec creation_fraction;
  static PTYPE domain_xmin;
  
  if (first_entry) {
    first_entry=FALSE;
    vel_in = new velsample[ndist];
    creation_fraction=FTYPEAVec(ndist) = 0.0;
    domain_xmin = nx*dx*(subdomain.id_number);
    PTYPE domain_xmax = nx*dx*(1+subdomain.id_number);
        
    for (int id2=0; id2<ndist; id2++) {
      //init_velsample(vel_in[id2],create_vth[id2],create_v0[id2], 6, 1, 32939+mpi_rank);
      init_velsample_ablating(vel_in[id2],create_vth[id2],create_v0[id2], 6, 1, 32939+mpi_rank);

      // convert the create positions so that they are in position units instead of
      // fractional
      create_posx[id2] *= nx*nsubdomains*dx;
      create_posy[id2] *= ny*dy;
      create_posz[id2] *= nz*dz;

      // Make sure this processor has particles being created within it.  
      PTYPE max_dist_traveled = (create_v0[id1]+create_vth[id1]*10.)*dt;
      FTYPE sphere_xmax = create_posx[id1] + create_radius[id1] + max_dist_traveled;
      FTYPE sphere_xmin = create_posx[id1] - create_radius[id1] - max_dist_traveled; 
      if ( ( domain_xmin < sphere_xmin && sphere_xmin < domain_xmax ) ||
           ( domain_xmin < sphere_xmax && sphere_xmax < domain_xmax ) ||
           ( domain_xmin > sphere_xmin && domain_xmax < sphere_xmax ) ){ 
        // See what fraction of the sphere is in this processor
        // Assume the entire sphere is in the domain and then subtract off areas (ignore pi*r^2)
        FTYPE area_in_domain = 4.0;
        FTYPE radius = (sphere_xmax-sphere_xmin)/2.0;
        if (sphere_xmin < domain_xmin){
          // Part of the sphere is outside of the domain
          area_in_domain -= 2.0*(1.0-(radius-(domain_xmin-sphere_xmin))/radius);
        }
        if (sphere_xmax > domain_xmax){
          area_in_domain -= 2.0*(1.0-(radius-(sphere_xmax-domain_xmax))/radius);
        }
        creation_fraction[id2] = area_in_domain/4.0;
      }
    }

  }// End first_entry
  
  if (id1 == -1 || id1 >= ndist) 
    terminate(-1,"create_particles_ablating does not give appropriate id");
  
  
  //initialize the random number generator:
  const unsigned long long int iran = 1038830+mpi_rank;
  static RanFast rnd(iran);
  
  static FTYPE fnx = (FTYPE) nx;
  PTYPE smallest_float_below(PTYPE x);
  static FTYPE fnx_minus_epsilon=smallest_float_below(fnx);

  // Make sure this processor has particles being created within it.  
  if ( creation_fraction[id1] > 0){
    // Number of particles to create this time step is:
    FTYPE n_part_create = creation_rate[id1]*dt*creation_fraction[id1];
    int int_part_create=static_cast<int>(n_part_create);
    // The fractional particle will contribute:
    if ( (n_part_create-int_part_create) > (rnd.dbl()) ) int_part_create += 1;
    
    if (int_part_create > 0) {        
      int ndim_place=ndim+vel_dim[id1];
      PTYPEAVec place(int_part_create*ndim_place);
      int nplace = 0;
      for (int i=0; i<int_part_create; i++) {
        // Particles emerge from a small sphere of radius, create_radius, centered on create_pos
        // They will have a gaussian velocity flux (meaning they need a gaussian*v distribution)
        // they will originate at any time between t-dt and t.
        // It will be randomly created on the sphere of radius create_radius, centered on create_pos
        // NOTE phi is the polar anlge, theta is azimuthal
        
        // Generate the position
        FTYPE theta=rnd.dbl()*2*PI;
        FTYPE cos_phi=rnd.dbl()*2 -1; // Def. cos_phi=cos(phi)
        FTYPE sin_phi=sqrt(1-Sqr(cos_phi));
        //FTYPE phi = rnd.dbl()*PI;
        //FTYPE cos_phi = cos(phi);
        //FTYPE sin_phi = sin(phi);
        PTYPE xp=create_radius[id1]*sin_phi*cos(theta);
        PTYPE yp=create_radius[id1]*sin_phi*sin(theta);
        PTYPE zp=create_radius[id1]*cos_phi;

        /*
        // Particles must be within a subdomain, flip location if not
        if (xp+create_posx[id1] > domain_xmax || xp+create_posx[id1] < domain_xmin){
          xp = -xp;
          yp = -yp;
          zp = -zp;
        }
        */

        // Generate the velocity
        PTYPE vr=abs(sampleVel(vel_in[id1])+create_v0[id1]);// The sampleVel code is terrible but...
        //PTYPE vr = sqrt(-2*create_vth[id1]*create_vth[id1]*log(1-rnd.dbl()));
        //PTYPE vr = sqrt(2)*create_vth[id1]*erfinv(rnd.dbl());
        //PTYPE vr = sqrt(-log(1-rnd.dbl())*2.0/3.0)*create_vth[id1];
        //PTYPE vr = sqrt(-2.0*log(1-rnd.dbl()))*create_vth[id1]; // pdf = (v/(vth*vth)*exp(-v*v/(2*vth*vth))
        PTYPE vtheta=rnd.dbl()*2*PI;  
        PTYPE vcos_phi=rnd.dbl()*2 -1; // Def. cos_phi=cos(phi)
        PTYPE vsin_phi=sqrt(1-Sqr(vcos_phi));
        //PTYPE vtheta = theta;
        //PTYPE vcos_phi = cos_phi;
        //PTYPE vsin_phi = sin_phi;
        
        //PTYPE vsin_phi = sqrt(rnd.dbl()); // assumes pdf of phi = 2*cos(phi)*sin(phi)
        //PTYPE vcos_phi = sqrt(1.0-vsin_phi*vsin_phi);
        PTYPE vxp=vr*vsin_phi*cos(vtheta);
        PTYPE vyp=vr*vsin_phi*sin(vtheta);
        PTYPE vzp=vr*vcos_phi;
        
        /*    // Rotate the velocity vector into physical space
              FTYPE vxpf = vxp;
              FTYPE vypf = cos_phi*vyp+sin_phi*vzp;
              FTYPE vzpf = cos_phi*vzp-sin_phi*vyp;
              vyp = vypf;
              vzp = vzpf;
              // Rotate by azimuth
              vxpf = cos(vtheta)*vxp+sin(vtheta)*vyp;
              vypf = cos(vtheta)*vyp-sin(vtheta)*vxp;
              vxp=vxpf;
              vyp=vypf;
        */
        // Generate vxyz randomly
        //FTYPE vxp = sqrt(2)*create_vth[id1]*erfinv(2*rnd.dbl()-1);
        //FTYPE vyp = sqrt(2)*create_vth[id1]*erfinv(2*rnd.dbl()-1);
        //FTYPE vzp = sqrt(2)*create_vth[id1]*erfinv(2*rnd.dbl()-1);
        
        // Particles should only be traveling outward, away from meteoroid so v dot x >0
        PTYPE vdotx=vxp*xp+vyp*yp+vzp*zp;
        if (vdotx < 0) {
          // Reflect the particle so it heads out of the ablating sphere instead of into it:
          vxp*=-1.;
          vyp*=-1.;
          vzp*=-1.;
        }
        
        // Move the center of the particle distribution to the correct location
        xp += create_posx[id1];
        yp += create_posy[id1];
        zp += create_posz[id1];
        
        // Subtracting domain_xmin translates this from global to local coordinates
        xp -= domain_xmin; 
        
        // The particle will be randomly created at a time between t-dt and t.
        PTYPE dtp=rnd.dbl()*dt;
        xp += vxp*dtp;
        yp += vyp*dtp;
        zp += vzp*dtp;
        
        // Put the x and v values into the appropriate place array
        if (xp >= 0 && xp < nx*dx){
          // fill_particle_holes will put them into the appropriate arrays.
          place(nplace++)=xp/dx; 
          if(ndim >= 2) place(nplace++)=yp/dy;
          if(ndim >= 3) place(nplace++)=zp/dz;
          place(nplace++)=vxp*dt/dx;
          if (vel_dim[id1] >= 2) place(nplace++)=vyp*dt/dy;
          if (vel_dim[id1] >= 3) place(nplace++)=vzp*dt/dz;
        }
        else{
          // We need to flip the velocity and spacial vectors of the particle to keep it inside the subdomain
          place(nplace++) = (create_posx[id1]-vxp*dtp-domain_xmin-create_radius[id1]*sin_phi*cos(theta))/dx;
          if(ndim >= 2) place(nplace++)=(create_posy[id1]-vyp*dtp-create_radius[id1]*sin_phi*sin(theta))/dy;
          if(ndim >= 3) place(nplace++)=(create_posz[id1]-vzp*dtp-create_radius[id1]*cos_phi)/dz;
          place(nplace++)=-vxp*dt/dx;
          if (vel_dim[id1] >= 2) place(nplace++)=-vyp*dt/dy;
          if (vel_dim[id1] >= 3) place(nplace++)=-vzp*dt/dz;
          /*
          place_absent(nplace_absent++)=xp/dx;
          if(ndim >= 2) place_absent(nplace_absent++)=yp/dy;
          if(ndim >= 3) place_absent(nplace_absent++)=zp/dz;
          place_absent(nplace_absent++)=vxp*dt/dx;
          if (vel_dim[id1] >= 2) place_absent(nplace_absent++)=vyp*dt/dy;
          if (vel_dim[id1] >= 3) place_absent(nplace_absent++)=vzp*dt/dz;
          */
        }
      } // End of for int_part_create
  
      // Insert the new particles inside the domain
      fill_particle_holes(place, nplace, adv[id1], vel_dim(id1));  
      // Update np_all to account for the newly created particles
      adv[id1].np_all += int_part_create;
      
      /*
      // update nabsent_start
      nabsent_start = adv[id1].nabsent;
      // Insert the new particles outside of the domain
      fill_particle_holes(place_absent, nplace_absent, adv[id1], vel_dim(id1), true);
      */
      
    } // End if int_part_create > 0
  } // End if creation region is inside processor's subdomain

  /*
  // Pass particles between processors
  void pass_particles(int id, particle_dist  *adv, int nabsent_start, int vel_dim);
  pass_particles(id1, adv, nabsent_start, vel_dim(id1));
  */

} // End of function
