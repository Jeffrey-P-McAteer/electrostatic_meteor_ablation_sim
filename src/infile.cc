/* Reads the input data */
/* The format is not truely free - but close*/
/* A better parser is needed but not important */

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-times.h"

void infile(char *infile,int printstat)
{


  FILE *fp;
  char name[130],valstr[130], line[130], *fstatus, *afterequal, tmpCharVal[256];
  int dist=-1;
  int fdistn=-1;


    
  /* See below for the following routines */
  int readvalue(char *name, char *identifier, char *valstr, bool *value,
		int printstat,int id=-1); 
  int readvalue(char *name, char *identifier, char *valstr, int *value,
                int printstat,int id=-1); 
  int readvalue(char *name, char *identifier, char *valstr, FTYPE *value,
		int printstat,int id=-1);
  int readvalue(char *name, char *identifier, char *valstr, char *value,
		int printstat,int id=-1);

  int tmp_bc_type=1;
  
  if (infile == NULL) {
    infile="../data/eppic.i";
    strncpy(outdir,"../data/\0",8);
  } else {
    int len=0;
    char *rest;
    rest=strrchr(infile,'/');
    if (rest != 0 ) {
      len=(int) (rest-infile+1);
      strncpy(outdir,infile,len);
      outdir[len]='\0';
    } else {
      strncpy(outdir,"./\0",3);
    }
    /* initialize outdir */
    //    strncpy(outdir,"./\0",3); 
  }
  if ((fp=fopen(infile,"r"))==NULL) {
    terminate(-1,"Cannot open input file for reading");
  };

  if ((mpi_rank == 0) && printstat) {
      printf("; Version: " PACKAGE_VERSION "\n");
      printf("; SVN Revision: " PACKAGE_REVISION "\n");
      printf("; Please report bugs or send comments to:\n;" 
	     PACKAGE_BUGREPORT "\n\n");
      printf("; Echoing EPPIC input file:\n\n");
  }
  
  /* Read line*/
  fstatus=fgets(line,130,fp);


  while (fstatus != (char *) EOF && fstatus != (char *) NULL)
    {
      /* Parse line into name (before "=") and valstr (after "=")*/
      if (strncmp(fstatus,"\n",1) != 0) sscanf(line,"%s",name);
      else strcpy(name,";\n");
      afterequal=strpbrk(line,"=");
      if (afterequal == NULL) strcpy(valstr,"\n");
      else {
	++afterequal;
	strncpy(valstr,afterequal,130);
      }

      /* search name for recognized components */

      /* First: Is it a comment? */
      if (name[0] == ';' || 
	  name[0] == '#' || 
	  name[0] == '%' ) {  
	if ((mpi_rank == 0)&&printstat) printf("%s",line); 
      }

      // GLOBAL PARAMETERS

      else if (readvalue(name,"title",valstr,title,printstat)== 0);
      else if (readvalue(name,"outdir",valstr,outdir,printstat)== 0) {
	restart_nonlocal=true;
      }
      else if (readvalue(name,"unscale_density",valstr,&unscale_density,printstat)==0);
      else if (readvalue(name,"inject_prop_dc",valstr,&inject_prop_dc,printstat)==0);
      else if (readvalue(name,"boundary_type_x",valstr,&boundary_type[0],printstat)== 0);
      else if (readvalue(name,"boundary_type_y",valstr,&boundary_type[1],printstat)== 0);
      else if (readvalue(name,"boundary_type_z",valstr,&boundary_type[2],printstat)== 0);
      else if (readvalue(name,"boundary_type",valstr,&boundary_type[0],printstat)== 0);
      // field_boundary_types and values added by Glenn for multigrid solver
      else if (readvalue(name,"field_boundary_type_xl",valstr,&field_boundary_type[0][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_type_xh",valstr,&field_boundary_type[0][1],printstat)== 0);
      else if (readvalue(name,"field_boundary_type_yl",valstr,&field_boundary_type[1][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_type_yh",valstr,&field_boundary_type[1][1],printstat)== 0);
      else if (readvalue(name,"field_boundary_type_zl",valstr,&field_boundary_type[2][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_type_zh",valstr,&field_boundary_type[2][1],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_xl",valstr,&field_boundary_values[0][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_xh",valstr,&field_boundary_values[0][1],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_yl",valstr,&field_boundary_values[1][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_yh",valstr,&field_boundary_values[1][1],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_zl",valstr,&field_boundary_values[2][0],printstat)== 0);
      else if (readvalue(name,"field_boundary_value_zh",valstr,&field_boundary_values[2][1],printstat)== 0);
      else if (readvalue(name,"ndim_space",valstr,&ndim_space,printstat)== 0);
      else if (readvalue(name,"nx",valstr,&nx,printstat)== 0);
      else if (readvalue(name,"ny",valstr,&ny,printstat)== 0);
      else if (readvalue(name,"nz",valstr,&nz,printstat)== 0);
      else if (readvalue(name,"dt",valstr,&dt,printstat)== 0);
      else if (readvalue(name,"dx",valstr,&dx,printstat)== 0);
      else if (readvalue(name,"dy",valstr,&dy,printstat)== 0);
      else if (readvalue(name,"dz",valstr,&dz,printstat)== 0);
      else if (readvalue(name,"nt",valstr,&nt,printstat)== 0);
      else if (readvalue(name,"Ermag",valstr,&Ermag,printstat)== 0);
      else if (readvalue(name,"Ex0_external",valstr,&Ex0_external,printstat)== 0);
      else if (readvalue(name,"Ex0_rate",valstr,&Ex0_rate,printstat)== 0);
      else if (readvalue(name,"Ey0_external",valstr,&Ey0_external,printstat)== 0);
      else if (readvalue(name,"Ey0_rate",valstr,&Ey0_rate,printstat)== 0);
      else if (readvalue(name,"Ez0_external",valstr,&Ez0_external,printstat)== 0);
      else if (readvalue(name,"Ez0_rate",valstr,&Ez0_rate,printstat)== 0);
      // Changed by Glenn 20180413 field_boundary_type[0][0,1] deals with this.
      //else if (readvalue(name,"rhs_boundary",valstr,&tmp_bc_type,printstat)==0) 
      //{  rhs_boundary = bc_type(tmp_bc_type);}
      //else if (readvalue(name,"lhs_boundary",valstr,&tmp_bc_type,printstat)==0) 
      //{ lhs_boundary = bc_type(tmp_bc_type);}
      else if (readvalue(name,"lhs_boundary",valstr,&field_boundary_type[0][0],printstat)==0)
        {printf("Warning! lhs_boundary input variable is depreciated.  Use field_boundary_type_xl.\n");}
      else if (readvalue(name,"rhs_boundary",valstr,&field_boundary_type[0][1],printstat)==0)
        {printf("Warning! rhs_boundary input variable is depreciated.  Use field_boundary_type_xh.\n");}
      else if (readvalue(name,"nEext",valstr,&nEext,printstat)== 0);
      else if (readvalue(name,"eps",valstr,&eps,printstat)== 0);
      else if (readvalue(name,"Bx",valstr,&Bx,printstat)== 0);
      else if (readvalue(name,"By",valstr,&By,printstat)== 0);
      else if (readvalue(name,"Bz",valstr,&Bz,printstat)== 0);
      else if (readvalue(name,"Bmag",valstr,&Bmag,printstat)== 0);
      else if (readvalue(name,"Binc",valstr,&Binc,printstat)== 0);
      else if (readvalue(name,"Bdec",valstr,&Bdec,printstat)== 0);
      else if (readvalue(name,"damp_nu",valstr,&damp_nu,printstat)== 0);
      else if (readvalue(name,"damp_nx",valstr,&damp_nx,printstat)== 0);
      else if (readvalue(name,"w0",valstr,&w0,printstat)== 0);
      else if (readvalue(name,"fwidth",valstr,&fwidth,printstat)== 0);
      else if (readvalue(name,"fsteep",valstr,&fsteep,printstat)== 0);
      else if (readvalue(name,"v0x_neutral",valstr,&v0_neutral[0],printstat)== 0);
      else if (readvalue(name,"v0y_neutral",valstr,&v0_neutral[1],printstat)== 0);
      else if (readvalue(name,"v0z_neutral",valstr,&v0_neutral[2],printstat)== 0);
      else if (readvalue(name,"vth_neutral",valstr,&vth_neutral,printstat)== 0);
      else if (readvalue(name,"m_neutral",valstr,&m_neutral,printstat)== 0);
      else if (readvalue(name,"cc_rand",valstr,&cc_rand,printstat)== 0); 
      else if (readvalue(name,"local",valstr,&local,printstat)== 0);
      else if (readvalue(name,"no_parallel_efields",valstr,
			 &no_parallel_efields,printstat)== 0);
      else if (readvalue(name,"kill_modes_ang",valstr,
			 &kill_modes_ang,printstat)== 0);
      else if (readvalue(name,"iontype",valstr,&iontype,printstat)== 0);
      else if (readvalue(name,"nsubdomains",valstr,&nsubdomains,printstat)== 0);
#if HAVE_PETSC
      else if (readvalue(name,"petsc_np",valstr,&petsc_np,printstat)== 0);
#endif
      else if (readvalue(name,"efield_algorithm",valstr,&efield_algorithm,printstat)== 0);
      else if (readvalue(name,"efield_zero_dc",valstr,&efield_zero_dc,printstat)== 0);
      else if (readvalue(name,"electron_dist",valstr,&electron_dist,printstat)== 0);
      else if (readvalue(name,"dust_shape",valstr,&dust_shape,printstat)== 0);
      else if (readvalue(name,"dust_charge",valstr,&dust_charge,printstat)== 0);
      else if (readvalue(name,"dust_den",valstr,&dust_den,printstat)== 0);
      else if (readvalue(name,"dust_sigma_x",valstr,&dust_sigma_x,printstat)== 0);
      else if (readvalue(name,"dust_sigma_y",valstr,&dust_sigma_y,printstat)== 0);
      else if (readvalue(name,"dust_sigma",valstr,&dust_sigma,printstat)== 0);
      else if (readvalue(name,"dust_mid",valstr,&dust_mid,printstat)== 0);
      else if (readvalue(name,"nullspace_method",valstr,&nullspace_method,printstat)== 0);
      else if (readvalue(name,"check_solution",valstr,&check_solution,printstat)== 0);
      else if (readvalue(name,"max_tridiag_solve",valstr,&MAX_TRIDIAG_SOLVE,printstat)== 0);
      else if (readvalue(name,"phift_out_kmax",valstr,&phift_out_kmax,printstat)== 0);
      else if (readvalue(name,"phift_out_min",valstr,&phift_out_min,printstat)== 0);
      //Allow user to tell the system to increase the number of domains per region
/* Developing unequal numbers of processor per domain
      else if (readvalue(name,"ndomain_mult0",valstr,ndomain_mult[0],printstat)== 0);
      else if (readvalue(name,"ndomain_index0",valstr,ndomain_index[0],printstat)== 0);
      else if (readvalue(name,"ndomain_mult1",valstr,ndomain_mult[1],printstat)== 0);
      else if (readvalue(name,"ndomain_index1",valstr,ndomain_index[1],printstat)== 0);
      else if (readvalue(name,"ndomain_mult2",valstr,ndomain_mult[2],printstat)== 0);
      else if (readvalue(name,"ndomain_index2",valstr,ndomain_index[2],printstat)== 0);
      else if (readvalue(name,"ndomain_mult3",valstr,ndomain_mult[3],printstat)== 0);
      else if (readvalue(name,"ndomain_index3",valstr,ndomain_index[3],printstat)== 0);
      else if (readvalue(name,"ndomain_mult4",valstr,ndomain_mult[4],printstat)== 0);
      else if (readvalue(name,"ndomain_index4",valstr,ndomain_index[4],printstat)== 0);
      else if (readvalue(name,"ndomain_mult5",valstr,ndomain_mult[5],printstat)== 0);
      else if (readvalue(name,"ndomain_index5",valstr,ndomain_index[5],printstat)== 0);
      else if (readvalue(name,"ndomain_mult6",valstr,ndomain_mult[6],printstat)== 0);
      else if (readvalue(name,"ndomain_index6",valstr,ndomain_index[6],printstat)== 0);
      else if (readvalue(name,"ndomain_mult7",valstr,ndomain_mult[7],printstat)== 0);
      else if (readvalue(name,"ndomain_index7",valstr,ndomain_index[7],printstat)== 0);
      else if (readvalue(name,"ndomain_mult8",valstr,ndomain_mult[8],printstat)== 0);
      else if (readvalue(name,"ndomain_index8",valstr,ndomain_index[8],printstat)== 0);
      else if (readvalue(name,"ndomain_mult9",valstr,ndomain_mult[9],printstat)== 0);
      else if (readvalue(name,"ndomain_index9",valstr,ndomain_index[9],printstat)== 0);
*/

      // END GLOBAL PARAMETERS


      // IO CONTROLS

      else if (readvalue(name,"nout_avg",valstr,&nout_avg,printstat)== 0);
      else if (readvalue(name,"nout",valstr,&nout,printstat)== 0);
      else if (readvalue(name,"full_array_nout",valstr,&full_array_nout,printstat)== 0);
      else if (readvalue(name,"npout",valstr,&npout,printstat)== 0);
      else if (readvalue(name,"iwrite",valstr,&iwrite,printstat)== 0);
      else if (readvalue(name,"time_limit",valstr,&time_limit,printstat)== 0);
      else if (readvalue(name,"iread",valstr,&iread,printstat)== 0);
      else if (readvalue(name,"phi_out_subcycle",valstr,&divj_out_subcycle,printstat)== 0);
      else if (readvalue(name,"divj_out_subcycle",valstr,&divj_out_subcycle,printstat)== 0);
      else if (readvalue(name,"charge_out_subcycle",valstr,&charge_out_subcycle,printstat)== 0);
      else if (readvalue(name,"hybrid_diagnostic_subcycle",valstr,&hybrid_diagnostic_subcycle,printstat)== 0);
      else if (readvalue(name,"hdf_output_arrays",valstr,&hdf_output_arrays,printstat)== 0);
      else if (readvalue(name,"ft_output_arrays",valstr,&ft_output_arrays,printstat)== 0);
  // END IO CONTROLS

  // Distribution specific information

 else if (readvalue(name,"ndist",valstr,&ndist,printstat)== 0){
   /* initialize the distribution parameter arrays */
   // DEFAULT VALUES
   dist=0;
   npd=intAVec(ndist) = 4096;
   // n0peak=FTYPEAVec(ndist) = 1e6;
   n0d=FTYPEAVec(ndist) = 1e6;
   n0lhsd=FTYPEAVec(ndist) = -1;
   n0rhsd=FTYPEAVec(ndist) = -1;
   vx0lhsd=FTYPEAVec(ndist) = 0;
   vx0rhsd=FTYPEAVec(ndist) = 0;
   vy0lhsd=FTYPEAVec(ndist) = 0;
   vy0rhsd=FTYPEAVec(ndist) = 0;
   vz0lhsd=FTYPEAVec(ndist) = 0;
   vz0rhsd=FTYPEAVec(ndist) = 0;
   vxthd=FTYPEAVec(ndist) = 1.;
   vythd=FTYPEAVec(ndist) = 1.;
   vzthd=FTYPEAVec(ndist) = 1.;
   vth_gb=FTYPEAVec(ndist) = 0;
   vx0d=FTYPEAVec(ndist) = 0.;
   vy0d=FTYPEAVec(ndist) = 0.;
   vz0d=FTYPEAVec(ndist) = 0.;
   n0b=FTYPEAVec(ndist) = 0.;
   vxthb=FTYPEAVec(ndist) = 0.;
   vythb=FTYPEAVec(ndist) = 0.;
   vzthb=FTYPEAVec(ndist) = 0.;
   vx0b=FTYPEAVec(ndist) = 0.;
   vy0b=FTYPEAVec(ndist) = 0.;
   vz0b=FTYPEAVec(ndist) = 0.;
   qd=FTYPEAVec(ndist) = QE;
   md=FTYPEAVec(ndist) = ME;
   species_Bx=FTYPEAVec(ndist) = Bx;
   species_By=FTYPEAVec(ndist) = By;
   species_Bz=FTYPEAVec(ndist) = Bz;
   chargeon=intAVec(ndist) = 0;
   init_dist=intAVec(ndist) = 0;
   coll_rate=FTYPEAVec(ndist) = 0.;
   coll_type=intAVec(ndist) = 0;
   coll_create_id=intAVec(ndist) = 0;
   coll_create_a_id=intAVec(ndist) = 0;
   coll_create_b_id=intAVec(ndist) = 0;
   coll_create_vthb=FTYPEAVec(ndist) = 0;
   coll_start_time=FTYPEAVec(ndist) = 0.;
   crosssec=FTYPEAVec(ndist) = 0.;
   crosssec_m_model=intAVec(ndist) = 0;
   beta_model=intAVec(ndist) = 0;
   e_collision_model=intAVec(ndist) = 0;
   vsim_to_kmsec= 1.;
   lsqsim_to_msq = 1.;
   creation_rate=FTYPEAVec(ndist) = 0.;
   create_pair_dist=intAVec(ndist) = -1;
   annihilation_rate=FTYPEAVec(ndist) = 0.;
   vth_shell=FTYPEAVec(ndist) = 0.;
   v0_shell=FTYPEAVec(ndist) = 0.;
   n0_shell=FTYPEAVec(ndist) = 0.;
   background_neutral_dens=FTYPEAVec(ndist) = 0.;
   coll_cross_section=FTYPEAVec(ndist) = 0.;
   //
   boost_rate=FTYPEAVec(ndist) = 0.;
   vth_boost=FTYPEAVec(ndist) = 0.;
   v0_boost=FTYPEAVec(ndist) = 0.;
   //
   create_vth=FTYPEAVec(ndist) = 0.;
   create_v0=FTYPEAVec(ndist) = 0.;
   create_radius=FTYPEAVec(ndist) = 0.;
   create_posx=FTYPEAVec(ndist) = -1;
   create_posy=FTYPEAVec(ndist) = -1;
   create_posz=FTYPEAVec(ndist) = -1;
   //
   part_pad=FTYPEAVec(ndist) = 1.2;
   vrelmax=FTYPEAVec(ndist) = 0.;
   massd_neutral=FTYPEAVec(ndist) = -1;
   vx0d_neutral=FTYPEAVec(ndist) = v0_neutral[0];
   vy0d_neutral=FTYPEAVec(ndist) = v0_neutral[1];
   vz0d_neutral=FTYPEAVec(ndist) = v0_neutral[2];
   vxthd_neutral=FTYPEAVec(ndist) = -1;
   vythd_neutral=FTYPEAVec(ndist) = -1;
   vzthd_neutral=FTYPEAVec(ndist) = -1;
   method=intAVec(ndist) = 0;
   pnvx=intAVec(ndist) = 1;
   pnvy=intAVec(ndist) = 1;
   pnvz=intAVec(ndist) = 1;
   pvxmin=FTYPEAVec(ndist) = 0.;
   pvymin=FTYPEAVec(ndist) = 0.;
   pvzmin=FTYPEAVec(ndist) = 0.;
   pvxmax=FTYPEAVec(ndist) = 0.;
   pvymax=FTYPEAVec(ndist) = 0.;
   pvzmax=FTYPEAVec(ndist) = 0.;
   pdamp_nu=FTYPEAVec(ndist) = -1.0;
   den_file = stringVec(ndist,"");
   param1=FTYPEAVec(ndist) = 0.;
   param2=FTYPEAVec(ndist) = 0.;
   param3=FTYPEAVec(ndist) = 0.;
   param4=FTYPEAVec(ndist) = 0.;
   param5=FTYPEAVec(ndist) = 0.;
   param6=FTYPEAVec(ndist) = 0.;
   param7=FTYPEAVec(ndist) = 0.;
   param8=FTYPEAVec(ndist) = 0.;
   param9=FTYPEAVec(ndist) = 1.;
   param10=FTYPEAVec(ndist) = 0.;
   start_col=FTYPEAVec(ndist) = 0.;
   stop_col=FTYPEAVec(ndist) = 1.;
   thermal_gamma=FTYPEAVec(ndist) = 1.;
   diffc=FTYPEAVec(ndist) = 0;
   subcycle=intAVec(ndist) = 1;
   coulomb_subcycle=intAVec(ndist) = 1;
   cc_den=FTYPEAVec(ndist) = -1.0;
   supercycle=intAVec(ndist) = 1;
   species_dim=intAVec(ndist) = ndim;
   vel_dim=intAVec(ndist) = ndim;
   den_out_subcycle=intAVec(ndist) = 1;
   part_out_subcycle=intAVec(ndist) = 1;
   vdist_out_subcycle=intAVec(ndist) = 1;
   flux_out_subcycle=intAVec(ndist) = 1;
   nvsqr_out_subcycle=intAVec(ndist) = 1;
   phistore_out_subcycle=intAVec(ndist) = -1;
   phift_out_min=1E-3;
   phift_out_kmax=0.0;
   Eft_out_min=1E-3;
   Eft_out_kmax=0.0;
   denft_out_min=FTYPEAVec(ndist) = 1E-3;
   denft_out_kmax=FTYPEAVec(ndist) = 0.0;
   // END DEFAULT VALUES
 }
 else if (readvalue(name,"dist",valstr,&dist,printstat)== 0) {
   if (dist >= ndist) {
     terminate(-1," Error: dist >= ndist ... Aborting");
   } }

// COLLISION RELATED INPUT
 else if(ndist>0&&dist<ndist&&readvalue(name,"massd_neutral",valstr,&massd_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vx0d_neutral",valstr,&vx0d_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vy0d_neutral",valstr,&vy0d_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vz0d_neutral",valstr,&vz0d_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vthd_neutral",valstr,&vxthd_neutral[dist],printstat,dist)== 0){
   vythd_neutral[dist]=vxthd_neutral[dist];
   vzthd_neutral[dist]=vxthd_neutral[dist]; 
 }
 else if(ndist>0&&dist<ndist&&readvalue(name,"vxthd_neutral",valstr,&vxthd_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vythd_neutral",valstr,&vythd_neutral[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vzthd_neutral",valstr,&vzthd_neutral[dist],printstat,dist)== 0);

// END COLLISION RELATED INPUT

// Transform Scatter Collision Related input
 else if(ndist>0&&dist<ndist&&readvalue(name,"background_neutral_dens",valstr,&background_neutral_dens[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_cross_section",valstr,&coll_cross_section[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_create_a_id",valstr,&coll_create_a_id[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_create_b_id",valstr,&coll_create_b_id[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_create_vthb",valstr,&coll_create_vthb[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_create_id",valstr,&coll_create_id[dist],printstat,dist)== 0);

// General input
 else if(ndist>0&&dist<ndist&&readvalue(name,"npd",valstr,&npd[dist],printstat,dist)== 0);
 // else if(ndist>0&&dist<ndist&&readvalue(name,"n0peakd",valstr,&n0peak[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"n0lhsd",valstr,&n0lhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"n0rhsd",valstr,&n0rhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vx0lhsd",valstr,&vx0lhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vx0rhsd",valstr,&vx0rhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vy0lhsd",valstr,&vy0lhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vy0rhsd",valstr,&vy0rhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vz0lhsd",valstr,&vz0lhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vz0rhsd",valstr,&vz0rhsd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"n0d",valstr,&n0d[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"vthd",valstr,&vxthd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vxthd",valstr,&vxthd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vythd",valstr,&vythd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vzthd",valstr,&vzthd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"v0d",valstr,&vx0d[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vx0d",valstr,&vx0d[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vy0d",valstr,&vy0d[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vz0d",valstr,&vz0d[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"n0b",valstr,&n0b[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vthb",valstr,&vxthb[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vxthb",valstr,&vxthb[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vythb",valstr,&vythb[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vzthb",valstr,&vzthb[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"v0b",valstr,&vx0b[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vx0b",valstr,&vx0b[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vy0b",valstr,&vy0b[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vz0b",valstr,&vz0b[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"qd",valstr,&qd[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"md",valstr,&md[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"species_Bx",valstr,&species_Bx[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"species_By",valstr,&species_By[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"species_Bz",valstr,&species_Bz[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"chargeon",valstr,&chargeon[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"init_dist",valstr,&init_dist[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"part_pad",valstr,&part_pad[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_rate",valstr,&coll_rate[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_start_time",valstr,&coll_start_time[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_type",valstr,&coll_type[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"crosssec_m_model",valstr,&crosssec_m_model[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"crosssec",valstr,&crosssec[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"beta_model",valstr,&beta_model[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"e_collision_model",valstr,&e_collision_model[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"vsim_to_kmsec",valstr,&vsim_to_kmsec,printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"lsqsim_to_msq",valstr,&lsqsim_to_msq,printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_type",valstr,&coll_type[dist],printstat,dist)== 0);
      // Coulomb collision parameters
 else if(ndist>0&&dist<ndist&&readvalue(name,"coul_coll", valstr, &coulomb_type[dist][0],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"coll_create_id",valstr,&coll_create_id[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"vrelmax",valstr,&vrelmax[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"creation_rate",valstr,&creation_rate[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_pair_dist",valstr,&create_pair_dist[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"annihilation_rate",valstr,&annihilation_rate[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vth_shell",valstr,&vth_shell[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"v0_shell",valstr,&v0_shell[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"n0_shell",valstr,&n0_shell[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"boost_rate",valstr,&boost_rate[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vth_boost",valstr,&vth_boost[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"v0_boost",valstr,&v0_boost[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"create_vth",valstr,&create_vth[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_v0",valstr,&create_v0[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_radius",valstr,&create_radius[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_posx",valstr,&create_posx[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_posy",valstr,&create_posy[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"create_posz",valstr,&create_posz[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"method",valstr,&method[dist],printstat,dist)== 0);

 else if(ndist>0&&dist<ndist&&readvalue(name,"pnvx",valstr,&pnvx[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pnvy",valstr,&pnvy[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pnvz",valstr,&pnvz[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvmin",valstr,&pvxmin[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvxmin",valstr,&pvxmin[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvymin",valstr,&pvymin[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvzmin",valstr,&pvzmin[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvmax",valstr,&pvxmax[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvxmax",valstr,&pvxmax[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvymax",valstr,&pvymax[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pvzmax",valstr,&pvzmax[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"pdamp_nu",valstr,&pdamp_nu[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param2",valstr,&param2[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param3",valstr,&param3[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param4",valstr,&param4[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param5",valstr,&param5[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param6",valstr,&param6[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param7",valstr,&param7[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param8",valstr,&param8[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param9",valstr,&param9[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"param10",valstr,&param10[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"start_col",valstr,&start_col[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"stop_col",valstr,&stop_col[dist],printstat,dist)== 0);
 else if(ndist>0
	 &&dist<ndist
	 &&readvalue(name,"den_file",valstr,tmpCharVal,printstat,dist)== 0) {
   std::string tmp_str = tmpCharVal;
   void trim_strip_quote(std::string &instr);
   /* trims, removes quotes */
   trim_strip_quote(tmp_str); 
   /* again to trim extra white space just inside quote */
   trim_strip_quote(tmp_str); 
   den_file[dist] = tmp_str;
 }
 else if(ndist>0&&dist<ndist&&readvalue(name,"thermal_gamma",
					valstr,&thermal_gamma[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"diffc",valstr,&diffc[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"subcycle",valstr,&subcycle[dist],printstat,dist)== 0);
      // Subcycle for coulomb collisions
 else if(ndist>0&&dist<ndist&&readvalue(name,"coulomb_subcycle",valstr,&coulomb_subcycle[dist],printstat,dist)== 0); 
 else if(ndist>0&&dist<ndist&&readvalue(name,"cc_den",valstr,&cc_den[dist],printstat,dist)== 0); 
      
 else if(ndist>0&&dist<ndist&&readvalue(name,"supercycle",valstr,&supercycle[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"species_dim",valstr,
					&species_dim[dist],printstat,dist)== 0);
// Distribution specific IO:
 else if(ndist>0&&dist<ndist&&readvalue(name,"den_out_subcycle",valstr,
					&den_out_subcycle[dist],printstat,dist)== 0);
 else if (ndist>0&&dist<ndist&&readvalue(name,"denft_out_kmax",valstr,
					 &denft_out_kmax[dist],printstat)== 0);
 else if (ndist>0&&dist<ndist&&readvalue(name,"denft_out_min",valstr,
					 &denft_out_min[dist],printstat)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"part_out_subcycle",valstr,
					&part_out_subcycle[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"vdist_out_subcycle",valstr,
					&vdist_out_subcycle[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"flux_out_subcycle",valstr,
					&flux_out_subcycle[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"nvsqr_out_subcycle",valstr,
					&nvsqr_out_subcycle[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"denft_out_min",valstr,
					&denft_out_min[dist],printstat,dist)== 0);
 else if(ndist>0&&dist<ndist&&readvalue(name,"denft_out_kmax",valstr,
					&denft_out_kmax[dist],printstat,dist)== 0);

/* input the parameters to define the output of complete distribution 
   functions, f */
 else if (readvalue(name,"fndist",valstr,&fndist,printstat)== 0){
   /* initialize the f distribution output parameter arrays */
   fdistn=0;
   fnmin=ArrayNd<FTYPE,2>(fndist,2*MAX_NDIM);
   fnmax=ArrayNd<FTYPE,2>(fndist,2*MAX_NDIM);
   fn=ArrayNd<int,2>(fndist,2*MAX_NDIM);
   fn=1;
   fnout=intAVec(fndist) = 1;
   for (int i=0; i<fndist; i++) fof[i]=intAVec(MAXDIST) = 0;
   fnof=intAVec(fndist) = 0;
 }
 else if (readvalue(name,"fdistn",valstr,&fdistn,printstat)== 0) {
   if (fdistn >= fndist) {
     terminate(-1," Error: fdistn >= fndist ... Aborting");
   } }
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fxmin",valstr,&fnmin(fdistn,0),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fxmax",valstr,&fnmax(fdistn,0),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnx",valstr,&fn(fdistn,0),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fymin",valstr,&fnmin(fdistn,1),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fymax",valstr,&fnmax(fdistn,1),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fny",valstr,&fn(fdistn,1),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fzmin",valstr,&fnmin(fdistn,2),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fzmax",valstr,&fnmax(fdistn,2),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnz",valstr,&fn(fdistn,2),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvxmin",valstr,&fnmin(fdistn,3),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvxmax",valstr,&fnmax(fdistn,3),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnvx",valstr,&fn(fdistn,3),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvymin",valstr,&fnmin(fdistn,4),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvymax",valstr,&fnmax(fdistn,4),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnvy",valstr,&fn(fdistn,4),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvzmin",valstr,&fnmin(fdistn,5),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fvzmax",valstr,&fnmax(fdistn,5),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnvz",valstr,&fn(fdistn,5),printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&readvalue(name,"fnout",valstr,&fnout[fdistn],printstat,fdistn)== 0);
 else if(fndist>0&&fdistn<fndist&&
	 readvalue(name,"fof",valstr,&fof[fdistn][fnof[fdistn]++],printstat,fdistn)== 0);
  
 else 
   if ((mpi_rank == 0) && printstat)
     printf(";%s; WARNING: PREVIOUS LINE NOT UNDERSTOOD\n",line);

/* Read next line*/
fstatus=fgets(line,130,fp);
}

fclose(fp);
fflush(NULL);

//  if (mpi_np > 1 && mpi_rank == 0) {
//    char command[132];
//    sprintf(command,
//	    "mailx -s \"EPPIC run started- mpi_np= %d\" meerso@bu.edu <%s",
//	    mpi_np,infile);
//	    system(command); 
//  }


} /* End of routine */
