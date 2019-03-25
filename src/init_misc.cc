/* Initialize the various arrays */
#include "eppic.h"
#include "eppic-mpi.h"
#include <stdio.h> 
#include <sys/stat.h>
#include <string.h>
#include <hdf5.h>

extern hid_t H5_FID;

void init_misc(char *arg1)
{
  /* Initialize misc. data */

  

  /*Local variables: */
  int id;

  /* Establish some run charateristics: */

  if (ndim == 3) for (id=0; id<ndist; ++id) vel_dim[id]=3;

  /* Magnetization: */
  if (species_Bx.max() != 0 && species_Bz.max() != 0) {;
    //    terminate(-1,"only Bx !=0 or Bz!=0 (not both) implemented");
  }

  if (Bx != 0) {
    if (Bz != 0){;} //terminate(-1,"only Bx !=0 or Bz!=0 (not both) implemented");
    for (id=0; id<ndist; ++id) {
      if (method[id] == 0 && species_dim[id] >1) {
	vel_dim[id]=3;
	if (ndim == 2) {
	  vzthd[id]=vythd[id];
	  dz=dy;
	}
      }
      else if (method[id] != 0  && method[id] != -2 && method[id] != -3 &&
	       method[id] != -4) {
	if (mpi_rank == 0) 
	  printf("With method = %d, for dist = %d\n",method[id],id);
	terminate(-1,"Requesting unimplemented method");
      }
    }
  }

  if (Bz != 0) {
    for (id=0; id<ndist; ++id) {
      if (method[id] == 0 && species_dim[id] >1) {
	vel_dim[id]=3;
	if (ndim == 2) {
	  vzthd[id]=vythd[id];
	  dz=dy;
	}
      }
      else if (method[id] > 0)
	terminate(-1,"Magnatized plasmas only implemented for straight PIC");
    }
  }

  /* use a C++ temp string to process outdir */
  string temp_outdir (outdir);
  char name[256],infile[256];
  if (arg1 == NULL) 
    sprintf(infile,"%seppic.i","../data/");
  else
    strcpy(infile,arg1);


    
  void trim_strip_quote(std::string &instr);
  trim_strip_quote(temp_outdir); 
  trim_strip_quote(temp_outdir); /* again for wspace in quote */

  /* Make sure last char is slash */
  int find_pos  = temp_outdir.rfind('/');
  if (find_pos != temp_outdir.length()-1) temp_outdir.push_back('/');

  /* copy back */
  strcpy(outdir,temp_outdir.c_str());

  if (mpi_rank==0) {
    /* loop from first to last slash, making each directory */
    find_pos = temp_outdir.find_first_of('/');
    int last_pos = 0;
    while (find_pos>=0) {
      int mkerr=mkdir((temp_outdir.substr(0,
					  find_pos+last_pos)).c_str()
		      ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      last_pos += find_pos+1;
      find_pos = (temp_outdir.substr(last_pos)).find_first_of('/');
    }
    /* copy input to outdir/eppic.i, if it doesn't already exist */
    sprintf(name,"%seppic.i",outdir);
    struct stat     statbuf;

    // if it does not exist, copy input to outdir/eppic.i
    if (stat(name,&statbuf)==-1) // test if it does not exist
      {
	printf(";Creating '%s' \n",name);
	FILE *fpin,*fpout;
	if ((fpin=fopen(infile,"r"))==NULL) {
	  terminate(-1,"Cannot open input file for reading");
	};
	if ((fpout=fopen(name,"w"))==NULL) {
	  terminate(-1,"Cannot open copy of eppic.i for reading");
	};
	char line[130], *fstatus;
	fstatus=fgets(line,130,fpin);
	while (fstatus != (char *) EOF && fstatus != (char *) NULL) {
	  // if line has outdir, prepend a comment character: ';' 
	  string sline(line);
	  if (sline.find("outdir",0) != string::npos) 
	    fputs(";",fpout);
	  fputs(line,fpout);
	  fstatus=fgets(line,130,fpin);
	}
	fclose(fpin);
	fclose(fpout);
	
      }
    // if it exists, copy it to outdir/eppic.i_backup
    else
      {
	FILE *fpin,*fpout;
	if ((fpout=fopen(infile,"r"))==NULL) {
	  terminate(-1,"Cannot open input");
	};
	if ((fpin=fopen(name,"r"))==NULL) {
	  terminate(-1,"Cannot open eppic.i");
	};
	char ch1, ch2, same;
	same = 1;
	// compare the files
	while(!feof(fpin)) {
	  ch1 = fgetc(fpin);
	  if(ferror(fpin)) {
	    printf("Error reading input\n");
	    exit(1);
	  }
	  ch2 = fgetc(fpout);
	  if(ferror(fpout)) {
	    printf("Error reading eppic.i\n");
	    exit(1);
	  }
	  if(ch1 != ch2) {
	    same = 0;
	    break;
	  }
	}
	fclose(fpin);
	fclose(fpout);
	// test to see if input is not eppic.i, otherwise there's
	// no need to copy
	if (!same)
	  {
	    // first make a back up
	    char backup_name[128];
	    sprintf(backup_name,"%seppic.i_backup",outdir);
	    printf(";Creating/Overwriting '%s' \n",backup_name);
	    if ((fpout=fopen(backup_name,"w"))==NULL) {
	      terminate(-1,"Cannot open/create eppic.i_backup");
	    };
	    if ((fpin=fopen(name,"r"))==NULL) {
	      terminate(-1,"Cannot open copy of eppic.i for reading");
	    };
	    char line[130], *fstatus;
	    fstatus=fgets(line,130,fpin);
	    while (fstatus != (char *) EOF && fstatus != (char *) NULL) {
	      fputs(line,fpout);
	      fstatus=fgets(line,130,fpin);
	    }
	    fclose(fpin);
	    fclose(fpout);
	    
	    // then copy into eppic.i the current input
	    printf(";Overwriting '%s' \n",name);

	    if ((fpin=fopen(infile,"r"))==NULL) {
	      terminate(-1,"Cannot open input file for reading");
	    };
	    if ((fpout=fopen(name,"w"))==NULL) {
	      terminate(-1,"Cannot open copy of eppic.i for reading");
	    };

	    fstatus=fgets(line,130,fpin);
	    while (fstatus != (char *) EOF && fstatus != (char *) NULL) {
	      // if line has outdir, prepend a comment character: ';' 
	      string sline(line);
	      if (sline.find("outdir",0) != string::npos) 
		fputs(";",fpout);
	      fputs(line,fpout);
	      fstatus=fgets(line,130,fpin);
	    }
	    fclose(fpin);
	    fclose(fpout);
	
	    
	  } // if input is not eppic.i
      } // else eppic.i does exist

  }
  /* all processors need to wait for first to finish it's work */
  MPI_Barrier(MPI_COMM_WORLD);

  /* Make subdomain directories */
  if (nsubdomains > 0) 
    if (subdomain.rank==0) {
      //outdir[len-1]='\0';
      char temp[128];
      int len=sprintf(temp,"%sdomain%03d/",outdir,subdomain.id_number);
      if (len <0) terminate(-1,"Error in output: Failure to attach domain number to directory name\n");
      int mkerr=mkdir(temp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
  
  //HDF creation moved to output.cc

  // Open the directory used for the restart dumps:
  if (mpi_rank == 0) {
    if (restart_nonlocal) {
      int mkerr=mkdir(strcat(strcpy(name,outdir),"restart"),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      //      system("if [ ! -d restart ]; then mkdir restart; fi");
    } else {
      int mkerr=mkdir("restart",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    int mkerr=mkdir("parallel",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  }
    

  
  // Collision defaults:
  for (id=0; id<ndist; ++id) {
    if (coll_rate[id] > 0) {
      if (vxthd_neutral[id] < 0) vxthd_neutral[id]=vth_neutral;
      if (vythd_neutral[id] < 0) vythd_neutral[id]=vth_neutral;
      if (vzthd_neutral[id] < 0) vzthd_neutral[id]=vth_neutral;
      if (massd_neutral[id] < 0) massd_neutral[id]=m_neutral;
      // Collision algorithms currently require 3D:
	vel_dim[id]=3;
	if (ndim <= 2) {
	  vzthd[id]=vythd[id];
	  dz=dy;
	}
    }
  }

  // Coulomb collisions also need to be in 3D velocity
  for (id=0; id<ndist; ++id) {
    int ptr1, ptr2;
    ptr1 = strcspn(coulomb_type[id],"1"); // For electron-ion routine
    ptr2 = strcspn(coulomb_type[id],"2"); // For 2nd coulomb method (same species)
    if (ptr1 <= ndist || ptr2 <= ndist){
      vel_dim[id]=3;
    }
  }

  /* boundary conditions */
  for (int idim=0; idim<NDIM; idim++){
    if (boundary_type[idim] == INJECT) {
      for (int idist = 0;idist<ndist;idist++) {
        if (n0lhsd[idist]<0){
          printf("WARNING! Particle boundary type is INJECT but n0lhsd was not set for dist %d!\n", idist);
          printf("         We are setting it to %f\n",n0d[idist]);
          n0lhsd[idist]=n0d[idist];
        }
        if (n0rhsd[idist]<0){
          printf("WARNING! Particle boundary type is INJECT but n0rhsd was not set for dist %d!\n", idist);
          printf("         We are setting it to %f\n",n0d[idist]);
          n0rhsd[idist]=n0d[idist];
        }
      }
      //  } else { // changed by Glenn 3/1/2018 to deal with OPEN boundary_type
    } else if (boundary_type[idim] == PERIODIC){
      // added by Glenn 4/16/2018 to check for periodic/field bounary mismatch
      if (field_boundary_type[idim][0] != PERIODIC_FIELD){
        printf("WARNING! Particle boundary type is periodic but the field boundary is not!\n");
      }
      //rhs_boundary = periodic;
      //lhs_boundary = periodic;
      //field_boundary_type[0][0] = PERIODIC_FIELD;
      //field_boundary_type[0][1] = PERIODIC_FIELD;
    }
  }

  // Stuff for quasineutral PETSc solver 
  if (NDIM==2) stencil_size = 9;    // 2-D box stencil perp. to B
  if (NDIM==3) stencil_size = 11;   // + isotropic par. to B 
  global_length = nx*nsubdomains*ny*nz;
  if (nullspace_method == 1) global_length -= 1;
  
  flush(cout);

  return;

} /* End init */

