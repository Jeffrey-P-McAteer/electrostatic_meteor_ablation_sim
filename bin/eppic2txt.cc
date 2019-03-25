
#include "unistd.h"
#include<iostream>
#include "eppic.h"
#include "eppic-mpi.h"
#include "global_defs.h"
//============================================================================
// Global variables that input will set, or needs
// Copied from ../src/main.cc
//============================================================================

//============================================================================
// USAGE:
// Prints a nice usage message if asked for or needed
//============================================================================
void usage() {
  
  printf("USAGE: eppic2latex inputfile\n");
  exit(0);
}


//============================================================================
// MAIN
//============================================================================
main (int argc, char* argv []) {

  extern void infile(char *,int);
  //  extern void check(void);
  void print_infile(int);

  // get options 
  int style=0;
  char *optstring = "hs:";
  int opt = getopt( argc, argv, optstring );
  while( opt != -1 ) {
    switch( opt ) {
    case 'h':
      usage();
      break;
    case 's':
      style = atoi(optarg);
      break;
    case '?':
      usage();
      break;
    default:
      usage();
      break;
    }
    opt = getopt( argc, argv, optstring);
  }
  if (argc-optind!=1) usage();
  // parse infile and print
  infile(argv[optind],0); 
  //  check();
  print_infile(style);
  exit(1);
}

//============================================================================
// Some Global print formats
char *outfmt_int;
char *outfmt_ftype;
char *outfmt_char;
char *dist_outfmt_int;
char *dist_outfmt_ftype;
char *dist_outfmt_char;
//============================================================================

//============================================================================
// Style index
enum {
  idl = 0,
  latex = 1
};
//============================================================================

//============================================================================
// print_infile
// does what it says, takes a style that should a mach one in the above enum
// defaults to idl syntax
// prints everything
//============================================================================

void print_infile (int style) {
  void print_general(char *identifier,int *value);
  void print_general(char *identifier,FTYPE *value);
  void print_general(char *identifier,char *value);
  void print_general_dist(int idist,char *identifier,int *value);
  void print_general_dist(int idist,char *identifier,FTYPE *value);
  void print_general_dist(int idist,char *identifier, char *value);

  cout << "Style = " << style << "\n";
  switch (style)
    {
    case idl:outfmt_int="\t%s = %d\n";
      outfmt_ftype="\t%s = %g\n";
      outfmt_char="\t%s = %s\n";
      dist_outfmt_int="\t%s%d = %d\n";
      dist_outfmt_ftype="\t%s%d = %g\n";
      dist_outfmt_char="\t%s%d = %s\n";
      break;
    case latex:outfmt_int="%s & %d\\\\\n";
      outfmt_ftype="%s & %g\\\\\n";
      outfmt_char="%s & %s\\\\\n";
      dist_outfmt_int="%s%d & %d\\\\\n";
      dist_outfmt_ftype="%s%d & %g\\\\\n";
      dist_outfmt_char="%s%d & %s\\\\\n";
      break;
    default: outfmt_int="\t%s = %d\n";
      outfmt_ftype="\t%s = %g\n";
      outfmt_char="\t%s = %s\n";
      dist_outfmt_int="\t%s%d = %d\n";
      dist_outfmt_ftype="\t%s%d = %g\n";
      dist_outfmt_char="\t%s%d = %s\n";
      break;
    }
  //  print_general("something","yes");
  // Global Parameters
  switch (style) 
    {
    case idl: printf("; System Parameters\n");break;
    case latex:  printf("\\begin{tabular}{|c|c|}\\hline\n\\multicolumn{2}");
      printf("{|c|}{Global Paramters}\\\\\\hline\\hline\n");
      break;
    default: printf("; System Parameters\n");break;
    }
  //  print_general("title",title);
  print_general("ndim_space",&ndim_space);
  print_general("nx",&nx);
  print_general("ny",&ny);
  print_general("nz",&nz);
  print_general("dt",&dt);
  print_general("dx",&dx);
  print_general("dy",&dy);
  print_general("dz",&dz);
  print_general("nt",&nt);
  print_general("Ermag",&Ermag);
  print_general("Ey0_external",&Ey0_external);
  print_general("nEext",&nEext);
  print_general("eps",&eps);
  print_general("Bx",&Bx);
  print_general("Bz",&Bz);
  print_general("damp_nu",&damp_nu);
  print_general("damp_nx",&damp_nx);
  print_general("w0",&w0);
  print_general("fwidth",&fwidth);
  print_general("fsteep",&fsteep);
  print_general("v0x_neutral",&v0_neutral[0]);
  print_general("v0y_neutral",&v0_neutral[1]);
  print_general("v0z_neutral",&v0_neutral[2]);
  print_general("vth_neutral",&vth_neutral);
  print_general("m_neutral",&m_neutral);
  print_general("local",&local);
  print_general("no_parallel_efields",
		&no_parallel_efields);
  print_general("kill_modes_ang",
		&kill_modes_ang);
  print_general("iontype",&iontype);
  //dls /////////////////////////////////////////
  // parallel
  print_general("nsubdomains",&nsubdomains);
  //end dls /////////////////////////////////////

  switch (style) 
    {
    case idl: printf(";\n;\n");break;
    case latex:     printf("\\hline\\end{tabular}\n\n");break;
    default: printf(";\n;\n");break;
    }
  // IO Parameters
  switch (style) 
    {
    case idl: printf("; IO Parameters\n");break;
    case latex:  printf("\\begin{tabular}{|c|c|}\\hline\n\\multicolumn{2}");
      printf("{|c|}{Global Paramters}\\\\\\hline\\hline\n");
      break;
    default: printf("; System Parameters\n");break;
    }
  print_general("nout_avg",&nout_avg);
  print_general("nout",&nout);
  print_general("npout",&npout);
  print_general("iwrite",&iwrite);
  print_general("iread",&iread);
  print_general("divj_out_subcycle",&divj_out_subcycle);
  print_general("charge_out_subcycle",&charge_out_subcycle);

  switch (style) 
    {
    case idl: printf(";\n;\n");break;
    case latex:     printf("\\hline\\end{tabular}\n\n");break;
    default: printf(";\n;\n");break;
    }
  for (
       int idist=0;
       idist<ndist;
       idist++) {
    // idist Parameters
    switch (style) 
      {
      case idl: printf("; System Parameters\n");break;
      case latex:  printf("\\begin{tabular}{|c|c|}\\hline\n\\multicolumn{2}");
	printf("{|c|}{Distribution %d Parameters}\\\\\\hline\\hline\n",idist);
	break;
      default: printf("; System Parameters\n");break;
      }
    print_general_dist(idist,"massd_neutral",&massd_neutral[idist]);
    print_general_dist(idist,"vx0d_neutral",&vx0d_neutral[idist]);
    print_general_dist(idist,"vy0d_neutral",&vy0d_neutral[idist]);
    print_general_dist(idist,"vy0d_neutral",&vy0d_neutral[idist]);
    print_general_dist(idist,"vthd_neutral",&vxthd_neutral[idist]);
    print_general_dist(idist,"vxthd_neutral",&vxthd_neutral[idist]);
    print_general_dist(idist,"vythd_neutral",&vythd_neutral[idist]);
    print_general_dist(idist,"vythd_neutral",&vythd_neutral[idist]);
    // END COLLISION RELATED INPUT
    // General input
    print_general_dist(idist,"npd",&npd[idist]);
    print_general_dist(idist,"n0d",&n0peak[idist]);
    print_general_dist(idist,"vthd",&vxthd[idist]);
    print_general_dist(idist,"vxthd",&vxthd[idist]);
    print_general_dist(idist,"vythd",&vythd[idist]);
    print_general_dist(idist,"vzthd",&vzthd[idist]);
    print_general_dist(idist,"v0d",&vx0d[idist]);
    print_general_dist(idist,"vx0d",&vx0d[idist]);
    print_general_dist(idist,"vy0d",&vy0d[idist]);
    print_general_dist(idist,"vz0d",&vz0d[idist]);
    print_general_dist(idist,"n0b",&n0b[idist]);
    print_general_dist(idist,"vthb",&vxthb[idist]);
    print_general_dist(idist,"vxthb",&vxthb[idist]);
    print_general_dist(idist,"vythb",&vythb[idist]);
    print_general_dist(idist,"vzthb",&vzthb[idist]);
    print_general_dist(idist,"v0b",&vx0b[idist]);
    print_general_dist(idist,"vx0b",&vx0b[idist]);
    print_general_dist(idist,"vy0b",&vy0b[idist]);
    print_general_dist(idist,"vz0b",&vz0b[idist]);
    print_general_dist(idist,"qd",&qd[idist]);
    print_general_dist(idist,"md",&md[idist]);
    print_general_dist(idist,"species_Bx",&species_Bx[idist]);
    print_general_dist(idist,"species_Bz",&species_Bz[idist]);
    print_general_dist(idist,"chargeon",&chargeon[idist]);
    print_general_dist(idist,"init_dist",&init_dist[idist]);
    print_general_dist(idist,"coll_rate",&coll_rate[idist]);
    print_general_dist(idist,"method",&method[idist]);

    print_general_dist(idist,"pnvx",&pnvx[idist]);
    print_general_dist(idist,"pnvy",&pnvy[idist]);
    print_general_dist(idist,"pnvz",&pnvz[idist]);
    print_general_dist(idist,"pvmin",&pvxmin[idist]);
    print_general_dist(idist,"pvxmin",&pvxmin[idist]);
    print_general_dist(idist,"pvymin",&pvymin[idist]);
    print_general_dist(idist,"pvzmin",&pvzmin[idist]);
    print_general_dist(idist,"pvmax",&pvxmax[idist]);
    print_general_dist(idist,"pvxmax",&pvxmax[idist]);
    print_general_dist(idist,"pvymax",&pvymax[idist]);
    print_general_dist(idist,"pvzmax",&pvzmax[idist]);
    print_general_dist(idist,"pdamp_nu",&pdamp_nu[idist]);
    print_general_dist(idist,"param1",&param1[idist]);
    print_general_dist(idist,"param2",&param2[idist]);
    print_general_dist(idist,"param3",&param3[idist]);
    print_general_dist(idist,"param4",&param4[idist]);
    print_general_dist(idist,"param5",&param5[idist]);
    print_general_dist(idist,"thermal_gamma",
		       &thermal_gamma[idist]);
    print_general_dist(idist,"diffc",&diffc[idist]);
    print_general_dist(idist,"subcycle",&subcycle[idist]);
    print_general_dist(idist,"supercycle",&supercycle[idist]);
    print_general_dist(idist,"species_dim",
		       &species_dim[idist]);
    // Distribution specific IO:
    print_general_dist(idist,"den_out_subcycle",
		       &den_out_subcycle[idist]);
    print_general_dist(idist,"part_out_subcycle",
		       &part_out_subcycle[idist]);
    print_general_dist(idist,"vdist_out_subcycle",
		       &vdist_out_subcycle[idist]);
    print_general_dist(idist,"flux_out_subcycle",
		       &flux_out_subcycle[idist]);
    print_general_dist(idist,"nvsqr_out_subcycle",
		       &nvsqr_out_subcycle[idist]);


    switch (style) 
      {
      case idl: printf(";\n;\n");break;
      case latex:     printf("\\end{tabular}\n\n");break;
      default: printf(";\n;\n");break;
      }

  }
  for (int ifndist=0;
       ifndist<fndist;
       ifndist++) {
    
    switch (style) 
      {
      case idl: printf("; System Parameters\n");break;
      case latex:  printf("\\begin{tabular}{|c|c|}\\hline\n\\multicolumn{2}");
	printf("{|c|}{F Distribution %d Parameters}\\\\\\hline\\hline\n",ifndist);
	break;
      default: printf("; System Parameters\n");break;
      }

    /* input the parameters to define the output of complete distribution 
       functions, f */
    print_general_dist(ifndist,"fxmin",&fnmin(ifndist,0));
    print_general_dist(ifndist,"fxmax",&fnmax(ifndist,0));
    print_general_dist(ifndist,"fnx",&fn(ifndist,0));
    print_general_dist(ifndist,"fymin",&fnmin(ifndist,1));
    print_general_dist(ifndist,"fymax",&fnmax(ifndist,1));
    print_general_dist(ifndist,"fny",&fn(ifndist,1));
    print_general_dist(ifndist,"fzmin",&fnmin(ifndist,2));
    print_general_dist(ifndist,"fzmax",&fnmax(ifndist,2));
    print_general_dist(ifndist,"fnz",&fn(ifndist,2));
    print_general_dist(ifndist,"fvxmin",&fnmin(ifndist,3));
    print_general_dist(ifndist,"fvxmax",&fnmax(ifndist,3));
    print_general_dist(ifndist,"fnvx",&fn(ifndist,3));
    print_general_dist(ifndist,"fvymin",&fnmin(ifndist,4));
    print_general_dist(ifndist,"fvymax",&fnmax(ifndist,4));
    print_general_dist(ifndist,"fnvy",&fn(ifndist,4));
    print_general_dist(ifndist,"fvzmin",&fnmin(ifndist,5));
    print_general_dist(ifndist,"fvzmax",&fnmax(ifndist,5));
    print_general_dist(ifndist,"fnvz",&fn(ifndist,5));
    print_general_dist(ifndist,"fnout",&fnout[ifndist]);
    print_general_dist(ifndist,"fof",&fof[ifndist][fnof[ifndist]++]);

    switch (style) 
      {
      case idl: printf(";\n;\n");break;
      case latex:     printf("\\end{tabular}\n\n");break;
      default: printf(";\n;\n");break;
      }
  }
}

//============================================================================
// print_general and print_general_dist
// a set of funcitons that just prints things; It's overkill as there are now
// but it allows for more fancy tricks in the future, like search and replace
// for sepecial characters
//============================================================================

void print_general(char *identifier,int *value) {
  printf(outfmt_int, identifier, *value);
}

void print_general(char *identifier,FTYPE *value) {
  printf(outfmt_ftype, identifier, *value);
}

void print_general(char *identifier,char *value) {
  printf(outfmt_char, identifier, *value);
}

void print_general_dist(int idist, char *identifier,int *value) {
  printf(dist_outfmt_int, identifier, idist, *value);
}

void print_general_dist(int idist, char *identifier,FTYPE *value) {
  printf(dist_outfmt_ftype, identifier, idist, *value);
}

void print_general_dist(int idist, char *identifier,char *value) {
  printf(dist_outfmt_ftype, identifier, idist, *value);
}




