/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* overriding mpicc to link C programs (only for IBM) */
/* #undef CC */

/* arguments passed to configure script */
#define CONFIGURE_ARGS " '--enable-fftw' '--enable-intel' '--with-fftw=/opt/apps/intel13/mvapich2_1_9/fftw3/3.3.2' 'FC=mpif90' 'CC=mpicc' '--prefix=/home1/02418/moppenhe' '--enable-stride1'"

/* Define if you want to compile P3DFFT using CRAY compiler */
/* #undef CRAY */

/* Define if you want to enable C convention for processor dimensions */
/* #undef DIMS_C */

/* Define if you want to use the ESSL library instead of FFTW */
/* #undef ESSL */

/* Define whether you want to enable estimation */
/* #undef ESTIMATE */

/* Define if you want to use the FFTW library */
#define FFTW 1

/* Define if you want to compile P3DFFT using GNU compiler */
/* #undef GNU */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if you want to compile P3DFFT using IBM compiler */
/* #undef IBM */

/* Define if you want to compile P3DFFT using Intel compiler */
#define INTEL 1

/* Define if you want to enable the measure algorithm */
#define MEASURE 1

/* Define if you want to override the default value of NBL_X */
/* #undef NBL_X */

/* Define if you want to override the default value of NBL_Y1 */
/* #undef NBL_Y1 */

/* Define if you want to override the default value of NBL_Y2 */
/* #undef NBL_Y2 */

/* Define if you want to override the default value of NBL_Z */
/* #undef NBL_Z */

/* Define if you want 1D decomposition */
/* #undef ONED */

/* Name of package */
#define PACKAGE "p3dfft"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "dmitry@sdsc.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "P3DFFT"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "P3DFFT 2.6"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "p3dfft"

/* Define to the home page for this package. */
/* #undef PACKAGE_URL */

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.6"

/* Define if you want to enable the patient algorithm */
/* #undef PATIENT */

/* Define if you want to compile P3DFFT using PGI compiler */
/* #undef PGI */

/* Define if you want to compile P3DFFT in single precision */
/* #undef SINGLE_PREC */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if you want to enable stride-1 data structures */
#define STRIDE1 1

/* Define if you want to MPI_Alltoall instead of MPI_Alltotallv */
/* #undef USE_EVEN */

/* Version number of package */
#define VERSION "2.6"

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */
