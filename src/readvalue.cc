/* Reads the input data */
/* The format is not truely free - but close*/
/* A better parser is needed but not important */

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "eppic.h"
#include "eppic-mpi.h"


int readvalue(char *name, char *identifier, char *valstr, bool *value,
	      int printstat,int id)
{
  int value_int;
  int readstatus, status;

  if ((readstatus=strncmp(name,identifier,strlen(identifier))) == 0) {
    status=sscanf(valstr,"%d",&value_int);
    if (status == 0 || status == EOF) {
      char string[130];
      sprintf(string,"Error: Valid %s not found\n",identifier);
      terminate(-1,string);
    } else {
      if (value_int == 0){
        *value = false;
      }
      else{
        *value = true;
      }
      if ((mpi_rank == 0) && printstat) 
	if (id >= 0 ) {
	  printf("\t%s%1d = %d\n", identifier, id, value_int);
	} else {
	  printf("\t%s = %d\n", identifier, value_int);
	}
    }
  }

  return readstatus;
}


int readvalue(char *name, char *identifier, char *valstr, int *value,
	      int printstat,int id)
{
  int readstatus, status;

  if ((readstatus=strncmp(name,identifier,strlen(identifier))) == 0) {
    status=sscanf(valstr,"%d",value);
    if (status == 0 || status == EOF) {
      char string[130];
      sprintf(string,"Error: Valid %s not found\n",identifier);
      terminate(-1,string);
    } else {
      if ((mpi_rank == 0) && printstat) 
	if (id >= 0 ) {
	  printf("\t%s%1d = %d\n", identifier, id, *value);
	} else {
	  printf("\t%s = %d\n", identifier, *value);
	}
    }
  }
  return readstatus;
}

int readvalue(char *name, char *identifier, char *valstr, FTYPE *value, 
	      int printstat,int id)
{
  int readstatus, status;
  float fvalue;

  if ((readstatus=strncmp(name,identifier,strlen(identifier))) == 0) {
    status=sscanf(valstr,"%e",&fvalue);
    if (status == 0 || status == EOF) {
      char string[130];
      sprintf(string,"Error: Valid %s not found\n",identifier);
      terminate (-1,string);
    } else {
      *value = FTYPE(fvalue);
      if ((mpi_rank == 0) && printstat) 
	if (id >= 0 ) {
	  printf("\t%s%1d = %g\n", identifier, id, *value);
	} else {
	  printf("\t%s = %g\n", identifier, *value);
	}

    }
  }
  return readstatus;
}

int readvalue(char *name, char *identifier, char *valstr, char *value,
	      int printstat, int id)
{
  int readstatus;

  if ((readstatus=strncmp(name,identifier,strlen(identifier))) == 0) {
    strcpy(value,valstr);
    if ((mpi_rank == 0) && printstat) 
      if (id >= 0 ) {
	printf("\t%s%1d = %s", identifier, id, value);
      } else {
	printf("\t%s = %s", identifier, value);
      }

    /* If quotes surround the string, remove them - Later */
      
  }
  return readstatus;
}





