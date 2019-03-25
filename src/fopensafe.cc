/* Below are two versions of fopensafe:
   1. Opens file according to the mode
   2. Opens file either at beginning or after skipping a number of bytes 
*/
#include <stdio.h> 
#include <stdlib.h> 

void terminate(int, const char*);

FILE* fopensafe(char* filename, char* mode)
{ 

  FILE *fp;
  if ((fp = fopen(filename,mode)) == NULL ) {
    char string[130];
    sprintf(string,"fopen cannot open %s in mode %s\nABORTING\n", filename, 
	    mode);
    terminate(-1,string);
  }

  return fp;
}

FILE* fopensafe(char* filename, char* mode, unsigned long skip)
{ 


  FILE *fp;
  if (mode[0] != 'a') {
    /* Open files for a new simulation */
    if ((fp = fopen(filename,mode)) == NULL ) {
      char string[130];
      sprintf(string,"fopen cannot open %s in mode %s\nABORTING\n", filename, 
	      mode);
      terminate(-1,string);
    }
  } else {
    /* Open binary file for a preexisting simulation and 
       move forward in the file skip bytes */
    char *mode2;
    mode2="r+\0";
    if ((fp = fopen(filename,mode2)) == NULL ) {
      // Problem: Issue a warning and open the file at the beginning
      printf("Warning: file %s not found. Opening at beginning...\n",
	     filename);
      return fopensafe(filename,"w");
    } else {
      if (mode[1] == 'b') { // advance skip bytes in binary file
	if (fseek(fp, skip, SEEK_SET) != 0) {
	  char string[130];
	  sprintf(string,"fseek cannot advance %s %ld bytes\nABORTING\n", 
		  filename, skip);
	  terminate(-1,string);
	}
      } else { // Advance the file skip lines:
	char dummy[1024];
	for (int i=0;i<skip;i++) fgets(dummy,1024,fp);
	fflush(fp);
	fseek(fp,0,1);
      }
    }
  }

  return fp;
}











