/** returns the num_th pseudo-random number, from 0 to 1, using the 
    prime number n.

    \bug There is an error that crops up when this routine is used:
    " void value not ignored as it ought to be "

*/

#include "eppic-types.h"
FTYPE reverse(int num, int n)
{

    /* Local variables */
    int irem, inum;
    FTYPE power;
    int iquot;
    FTYPE rev;


/*     Function to reverse the digits of num in base n. */

    rev = 0.;
    inum = num;
    power = 1.;

    do {
	iquot = inum / n;
	irem = inum - n * iquot;
	power /= n;
	rev += irem * power;
	inum = iquot;
    } while(inum > 0);

    return rev;
} /* End reverse */


FTYPE reverse(long long int num, int n)
{

    /* Local variables */
    long long int irem, inum;
    FTYPE power;
    long long int iquot;
    FTYPE rev;


/*     Function to reverse the digits of num in base n. */

    rev = 0.;
    inum = num;
    power = 1.;

    do {
	iquot = inum / n;
	irem = inum - n * iquot;
	power /= n;
	rev += irem * power;
	inum = iquot;
    } while(inum > 0);

    return rev;
} /* End reverse */
