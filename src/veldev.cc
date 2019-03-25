#include <math.h>
#include "eppic-types.h"
#include "eppic.h"
#include <iostream>

FTYPE veldev(int &idum, int id, FTYPE Lim)
{
  /*     Generate a flux velocity distribution  */

  FTYPE ran3(int *);
  //Vectors where spicies information is saved.
  static FTYPEAVec v0,vth,v_pk,flux_pk,Seg,check;
  static int First_entry=TRUE;  
  
 //Initialization of vth. It will be used to check if the parameter of species has already been saved.
  if (First_entry){
    //v0=FTYPEAVec(ndist);
    //vth=FTYPEAVec(ndist);
    v_pk=FTYPEAVec(ndist);
    flux_pk=FTYPEAVec(ndist);
    Seg=FTYPEAVec(ndist);
    check=FTYPEAVec(ndist);
    for (int idist=0;idist<ndist;idist++){
      //vth[idist]=0.;
      //vth[id]=vxthd[id];
      //v0[id]=vx0d[id];

      //Species parameters are saved in the vectors. This is done just once for each species.
      
      //Maximum velocity value and the respective maximum flux is calculated for a flux distribution (vx*exp()). 
      v_pk[idist]=Lim*(vx0d[idist]+Lim*sqrt(Sqr(vx0d[idist])+4*Sqr(vxthd[idist])))/2.0;
      flux_pk[idist]=v_pk[idist]*exp(-Sqr(v_pk[idist]-Lim*vx0d[idist])/2.0/Sqr(vxthd[idist]));
    
      //Rejection method is used. A velocity segment is needed (Area=Seg*flux_pk) for the different cases.
      if (vx0d[idist]*Lim>0){
	Seg[idist]=8.0*vxthd[idist]; //normal case. 
      }else{
	FTYPE ratio=sqrt(Sqr(vx0d[idist]))/vxthd[idist];
	if (ratio<=5.0){
          Seg[idist]=(8.0-ratio)*vxthd[idist];//case when particles with small velocities are reinjected.
	}else{
          Seg[idist]=3.0*vxthd[idist];//case when particles with even smaller velocities are reinjected.
	}   
      }
      check[idist]=fabs(v_pk[idist]);
    }
    First_entry=FALSE;
   
  } 
  FTYPE v;
    //Case of drift = 0.
    if (vx0d[id]==0.){
       v=vxthd[id]*sqrt(-2.*log(ran3(&idum)));
    } else {
       //case of drift != 0. Rejection method.
      FTYPE Lminor=0;
      if (check[id]>4.*vxthd[id]) Lminor=check[id]-4.*vxthd[id];//Setting up the lower limit of the velocity segment.
       FTYPE TR=0;
          do {
	    v=Seg[id]*ran3(&idum)+Lminor;//velocity of the desired distribution.
             TR=v*exp(-Sqr(v-Lim*vx0d[id])/2.0/Sqr(vxthd[id]));//Threshold. Used to reject or accept v.  
          } while (TR < ran3(&idum)*flux_pk[id]);  
    }
    //cout<<id<<" "<<vxthd[id]<<" "<<vx0d[id]<<" "<<vth[id]<<" "<<v0[id]<<" "<<v_pk[id]<<" "<<flux_pk[id]<<" "<<Seg[id]<<" "<<check[id]<<" "<<Lminor<<" "<<Lim*v<<endl;
    return Lim*v;
}


