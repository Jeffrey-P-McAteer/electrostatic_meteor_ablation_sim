/* Return the n_th prime number */
int prime(int n)
{
  static int nprime=4, t=5;
  int i;

  if (n<=3) return n;
  
  /* If a continuous sequence is being generated, do not recalculate */
  if (n != nprime+1) {
    nprime=4;
    t=5;
  }

  for (;;) {
    int top=t/2+1;
    for (i=2; i<top; i++) {
      if (t%i == 0) break;
    }
        
    if (i==top) {
      if (nprime==n) return t;
      nprime++;
    }
    
    t++;
        
  }
}


	      
		    
