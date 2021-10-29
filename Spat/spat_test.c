/*
spat_test.c

Test code file for "spat" (spatio-temporal modelling with GDAGsim)

(C) 2002, Darren J Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

*/

#include "spat.h"

void memtest(void);
void varsurface(void);
void meansurface(void);

int main(int argc,char *argv[])
{
  
  /* varsurface(); */
  meansurface();
  /* memtest(); */

  return(EXIT_SUCCESS);
}

void meansurface(void)
{
  int i,j;
  Spat *spat;
  spat=spat_calloc(20,20,20);
  spat_set_latent(spat,10.0,0.8,1.0,1.0);
  spat_add_obs(spat,5,16,18,15.0);
  spat_add_obs(spat,15,15,2,5.0);
  spat_add_obs(spat,19,5,7,12.0);
  spat_process(spat);
  for (i=0;i<(spat->n);i++) {
    for (j=0;j<(spat->m);j++) {
      printf("%f\n",spat_mean(spat,(spat->T)-1,i,j));
      /* printf("%f\n",spat_var(spat,(spat->T)-1,i,j)); */
    }
  }
}

void varsurface(void)
{
  int i,j;
  Spat *spat;
  spat=spat_calloc(30,20,15);
  spat_set_latent(spat,10.0,0.8,1.0,1.0);
  spat_process(spat);
  for (i=0;i<(spat->n);i++) {
    for (j=0;j<(spat->m);j++) {
      /* printf("%f ",spat_mean(spat,(spat->T)-1,i,j)); */
      printf("%f\n",spat_var(spat,(spat->T)-1,i,j));
    }
  }
}

void memtest(void)
{
  int i,j;
  Spat *spat;
  for (i=0;i<10;i++) {
    spat=spat_calloc(20,10,10);
    for (j=0;j<5;j++) {
      spat_set_latent(spat,10.0,0.8,1.0,1.0);
      spat_prior_process(spat);
      spat_add_obs(spat,5,2,3,15.0);
      spat_add_obs(spat,10,3,2,5.0);
      spat_add_obs(spat,15,5,7,12.0);
      spat_process(spat);
      printf("%f %f %f\n",
	     spat_mloglik(spat),
	     spat_mean(spat,15,5,7),
	     spat_var(spat,15,5,7)
	     );
      spat_set_zero(spat);
    }
    spat_free(spat);
  }
}




/* eof */

