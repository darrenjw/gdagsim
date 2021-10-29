/*
sptkriging.c

Code file for "spat" (spatio-temporal modelling with GDAGsim)

(C) 2002, Darren J Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

*/

#include "spat.h"


void add_obs(Spat * spat,char * filename);
void printlastmean(Spat *spat);


int main(int argc,char *argv[])
{
  Spat *spat;
  spat=spat_calloc(31,16,16);
  spat_set_latent(spat,14.0,0.8,0.5,0.5);
  add_obs(spat,"largedata.tab");
  spat_process(spat);
  printlastmean(spat);
  return(EXIT_SUCCESS);
}

void add_obs(Spat * spat,char * filename)
{
  long i;
  int lon,lat,year;
  double temp;
  FILE *s;
  s=fopen(filename,"r");
  if (s == NULL) {
    perror("failed to open file");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<27699;i++) {
    fscanf(s,"%d",&lat);
    fscanf(s,"%d",&lon);
    fscanf(s,"%d",&year);
    fscanf(s,"%lf",&temp);
    spat_add_obs(spat,year-1955,(lat-20)/2,(lon+30)/2,temp);
  }
  fclose(s);
}

void printlastmean(Spat *spat)
{
  int i,j;
  for (i=0;i<(spat->n);i++) {
    for (j=0;j<(spat->m);j++) {
      printf("%f\n",spat_mean(spat,(spat->T)-1,i,j));
      /* printf("%f\n",spat_var(spat,(spat->T)-1,i,j)); */
    }
  }
}


/* eof */

