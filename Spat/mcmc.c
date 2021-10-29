/*
mcmc.c

Code file for MCMC in spat models

(C) 2002, Darren J Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

*/

#include <math.h>
#include "spat.h"

#define NUMOBS 27699

/* function prototypes */
void add_obs(Spat * spat,gsl_vector_int *t_vec,
	     gsl_vector_int *i_vec,gsl_vector_int *j_vec,
	     gsl_vector *temp_vec);
void read_data(char *filename,gsl_vector_int *t_vec,
	       gsl_vector_int *i_vec,gsl_vector_int *j_vec,
	       gsl_vector *temp_vec);


/* main function */
int main(int argc,char *argv[])
{
  Spat *spat;
  gsl_vector_int *t_vec,*i_vec,*j_vec;
  gsl_vector *temp_vec;
  gsl_rng *rng;
  long it,iters;
  double a,mll,mll_p,lambda_s,lambda_s_p,lambda_o,lambda_o_p,mu,mu_p;
  if (argc != 2) {
    fprintf(stderr,"Usage: %s <iters>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  iters=atoi(argv[1]);
  spat=spat_calloc(31,16,16);
  temp_vec=gsl_vector_calloc(NUMOBS);
  t_vec=gsl_vector_int_calloc(NUMOBS);
  i_vec=gsl_vector_int_calloc(NUMOBS);
  j_vec=gsl_vector_int_calloc(NUMOBS);
  read_data("largedata.dat",t_vec,i_vec,j_vec,temp_vec);
  rng=gsl_rng_alloc(gsl_rng_mt19937);

  mu=14.0; lambda_s=-0.5; lambda_o=-0.5; mll=-1e50;  /* inits */
  
  printf("Iter mu lambda_s lambda_o mll\n");
  for (it=0;it<iters;it++) {
    mu_p = mu + gsl_ran_gaussian(rng,0.1);
    lambda_s_p = lambda_s + gsl_ran_gaussian(rng,0.01);
    lambda_o_p = lambda_o + gsl_ran_gaussian(rng,0.01);
    spat_set_zero(spat);
    spat_set_latent(spat,mu_p,0.8,exp(lambda_s_p),exp(lambda_o_p));
    spat_prior_process(spat);
    add_obs(spat,t_vec,i_vec,j_vec,temp_vec);
    spat_process(spat);
    mll_p = spat_mloglik(spat);
    a = 1.0*(lambda_s_p-lambda_s + lambda_o_p-lambda_o)
      + mll_p - mll;
    if (gdag_accept_lp(rng,a)) {
      mu=mu_p; lambda_s=lambda_s_p; lambda_o=lambda_o_p; mll=mll_p;
    }
    printf("%ld %f %f %f %f\n",it,mu,lambda_s,lambda_o,mll);
  }

  return(EXIT_SUCCESS);
}


/* helper functions */

void read_data(char *filename,gsl_vector_int *t_vec,
	       gsl_vector_int *i_vec,gsl_vector_int *j_vec,
	       gsl_vector *temp_vec)
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
  for (i=0;i<NUMOBS;i++) {
    fscanf(s,"%d",&lat);
    fscanf(s,"%d",&lon);
    fscanf(s,"%d",&year);
    fscanf(s,"%lf",&temp);
    gsl_vector_int_set(t_vec,i,year-1955);
    gsl_vector_int_set(i_vec,i,(lat-20)/2);
    gsl_vector_int_set(j_vec,i,(lon+30)/2);
    gsl_vector_set(temp_vec,i,temp);
  }
  fclose(s);
}

void add_obs(Spat * spat,gsl_vector_int *t_vec,
	     gsl_vector_int *i_vec,gsl_vector_int *j_vec,
	     gsl_vector *temp_vec)
{
  long i;
  for (i=0;i<NUMOBS;i++) {
    spat_add_obs(spat,
		 gsl_vector_int_get(t_vec,i),
		 gsl_vector_int_get(i_vec,i),
		 gsl_vector_int_get(j_vec,i),
		 gsl_vector_get(temp_vec,i)
		 );
  }
}



/* eof */

