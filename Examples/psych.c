/*
psych.c

*/

#include <math.h>
#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <gdag.h>

/* global variables */
int n;
gsl_matrix *y;

/* helper function declarations */
gsl_matrix * read_data(void);
void add_latent(gdag * d,double mu,double sig0);
void add_obs(gdag * d,double sigma,double beta);

/* main program */
int main(int argc,char *argv[])
{
  gdag *d;
  gsl_rng *r;
  gsl_vector *s;
  long it,iters;
  double mu,sig0,mll,omll,a;
  double sigma,beta,sigma_p,beta_p;
  n=100;
  iters=10000;
  mu=1.0;
  sig0=100.0;
  r=gsl_rng_alloc(gsl_rng_mt19937);
  y=read_data();
  d=gdag_alloc(n+1);
  sigma=3.0;
  beta=2.5;
  omll=-1e100;
  printf("Iter beta sigma alpha x1\n");
  for (it=0;it<iters;it++) {
    /* update sigma and beta */
    sigma_p=exp(log(sigma)+gsl_ran_gaussian(r,0.1));
    beta_p=beta+gsl_ran_gaussian(r,0.1);
    gdag_set_zero(d);
    add_latent(d,mu,sig0);
    gdag_prior_process(d);
    add_obs(d,sigma_p,beta_p);
    gdag_process(d);
    mll=gdag_mloglik(d);
    a=mll-omll;
    if (gdag_accept_lp(r,a)) {
      sigma=sigma_p;
      beta=beta_p;
      omll=mll;
    }
    /* now simulate alpha and x */
    gdag_set_zero(d);
    add_latent(d,mu,sig0);
    add_obs(d,sigma,beta);
    gdag_process(d);
    s=gdag_sim(r,d);
    printf("%ld %f %f %f %f\n",it,beta,sigma,
	   gsl_vector_get(s,0),gsl_vector_get(s,1));
  }
  return(EXIT_SUCCESS);
}

/* helper functions */

void add_obs(gdag * d,double sigma,double beta)
{
  int i,j;
  gdag_usv *sv;
  for (i=0;i<n;i++) {
    for (j=0;j<3;j++) {
      sv=gdag_usv_alloc(2);
      gdag_usv_add(sv,0,1.0);
      gdag_usv_add(sv,i+1,beta);
      gdag_add_observation(d,sv,0.0,1.0/(sigma*sigma),
			   gsl_matrix_get(y,i,j));
      gdag_usv_free(sv);
    }
  }
}

void add_latent(gdag * d,double mu,double sig0)
{
  int i;
  gdag_add_root(d,0,mu,1.0/sig0);
  for (i=0;i<n;i++) {
    gdag_add_root(d,i+1,0,1);
  }
}

gsl_matrix * read_data(void)
{
  FILE *s;
  int i,j;
  float x;
  y=gsl_matrix_calloc(n,3);
  s=fopen("ydat.txt","r");
  if (s==NULL) {
    perror("failed to read ydat.txt");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<n;i++) {
    for (j=0;j<3;j++) {
      fscanf(s,"%f",&x);
      gsl_matrix_set(y,i,j,(double) x);
    }
  }
  fclose(s);
  return(y);
}


/* eof */

