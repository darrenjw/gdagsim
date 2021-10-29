/*
GDAGsim 0.3
(C) 2000-2002, Darren J Wilkinson

example_dlm.c
Example code for the GDAGsim library

Build and run with:

make example_dlm
./example_dlm 5000 > mcmc_output.tab

*/

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>

#include <gdag.h>

void add_latent(gdag *,int,double);
void add_obs(gdag *,gsl_vector *,double);

int main(int argc, char *argv[])
{
  gdag * d;
  gdag_usv * alpha;
  gsl_vector * s,* data, * scratch;
  gsl_vector_view vector_view;
  gsl_rng * r;
  int i,n;
  long it,iters;
  double state_p,obs_p,state_p_true,obs_p_true,sumsq;
  if (argc != 2) {
    printf("Usage: %s <iters>\n",argv[0]);
    exit(1);
  }
  iters=atoi(argv[1]);
  n=400;           /* number of time series observations */

  state_p_true=0.1; /* true precisions */
  obs_p_true=0.05;

  state_p=1.0;      /* sampled precisions */  
  obs_p=1.0;

  r=gsl_rng_alloc(gsl_rng_mt19937);
  
  /* Sampling from prior */
  d=gdag_calloc(2*n);
  add_latent(d,n,state_p_true);
  for (i=0;i<n;i++) {
    alpha=gdag_usv_basis(i);
    gdag_add_node(d,n+i,alpha,0.0,obs_p_true);
    gdag_usv_free(alpha);
  }
  gdag_process(d);
  s=gdag_sim(r,d);
  vector_view=gsl_vector_subvector(s,n,n);
  data=gsl_vector_alloc(n);
  gsl_vector_memcpy(data,&(vector_view.vector));
  gdag_free(d);
  /* Prior sampled */
  /* data stored in a vector called "data" */

  /* MCMC Loop... */
  scratch=gsl_vector_alloc(n);
  printf("Iter state_p obs_p\n");
  d=gdag_alloc(n);
  for (it=1;it<=iters;it++) {
    fprintf(stderr,"%ld ",it);
    gdag_set_zero(d);
    add_latent(d,n,state_p);
    add_obs(d,data,obs_p);
    gdag_process(d);
    s=gdag_sim(r,d);
    gsl_vector_memcpy(scratch,s);
    gsl_vector_sub(scratch,data);
    gsl_blas_ddot(scratch,scratch,&sumsq);
    obs_p=gdag_ran_gamma(r,0.1+0.5*n,0.1+0.5*sumsq);
    gsl_vector_memcpy(scratch,s);
    gdag_vector_diff(scratch);
    gsl_blas_ddot(scratch,scratch,&sumsq);
    state_p=gdag_ran_gamma(r,0.1+0.5*(n-1),0.1+0.5*sumsq);
    printf("%ld %f %f\n",it,state_p,obs_p);
  }
  fprintf(stderr,"\n");
  return(0);
}


void add_latent(gdag * d,int n,double prec)
{
  gdag_usv * alpha;
  int i;
  gdag_add_root(d,0,0.0,1.0);
  for (i=1;i<n;i++) {
    alpha=gdag_usv_basis(i-1);
    gdag_add_node(d,i,alpha,0.0,prec);
    gdag_usv_free(alpha);
  } 
}

void add_obs(gdag * d,gsl_vector * data,double prec)
{
  gdag_usv * alpha;
  int i;
  for (i=0;i < (data->size);i++) {
    alpha=gdag_usv_basis(i);
    gdag_add_observation(d,alpha,0.0,prec,gsl_vector_get(data,i));
    gdag_usv_free(alpha);
  }

}



/* eof */

