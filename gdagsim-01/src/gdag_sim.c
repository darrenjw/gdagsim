/*
GDAGsim 0.1
(C) 2000, Darren J Wilkinson

gdag_sim.c
Source file for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_blas.h>

#include <math.h>

#include <gdag.h>

static gsl_rng * gdag_r;   /* random number stream */

/* Simulation functions */

gsl_vector * gdag_sim(gdag * d)
{
  if (d->status < 3) {
    printf("Error: Status mismatch in gdag_sim.\n");
    exit(1);
  }
  gdag_vector_set_znorm(d->vv);
  gsl_blas_ddot(d->vv,d->vv,&d->ll);
  gsl_blas_dtrsv(CblasUpper,CblasNoTrans,CblasNonUnit,d->m,d->vv);
  gsl_blas_daxpy(1.0,d->v,d->vv);
  d->status = 4;
  return(d->vv);
}

void gdag_rng_init(void)
{
  gdag_r=gsl_rng_alloc(gsl_rng_mt19937);
}

void gdag_vector_set_znorm(gsl_vector * v)
{
  int i;
  for (i=0;i < v->size;i++) {
    gsl_vector_set(v,i,gsl_ran_gaussian(gdag_r,1));
  }
}

double gdag_ran_gamma(double a,double b)
{
  return(gsl_ran_gamma(gdag_r,a,1/b));
}

int gdag_accept_p(double p)
{
  return(gsl_ran_flat(gdag_r,0,1) < p);
}

int gdag_accept_lp(double lp)
{
  return(log(gsl_ran_flat(gdag_r,0,1)) < lp);
}




/* eof */

