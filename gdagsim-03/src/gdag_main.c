/*
GDAGsim 0.3
(C) 2000-2002, Darren J Wilkinson

gdag_main.c
Main source file for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <gsl_rng.h>

#include <math.h>

#include <gdag.h>

#define GDAG_MES_SIZE 50

/* Allocation functions */

gdag * gdag_alloc(size_t n)
{
  gdag * d;
  d=malloc(sizeof(gdag));
  if (d == NULL) { perror("gdag_alloc d"); }
  d->n=n;
  d->status=-1;
  d->v=gsl_vector_alloc(n);
  if (d->v == NULL) { perror("gdag_alloc v"); }
  d->V=gdag_gsltomes(d->v);
  d->vv=gsl_vector_alloc(n);
  if (d->vv == NULL) { perror("gdag_alloc vv"); }
  d->VV=gdag_gsltomes(d->vv);
  d->m=sp_get(n,n,GDAG_MES_SIZE);
  if (d->m == NULL) { perror("gdag_alloc m"); }
  d->mm=sp_get(n,n,GDAG_MES_SIZE);
  if (d->mm == NULL) { perror("gdag_alloc mm"); }
  return(d);
}

void gdag_set_zero(gdag * d)
{
  gsl_vector_set_zero(d->v);
  gsl_vector_set_zero(d->vv);
  sp_zero(d->m);
  sp_zero(d->mm);
  d->count = 0;
  d->obs_count = 0;
  d->status = 0;
  d->l = 0.0;
  d->ll = 0.0;
  d->pl = 0.0;
  d->cl = 0.0;
}

gdag * gdag_calloc(size_t n)
{
  gdag * d;
  d=gdag_alloc(n);
  gdag_set_zero(d);
  return(d);
}

void gdag_free(gdag * d)
{
  gsl_vector_free(d->v);
  gsl_vector_free(d->vv);
  sp_free(d->m);
  sp_free(d->mm);
  free(d->V);
  free(d->VV);
  free(d);
}


/* Building functions */

void gdag_add_root(gdag * d,size_t i,double mean, double precision)
{
  if ((d->status < 0) || (d->status > 1)) {
    printf("Error: Status mismatch in gdag_add_root.\n");
    exit(1);
  }
  gsl_vector_set(d->v,i,(precision*mean));
  sp_set_val(d->m,i,i,precision);
  d->status = 1;
  d->count++;
  if (d->count == d->n) {
    d->status = 2;
  }
}

void gdag_add_node(gdag * d,size_t n,gdag_usv * alpha,double b,double precision)
{
  if ((d->status < 0) || (d->status > 1)) {
    printf("Error: Status mismatch in gdag_add_node.\n");
    exit(1);
  }
  /* update precision matrix */
  gdag_dusger(precision,alpha,alpha,d->m);
  sp_set_val(d->m,n,n,precision);
  gdag_row_update(-precision,alpha,d->m,n);
  gdag_col_update(-precision,alpha,d->m,n);
  /* update location vector */
  gdag_dusaxpy(-precision*b,alpha,d->v);
  gsl_vector_set(d->v,n,precision*b);
  /* update count */
  d->count++;
  d->status = 1;
  if (d->count == d->n) {
    d->status = 2;
  }
}

void gdag_add_observation(gdag * d,gdag_usv * alpha,double b,double precision,double obs)
{
  if ((d->status < 1)||(d->status >2)) {
    printf("Error: Status mismatch in gdag_add_observation.\n");
    exit(1);
  }
  /* update structure */
  gdag_dusger(precision,alpha,alpha,d->m);
  gdag_dusaxpy(precision*(obs-b),alpha,d->v);
  /* update d->cl */
  d->cl += 0.5*(log(precision) - precision*(obs-b)*(obs-b));
  /* update count */
  d->obs_count++;
}

/* Processing functions */

void gdag_process(gdag * d)
{
  int i;
  if (d->status != 2) {
    printf("Error: Status mismatch in gdag_process.\n");
    exit(1);
  }
  /* decompose structure */
  spCHfactor(d->m);
  (d->V)=spCHsolve(d->m,d->V,d->V);
  /* compute d->l */
  for (i=0;i < (d->n);i++) {
    d->l += log(sp_get_val(d->m,i,i));
  }
  d->status = 3;
}

void gdag_prior_process(gdag * d)
{
  int i;
  if (d->status != 2) {
    printf("Error: Status mismatch in gdag_prior_process.\n");
    exit(1);
  }
  /* compute d->pl */
  gsl_vector_memcpy(d->vv,d->v);
  sp_copy2(d->m,d->mm);
  spCHfactor(d->mm);
  (d->VV)=spCHsolve(d->mm,d->VV,d->VV);
  gsl_blas_ddot(d->v,d->vv,&(d->pl));
  d->pl = -0.5*(d->pl);
  for (i=0;i < (d->n);i++) {
    d->pl += log(sp_get_val(d->mm,i,i));
  }
}


/* eof */
