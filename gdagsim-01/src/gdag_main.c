/*
GDAGsim 0.1
(C) 2000, Darren J Wilkinson

gdag_main.c
Main source file for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <math.h>

#include <gdag.h>


/* Allocation functions */

gdag * gdag_alloc(size_t n)
{
  gdag * d;
  d=malloc(sizeof(gdag));
  d->n=n;
  d->status=-1;
  d->v=gsl_vector_alloc(n);
  d->vv=gsl_vector_alloc(n);
  d->m=gsl_matrix_alloc(n,n);
  d->mm=gsl_matrix_alloc(n,n);
  return(d);
}

void gdag_set_zero(gdag * d)
{
  gsl_vector_set_zero(d->v);
  gsl_vector_set_zero(d->vv);
  gsl_matrix_set_zero(d->m);
  gsl_matrix_set_zero(d->mm);
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
  gsl_matrix_free(d->m);
  gsl_matrix_free(d->mm);
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
  gsl_matrix_set(d->m,i,i,precision);
  d->status = 1;
  d->count++;
  if (d->count == d->n) {
    d->status = 2;
  }
}

void gdag_add_node(gdag * d,size_t n,gdag_usv * alpha,double b,double precision)
{
  gsl_vector vec;
  if ((d->status < 0) || (d->status > 1)) {
    printf("Error: Status mismatch in gdag_add_node.\n");
    exit(1);
  }
  /* update precision matrix */
  gdag_dusger(precision,alpha,alpha,d->m);
  gsl_matrix_set(d->m,n,n,precision);
  vec=gsl_matrix_row(d->m,n);
  gdag_dusaxpy(-precision,alpha,&vec);
  vec=gsl_matrix_column(d->m,n);
  gdag_dusaxpy(-precision,alpha,&vec);
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
  gsl_linalg_cholesky_decomp(d->m);
  gsl_linalg_cholesky_svx(d->m,d->v);
  /* compute d->l */
  for (i=0;i < (d->n);i++) {
    d->l += log(gsl_matrix_get(d->m,i,i));
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
  gsl_matrix_memcpy(d->mm,d->m);
  gsl_linalg_cholesky_decomp(d->mm);
  gsl_linalg_cholesky_svx(d->mm,d->vv);
  gsl_blas_ddot(d->v,d->vv,&(d->pl));
  d->pl = -0.5*(d->pl);
  for (i=0;i < (d->n);i++) {
    d->pl += log(gsl_matrix_get(d->mm,i,i));
  }
}


/* eof */
