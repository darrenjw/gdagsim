/*
GDAGsim 0.1
(C) 2000, Darren J Wilkinson

gdag_utils.c
Source file for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gdag.h>


/* Accessor functions */

void gdag_dump(gdag * d)
{
  printf("n: %d\n",d->n);
  printf("count: %d\n",d->count);
  printf("status: %d\n",d->status);
  printf("v: ");
  gdag_vector_dump(d->v);
  printf("vv: ");
  gdag_vector_dump(d->vv);
  printf("m:\n");
  gdag_matrix_dump(d->m);
}

size_t gdag_size(gdag * d)
{
  return(d->n);
}

size_t gdag_count(gdag * d)
{
  return(d->count);
}

int gdag_status(gdag * d)
{
  return(d->status);
}

gsl_vector * gdag_mean(gdag * d)
{
  if (d->status < 3) {
    printf("Error: Status mismatch in gdag_mean.\n");
    exit(1);
  }
  return(d->v);
}

gsl_matrix * gdag_chol(gdag * d)
{
  if (d->status < 3) {
    printf("Error: Status mismatch in gdag_chol.\n");
    exit(1);
  }
  return(d->m);
}

double gdag_var(gdag * d,gdag_usv * alpha)
{
  double res;
  if (d->status < 3) {
    printf("Error: Status mismatch in gdag_var.\n");
    exit(1);
  }
  gsl_vector_set_zero(d->vv);
  gdag_dussc(alpha,d->vv);
  gsl_blas_dtrsv(CblasLower,CblasNoTrans,CblasNonUnit,d->m,d->vv);
  gsl_blas_ddot(d->vv,d->vv,&res);
  return(res);
}

double gdag_ex_sq(gdag * d,gdag_usv * alpha)
{
  double res;
  gdag_dusdot(alpha,gdag_mean(d),&res);
  return(res*res+gdag_var(d,alpha));
}

double gdag_loglik(gdag * d)
{
  if (d->status < 4) {
    printf("Error: Status mismatch in gdag_loglik.\n");
    exit(1);
  }
  return((d->l) - 0.5*(d->ll) - GDAG_LTPOT*(d->n) );
}  

double gdag_mloglik(gdag * d)
{
  double res;
  if (d->status < 3) {
    printf("Error: Status mismatch in gdag_mloglik.\n");
    exit(1);
  }
  /* gsl_blas_dgemv(CblasNoTrans,1.0,d->m,d->v,0.0,d->vv); */
  gsl_vector_memcpy(d->vv,d->v);
  gsl_blas_dtrmv(CblasUpper,CblasNoTrans,CblasNonUnit,d->m,d->vv);
  gsl_blas_ddot(d->vv,d->vv,&res);
  res=(d->l - 0.5*res);
  return( -GDAG_LTPOT*(d->obs_count) + (d->pl) + (d->cl) - res );
}

/* Helper functions */

void gdag_vector_dump(gsl_vector * v)
{
  size_t i;
  for (i=0;i< v->size;i++) {
    printf("%f ",gsl_vector_get(v,i));
  }
  printf("\n");
}

void gdag_matrix_dump(gsl_matrix * m)
{
  size_t i,j;
  for (i=0;i < m->size1;i++) {
    for (j=0;j < m->size2; j++) {
      printf("%f ",gsl_matrix_get(m,i,j));
    }
    printf("\n");
  }
  printf("\n");
}

void gdag_vector_diff(gsl_vector * v)
{
  int i,n;
  n=v->size;
  for (i=0;i<(n-1);i++) {
    gsl_vector_set(v,i,gsl_vector_get(v,i+1)-gsl_vector_get(v,i));
  }
  gsl_vector_set(v,n-1,0);
}

double gdag_sqr(double x)
{
  return(x*x);
}




/* eof */

