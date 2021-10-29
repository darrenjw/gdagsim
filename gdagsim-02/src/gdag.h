/*
GDAGsim 0.2
(C) 2000, Darren J Wilkinson

gdag.h
Header file for the GDAGsim library
*/

#include <sparse.h>
#include <sparse2.h>
#include <gdag_meschach.h>

typedef struct
{
  size_t n;         /*  num latent vars                */
  size_t count;     /*  num added so far               */
  size_t obs_count; /*  num obs added                  */
  int status;       /*  see below...                   */
  gsl_vector * v;   /*  main vector slot               */
  VEC * V;          /*  meschach view of v             */
  gsl_vector * vv;  /*  temporary vector               */
  VEC * VV;         /*  meschach view of vv            */
  SPMAT * m;        /*  main matrix slot               */
  SPMAT * mm;       /*  temporary matrix               */
  double l;         /*  (log|K'|)/2                    */
  double ll;        /*  z'z                            */
  double pl;        /*  (log|K|-h'K^{-1}h)/2           */
  double cl;        /*  (log|L|-(y-b)L(y-b))/2         */
} gdag;

/* 
Status:
-1 Allocated but not cleared
0 Allocated and cleared
1 Some node definitions, but not full
2 Complete node definitions
3 Processed
4 Sampled
*/

typedef struct
{
  int nz;                  /* number of non-zeros   */
  gsl_vector * x;          /* non-zero values       */
  gsl_vector_int * indx;   /* non-zero indices      */
  int count;               /* number added so far   */
} gdag_usv;

int gdag_accept_p(double);

int gdag_accept_lp(double);

void gdag_add_node(gdag *,size_t,gdag_usv *,double,double);

void gdag_add_observation(gdag *,gdag_usv *,double,double,double);

void gdag_add_root(gdag *,size_t,double, double);

gdag * gdag_alloc(size_t);

gdag * gdag_calloc(size_t);

SPMAT * gdag_chol(gdag *);

size_t gdag_count(gdag *);

void gdag_dump(gdag *);

void gdag_dusaxpy(double,gdag_usv *,gsl_vector *);

void gdag_row_update(double,gdag_usv *,SPMAT *,size_t);

void gdag_col_update(double,gdag_usv *,SPMAT *,size_t);

void gdag_dusdot(gdag_usv *,gsl_vector *,double *);

void gdag_dusger(double,gdag_usv *,gdag_usv *,SPMAT *);

void gdag_dussc(gdag_usv *,gsl_vector *);

double gdag_ex_sq(gdag *,gdag_usv *);

void gdag_free(gdag *);

double gdag_loglik(gdag *);

double gdag_mloglik(gdag *);

double gdag_vloglik(gdag *,gsl_vector *);

void gdag_matrix_dump(gsl_matrix *);

gsl_vector * gdag_mean(gdag *);

void gdag_prior_process(gdag *);

void gdag_process(gdag *);

double gdag_ran_gamma(double,double);

double gdag_ran_gaussian(double,double);

double gdag_ran_uniform(double,double);

void gdag_rng_init(void);

void gdag_set_zero(gdag *);

gsl_vector * gdag_sim(gdag *);

size_t gdag_size(gdag *);

double gdag_sqr(double);

int gdag_status(gdag *);

void gdag_usv_add(gdag_usv *,int,double);

gdag_usv * gdag_usv_alloc(int);

gdag_usv * gdag_usv_basis(int);

void gdag_usv_dump(gdag_usv *);

void gdag_usv_free(gdag_usv *);

double gdag_var(gdag *,gdag_usv *);

void gdag_vector_diff(gsl_vector *);

void gdag_vector_dump(gsl_vector *);

void gdag_vector_set_znorm(gsl_vector *);




/* log(2 pi)/2 */
#define GDAG_LTPOT 0.918938533205





/* eof */


