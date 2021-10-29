/*
spat.h

Header file for "spat" (spatio-temporal modelling with GDAGsim)

(C) 2002, Darren J Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

*/

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>

#include <gdag.h>

/* main structure declaration */

typedef struct
{
  size_t T;          /* time points */
  size_t n;          /* rows */
  size_t m;          /* cols */
  double mu;      /* mean */
  double alpha;   /* autoregressive parameter */
  double tauS;    /* state precision */
  double tauO;    /* observational precision */
  gdag * gdso;    /* underlying GDAGsim object */
} Spat;

/* external function prototypes */

Spat * spat_calloc(size_t T,size_t n,size_t m);
void spat_free(Spat *self);
void spat_set_latent(Spat *self,double mu,double alpha,
		       double tauS,double tauO);
void spat_set_zero(Spat *self);
void spat_add_obs(Spat *self,size_t t,size_t i,size_t j,double obs);
void spat_process(Spat *self);
void spat_prior_process(Spat *self);
double spat_mloglik(Spat *self);
double spat_mean(Spat *self,size_t t,size_t i,size_t j);
double spat_var(Spat *self,size_t t,size_t i,size_t j);

/* eof */

