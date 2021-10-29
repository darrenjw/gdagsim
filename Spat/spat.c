/*
spat.c

Code file for "spat" (spatio-temporal modelling with GDAGsim)

(C) 2002, Darren J Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

*/

#include "spat.h"

/* internal function prototypes */

size_t spat_index(Spat *self,size_t t,size_t i,size_t j);


/* function definitions */

Spat * spat_calloc(size_t T,size_t n,size_t m)
{
  Spat *self;
  self=malloc(sizeof(Spat));
  self->T = T;
  self->n = n;
  self->m = m;
  self->mu = 0.0;
  self->tauS = 0.0;
  self->tauO = 0.0;
  self->gdso = gdag_calloc(n*m*T);
  return(self);
}

void spat_free(Spat *self)
{
  gdag_free(self->gdso);
  free(self);
}

void spat_set_zero(Spat *self)
{
  self->mu = 0.0;
  self->alpha = 0.0;
  self->tauS = 0.0;
  self->tauO = 0.0;
  gdag_set_zero(self->gdso);
}

size_t spat_index(Spat *self,size_t t,size_t i,size_t j)
{
  return( t*(self->n)*(self->m) + i*(self->m) + j );
}

void spat_add_obs(Spat *self,size_t t,size_t i,size_t j,double obs)
{
  gdag_usv *sp=gdag_usv_basis(spat_index(self,t,i,j));
  gdag_add_observation(self->gdso,sp,0.0,self->tauO,obs);
  gdag_usv_free(sp);
}

void spat_process(Spat *self)
{
  gdag_process(self->gdso);
}

void spat_prior_process(Spat *self)
{
  gdag_prior_process(self->gdso);
}

double spat_mloglik(Spat *self)
{
  return(gdag_mloglik(self->gdso));
}

double spat_mean(Spat *self,size_t t,size_t i,size_t j)
{
  gsl_vector *v=gdag_mean(self->gdso);
  return( gsl_vector_get(v,spat_index(self,t,i,j)) );
}

double spat_var(Spat *self,size_t t,size_t i,size_t j)
{
  double var;
  gdag_usv *sv=gdag_usv_basis(spat_index(self,t,i,j));
  var=gdag_var(self->gdso,sv);
  gdag_usv_free(sv);
  return(var);
}

void spat_set_latent(Spat *self,double mu,double alpha,
		       double tauS,double tauO)
{
  size_t t,i,j;
  gdag_usv *sv;
  double b;
  t=0;
  b=(1.0-alpha)*mu;
  /* record parameters in the object for future use */
  self->mu = mu;
  self->alpha = alpha;
  self->tauS = tauS;
  self->tauO = tauO;
  /* time zero layer */
  for (i=0 ; i < self->n ; i++) {
    for (j=0 ; j < self->m ; j++) {
      gdag_add_root(self->gdso,spat_index(self,t,i,j),mu,tauS);
    }
  }
  /* all other times */
  for (t=1 ; t < self->T ; t++) {
    for (i=0 ; i < self->n ; i++) {
      for (j=0 ; j < self->m ; j++) {
	/* corners */
	if ( (i==0)&&(j==0) ) {
	  sv=gdag_usv_alloc(4);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j+1),alpha*1.0/9);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	} 
	else if ( (i==0)&&(j==(self->m)-1) ) {
	  sv=gdag_usv_alloc(4);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j-1),alpha*1.0/9);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	} 
	else if ( (i==(self->n)-1)&&(j==0) ) {
	  sv=gdag_usv_alloc(4);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j+1),alpha*1.0/9);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	} 
	else if ( (i==(self->n)-1)&&(j==(self->m)-1) ) {
	  sv=gdag_usv_alloc(4);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*2.0/9);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j-1),alpha*1.0/9);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);	  
	}
	/* edges */
	else if (i==0) {
	  sv=gdag_usv_alloc(6);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j+1),alpha*1.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j-1),alpha*1.0/12);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	}
	else if (j==0) {
	  sv=gdag_usv_alloc(6);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j+1),alpha*1.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j+1),alpha*1.0/12);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	}
	else if (i==(self->n)-1) {
	  sv=gdag_usv_alloc(6);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j+1),alpha*1.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j-1),alpha*1.0/12);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	}
	else if (j==(self->m)-1) {
	  sv=gdag_usv_alloc(6);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*4.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*2.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j-1),alpha*1.0/12);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j-1),alpha*1.0/12);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);
	}
	/* interior */
	else {
	  sv=gdag_usv_alloc(9);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j),alpha*0.25);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j),alpha*0.125);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j),alpha*0.125);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j+1),alpha*0.125);
	  gdag_usv_add(sv,spat_index(self,t-1,i,j-1),alpha*0.125);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j+1),alpha*0.0625);
	  gdag_usv_add(sv,spat_index(self,t-1,i+1,j-1),alpha*0.0625);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j+1),alpha*0.0625);
	  gdag_usv_add(sv,spat_index(self,t-1,i-1,j-1),alpha*0.0625);
	  gdag_add_node(self->gdso,spat_index(self,t,i,j),sv,b,tauS);
	  gdag_usv_free(sv);	  
	}
      }
    }
  }
}



/* eof */
