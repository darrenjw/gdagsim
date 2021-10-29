/*
int-map.c

test code for inferring ingrated physical maps from multiple linkage maps

25/4/2005

*/

#include <math.h>
#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <gdag.h>

#define GVIG gsl_vector_int_get
#define GVG gsl_vector_get
#define GVS gsl_vector_set

/* global variables */
int n,m;
gdag_usv **y;  /* data is an array of sparse vectors */

/* helper function declarations */
void output_header(void);
void update_parameters(gsl_rng *r,gsl_vector *p,gsl_vector *pp);
void add_latent(gdag *d,gsl_vector *p);
void read_data(void);
void add_obs(gdag * d,gsl_vector *p);
void output_state(long it,gsl_vector *p,gsl_vector *s);
int in_range(gsl_vector *p);
void init_param(gsl_vector *p);


/* main program */
int main(int argc,char *argv[])
{
  gdag *d;
  gsl_rng *r;
  gsl_vector *s,*param,*param_prop;
  double mll,omll,a;
  long it,iters;

  n=53; m=2; /* distinct markers and number of maps */
  if (argc!=2) {
    fprintf(stderr,"Wrong number of args.\n");
    return(EXIT_FAILURE);
  }
  iters=atoi(argv[1]);
  param=gsl_vector_calloc(2*m);
  param_prop=gsl_vector_calloc(2*m);
  r=gsl_rng_alloc(gsl_rng_mt19937);
  d=gdag_alloc(n+m);
  omll=-1e100;
  init_param(param);
  read_data();
  output_header();
  for (it=0;it<iters;it++) {
    update_parameters(r,param,param_prop);
    if (in_range(param_prop)) {
      gdag_set_zero(d);
      add_latent(d,param_prop);
      gdag_prior_process(d);
      add_obs(d,param_prop);
      gdag_process(d);
      mll=gdag_mloglik(d);
      a=mll-omll;
      if (gdag_accept_lp(r,a)) {
	gsl_vector_memcpy(param,param_prop);
	omll=mll;
      }
    }
    /* now simulate the mu */
    gdag_set_zero(d);
    add_latent(d,param);
    add_obs(d,param);
    gdag_process(d);
    s=gdag_sim(r,d);
    output_state(it,param,s);
  }
  return(EXIT_SUCCESS);
}

/* helper functions */

void update_parameters(gsl_rng *r,gsl_vector *p,gsl_vector *pp)
{
  int i;
  for (i=0;i<m;i++) {
    GVS(pp,i,GVG(p,i)+gsl_ran_gaussian(r,0.1));
    GVS(pp,m+i,GVG(p,m+i)+gsl_ran_gaussian(r,0.2));
  }
}

void init_param(gsl_vector *p)
{
  int i;
  for (i=0;i<m;i++)
    GVS(p,i,0); /* nu_i init */
  for (i=m;i<2*m;i++)
    GVS(p,i,2); /* lambda_i init */
}

int in_range(gsl_vector *p)
{
  int i,ok;
  ok=1;
  /* finite uniform priors on the parameters */
  for (i=0;i<m;i++) {
    if ( (GVG(p,i)<-1) || (GVG(p,i)>2) ) /* nu_i prior */
      ok=0;
  }
  for (i=m;i<2*m;i++) {
    if ( (GVG(p,i)<1) || (GVG(p,i)>4) ) /* lambda_i prior */
      ok=0;
  }
  return ok;
}

void add_latent(gdag * d,gsl_vector *p)
{
  int i;
  for (i=0;i<n;i++) {
    if (i==9)
      gdag_add_root(d,i,50,100); /* marker 9 around 50 */
    else if (i==31) 
      gdag_add_root(d,i,120,100); /* marker 31 around 120 */
    else
      gdag_add_root(d,i,100,0.0001); /* prior for mu_i */
  }
  for (i=n;i<m+n;i++) {
    gdag_add_root(d,i,0,0.0001); /* prior for beta_i */
  }
}

void output_state(long it,gsl_vector *p,gsl_vector *s)
{
  int i;
  printf("%ld",it);
  for (i=0;i<m;i++)
    printf(" %f",GVG(p,i));
  for (i=0;i<m;i++)
    printf(" %f",GVG(p,m+i));
  for (i=0;i<m;i++)
    printf(" %f",GVG(s,n+i));
  for (i=0;i<n;i++)
    printf(" %f",GVG(s,i));
  printf("\n");
}

void output_header(void)
{
  int i;
  printf("Iter");
  for (i=0;i<m;i++)
    printf(" nu[%d]",i);
  for (i=0;i<m;i++)
    printf(" lambda[%d]",i);
  for (i=0;i<m;i++)
    printf(" beta[%d]",i);
  for (i=0;i<n;i++)
    printf(" mu[%d]",i);
  printf("\n");
}

void add_obs(gdag * d,gsl_vector *p)
{
  int i,j;
  gdag_usv *sv;
  for (i=0;i<m;i++) {
    for (j=0;j<(y[i]->nz);j++) {
      sv=gdag_usv_alloc(2);
      gdag_usv_add(sv,GVIG(y[i]->indx,j),exp(-2*GVG(p,i)));
      gdag_usv_add(sv,n+i,1.0);
      gdag_add_observation(d,sv,0.0,exp(-2*GVG(p,m+i)),
			   GVG(y[i]->x,j));
      gdag_usv_free(sv);
    }
  }
}

void read_data(void)
{
  y=malloc(m*sizeof(gdag_usv *));
  y[0]=gdag_usv_alloc(32);
  gdag_usv_add(y[0],1,0.0);
  gdag_usv_add(y[0],2,3.3);
  gdag_usv_add(y[0],3,4.4);
  gdag_usv_add(y[0],4,24.4);
  gdag_usv_add(y[0],5,34.2);
  gdag_usv_add(y[0],6,38.8);
  gdag_usv_add(y[0],7,43.8);
  gdag_usv_add(y[0],8,43.8);
  gdag_usv_add(y[0],9,48.3);
  gdag_usv_add(y[0],10,48.3);
  gdag_usv_add(y[0],11,59.8);
  gdag_usv_add(y[0],12,62.2);
  gdag_usv_add(y[0],13,62.2);
  gdag_usv_add(y[0],14,65.5);
  gdag_usv_add(y[0],15,68.8);
  gdag_usv_add(y[0],16,71.8);
  gdag_usv_add(y[0],17,76.5);
  gdag_usv_add(y[0],18,76.5);
  gdag_usv_add(y[0],19,76.5);
  gdag_usv_add(y[0],20,79.9);
  gdag_usv_add(y[0],21,81.5);
  gdag_usv_add(y[0],22,83.0);
  gdag_usv_add(y[0],23,90.7);
  gdag_usv_add(y[0],24,95.9);
  gdag_usv_add(y[0],25,96.5);
  gdag_usv_add(y[0],26,97.6);
  gdag_usv_add(y[0],27,98.3);
  gdag_usv_add(y[0],28,99.3);
  gdag_usv_add(y[0],29,100.3);
  gdag_usv_add(y[0],30,104.6);
  gdag_usv_add(y[0],31,110.1);
  gdag_usv_add(y[0],32,120.8);
  y[1]=gdag_usv_alloc(29);
  gdag_usv_add(y[1],33,0.0);
  gdag_usv_add(y[1],34,3.7);
  gdag_usv_add(y[1],35,10.4);
  gdag_usv_add(y[1],36,26.0);
  gdag_usv_add(y[1],37,61.9);
  gdag_usv_add(y[1],38,80.6);
  gdag_usv_add(y[1],39,86.1);
  gdag_usv_add(y[1],9,99.8);
  gdag_usv_add(y[1],40,106.2);
  gdag_usv_add(y[1],12,110.7);
  gdag_usv_add(y[1],41,112.5);
  gdag_usv_add(y[1],42,113.9);
  gdag_usv_add(y[1],43,114.7);
  gdag_usv_add(y[1],11,115.8);
  gdag_usv_add(y[1],44,116.8);
  gdag_usv_add(y[1],45,117.2);
  gdag_usv_add(y[1],46,120.2);
  gdag_usv_add(y[1],47,126.4);
  gdag_usv_add(y[1],48,130.0);
  gdag_usv_add(y[1],16,132.7);
  gdag_usv_add(y[1],23,157.4);
  gdag_usv_add(y[1],29,168.4);
  gdag_usv_add(y[1],32,179.8);
  gdag_usv_add(y[1],31,184.4);
  gdag_usv_add(y[1],49,188.4);
  gdag_usv_add(y[1],50,226.7);
  gdag_usv_add(y[1],51,231.0);
  gdag_usv_add(y[1],52,238.5);
  gdag_usv_add(y[1],0,244.2);
}


void read_data_test(void)
{
  y=malloc(m*sizeof(gdag_usv *));
  y[0]=gdag_usv_alloc(9);
  gdag_usv_add(y[0],0,0);
  gdag_usv_add(y[0],2,2);
  gdag_usv_add(y[0],3,4);
  gdag_usv_add(y[0],5,6);
  gdag_usv_add(y[0],7,8);
  gdag_usv_add(y[0],8,9);
  gdag_usv_add(y[0],9,11);
  gdag_usv_add(y[0],11,13);
  gdag_usv_add(y[0],13,15);
  y[1]=gdag_usv_alloc(10);
  gdag_usv_add(y[1],1,0);
  gdag_usv_add(y[1],3,2);
  gdag_usv_add(y[1],4,3);
  gdag_usv_add(y[1],6,4);
  gdag_usv_add(y[1],7,6);
  gdag_usv_add(y[1],8,5);
  gdag_usv_add(y[1],10,8);
  gdag_usv_add(y[1],11,10);
  gdag_usv_add(y[1],12,12);
  gdag_usv_add(y[1],14,13);
}



/* OLD FUNCTIONS */

/*
gsl_matrix * read_data_old(void)
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
*/

/* eof */

