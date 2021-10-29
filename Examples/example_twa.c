/*
GDAGsim 0.1
(C) 2000, Darren J Wilkinson

example_twa.c
Example code for the GDAGsim library

Build and run with:

make example_twa
./example_twa 5000 > mcmc_output.tab

*/

#include <math.h>
#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <gdag.h>

/* global variables */
int p,q,**r;

/* helper function declarations */
void add_latent(gdag *,double,double,double,double,double);
void add_obs(gdag *,double ***,double);
gdag_usv * dep_usv(int,int);

/* data reading functions */
void read_inf(double *,double *,double *,double *
	      ,double *,double *,double *,double *, int *);
double *** read_data();
void output_data(double ***);


/* main program */
int main(int argc,char *argv[])
{
  gdag *d;
  gdag_usv *sv;
  gsl_vector *s;
  gsl_vector vec_view;
  int i,j,k, numrs;
  long it,iters;
  double ***y,sumsq,res;
  double tau_alpha,tau_beta,tau_gamma,tau_eps;
  double a_alpha,b_alpha,a_beta,b_beta,a_gamma,b_gamma,a_eps,b_eps;
  if (argc != 2) {
    printf("Usage: %s <iters>\n",argv[0]);
    exit(1);
  }
  iters=atoi(argv[1]);
  gdag_rng_init();
  read_inf(&a_alpha,&b_alpha,&a_beta,&b_beta,&a_gamma,&b_gamma
	   ,&a_eps,&b_eps, &numrs);
  tau_alpha=a_alpha/b_alpha;
  tau_beta=a_beta/b_beta;
  tau_gamma=a_gamma/b_gamma;
  tau_eps=a_eps/b_eps;
  y=read_data();

  /* MCMC Loop... */
  printf("Iter lv_alpha lv_beta lv_gamma lv_eps mu\n");
  d=gdag_alloc(1+p+q+p*q);
  for (it=1;it<=iters;it++) {
    fprintf(stderr,"%ld ",it);
    gdag_set_zero(d);
    add_latent(d,0.0,0.0,tau_alpha,tau_beta,tau_gamma); /* flat prior on mu */
    add_obs(d,y,tau_eps);
    gdag_process(d);
    s=gdag_sim(d);
    /* update tau_alpha */
    vec_view=gsl_vector_subvector(s,1,p);
    gsl_blas_ddot(&vec_view,&vec_view,&sumsq);
    tau_alpha=gdag_ran_gamma(a_alpha+0.5*p,b_alpha+0.5*sumsq);
    /* update tau_beta */
    vec_view=gsl_vector_subvector(s,p+1,q);
    gsl_blas_ddot(&vec_view,&vec_view,&sumsq);
    tau_beta=gdag_ran_gamma(a_beta+0.5*q,b_beta+0.5*sumsq);
    /* update tau_gamma */
    vec_view=gsl_vector_subvector(s,p+q+1,p*q);
    gsl_blas_ddot(&vec_view,&vec_view,&sumsq);
    tau_gamma=gdag_ran_gamma(a_gamma+0.5*p*q,b_gamma+0.5*sumsq);
    /* update tau_eps */
    sumsq=0.0;
    for (i=0;i<p;i++) {
      for (j=0;j<q;j++) {
	sv=dep_usv(i,j);
	gdag_dusdot(sv,s,&res);
	for (k=0;k<r[i][j];k++) {
	  sumsq+=gdag_sqr(y[i][j][k]-res);
	}
	gdag_usv_free(sv);
      }
    }
    tau_eps=gdag_ran_gamma(a_eps+0.5*numrs,b_eps+0.5*sumsq);
    /* output iteration */
    printf("%ld %f %f %f %f %f\n",it,log(1/tau_alpha),log(1/tau_beta),
	   log(1/tau_gamma),log(1/tau_eps),gsl_vector_get(s,0));
  }
  fprintf(stderr,"\n");
  return(0);
}


/* helper functions */
void add_latent(gdag * d,double mu,double tau,double tau_alpha,double tau_beta,double tau_gamma)
{
  int i;
  gdag_add_root(d,0,mu,tau);
  for (i=1;i<=p;i++) {
    gdag_add_root(d,i,0.0,tau_alpha);
  } 
  for (i=1;i<=q;i++) {
    gdag_add_root(d,p+i,0.0,tau_beta);
  }
  for (i=1;i<=(p*q);i++) {
    gdag_add_root(d,p+q+i,0.0,tau_gamma);
  }
}

void add_obs(gdag * d,double *** y,double prec)
{
  gdag_usv * sv;
  int i,j,k;
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      sv=dep_usv(i,j);
      for (k=0;k<r[i][j];k++) {
	gdag_add_observation(d,sv,0.0,prec,y[i][j][k]);
      }
      gdag_usv_free(sv);
    }
  }
}

gdag_usv * dep_usv(int i,int j)
{
  gdag_usv * sv;
	sv=gdag_usv_alloc(4);
	gdag_usv_add(sv,0,1.0);
	gdag_usv_add(sv,1+i,1.0);
	gdag_usv_add(sv,1+p+j,1.0);
	gdag_usv_add(sv,1+p+q+q*i+j,1.0);
	return(sv);
}




/* functions for reading in the data */

void read_inf(double *a_alpha,double *b_alpha,double *a_beta,double *b_beta
	      ,double *a_gamma,double *b_gamma,double *a_eps,double *b_eps, 
	      int *numrs)
{
  FILE *file;
  char *data;
  int i,j;
  *numrs = 0;
  data=malloc(100);
  if ((file=fopen("foo.inf","r")) == NULL) {
    printf("Cannot open foo.inf\n");
    exit(1);
  }	
  fscanf(file,"%s",data); p=atoi(data);
  fscanf(file,"%s",data); q=atoi(data);
  r=malloc(p*sizeof(int *));
  for (i=0;i<p;i++) {
    r[i]=malloc(q*sizeof(int));
    for (j=0;j<q;j++) {
      fscanf(file,"%s",data); r[i][j]=atoi(data);
      *numrs += r[i][j];
    }
  }
  fscanf(file,"%s",data); *a_alpha=atof(data);
  fscanf(file,"%s",data); *b_alpha=atof(data);
  fscanf(file,"%s",data); *a_beta=atof(data);
  fscanf(file,"%s",data); *b_beta=atof(data);
  fscanf(file,"%s",data); *a_gamma=atof(data);
  fscanf(file,"%s",data); *b_gamma=atof(data);
  fscanf(file,"%s",data); *a_eps=atof(data);
  fscanf(file,"%s",data); *b_eps=atof(data);
  fclose(file);
  free(data);
}

double *** read_data()
{
  FILE *file;
  char *data;
  int i,j,k;
  double ***y;
  data=malloc(100);
  if ((file=fopen("foo.dat","r")) == NULL) {
    printf("Cannot open foo.dat\n");
    exit(1);
  }
  y=malloc(p*sizeof(double **));
  for (i=0;i<p;i++) {
    y[i]=malloc(q*sizeof(double *));
    for (j=0;j<q;j++) {
      y[i][j]=malloc(r[i][j]*sizeof(double));
      for (k=0;k<r[i][j];k++) {
	fscanf(file,"%s",data); y[i][j][k]=atof(data);
      }
    }
  }
  fclose(file);
  free(data);
  return(y);
}

void output_data(double ***y)
{
  int i,j,k;
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      for (k=0;k<r[i][j];k++) {
	printf("%f ",y[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
}




/* eof */

