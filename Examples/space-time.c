/* now interested in generating the structure for AR1 process through time.
Using an index function which is indaex(i,j), i=location in space, j is location in time, and the function returns the associated point on the lattice */

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <math.h>
#include <gdag.h>

void add_latent(gdag *, int,int ,double, float);
void add_obs(gdag *, double, float);
int index(int,int);
int points;

int main(int argc, char *argv[])
{
  gdag * d;
  gsl_vector * s;
  int m,p,i;
  double state_p_true,obs_p_true;
  float a;
  
  gdag_rng_init();
  state_p_true=1; /* true precisions */
  obs_p_true=10;
  
  a=0.3;            /* value of alpha */
  
  points=15;        /* number of spatial obs */
  m=20;             /* number of time obs */
  p=points*m;       /* the size of grid */
 
  /* set up the latent structure */
  d=gdag_calloc(p);
  add_latent(d,m,points,state_p_true,a);
  
  /*add a few observations to the stucture*/
   add_obs(d,obs_p_true,a);  

  gdag_process(d);
  s=gdag_sim(d);
  for (i=0;i<p;i++)
    {
      printf("%3.4f     ",gsl_vector_get(s,i));
    }
  printf("\n");
  return(0);
}

int index(int i,int j)   
{ 
  int index;
  index=j*points+i;
  return(index);   
}  

void add_latent(gdag * d,int m, int points,double prec,float a)
{
  int i,j,x;
  gdag_usv * alpha;

  for(j=0;j<m;j++)
    {
      for(i=0;i<points;i++)
	{
	  x=index(i,j);

          if(x==0)                 /* root node*/
	    {
	      gdag_add_root(d,0,0.0,(1-a*prec));
	    }
	 	  
	  else if(j==0)            /* time=0 */
	    { 
	      alpha=gdag_usv_alloc(1);
	      gdag_usv_add(alpha,x-1,a);
	      gdag_add_node(d,x,alpha,0,prec);
	      gdag_usv_free(alpha);
	    }
	  
	  else if(i==0)            /* first column */
	    {   
	      alpha=gdag_usv_alloc(2); 
	      gdag_usv_add(alpha,(x-points+1),a);
	      gdag_usv_add(alpha,(x-points),a); 
	      gdag_add_node(d,x,alpha,0,prec); 
	      gdag_usv_free(alpha); 
	    }

	  
	  else if(i==(points-1))   /* last column */
	    {
	      alpha=gdag_usv_alloc(2);
	      gdag_usv_add(alpha,x-points,a);
	      gdag_usv_add(alpha,x-points-1,a);
	      gdag_add_node(d,x,alpha,0,prec);
	      gdag_usv_free(alpha);
	    }
	  
	  else                     /* middle */
	    {
	      alpha=gdag_usv_alloc(3);
	      gdag_usv_add(alpha,(x-points+1),a);
	      gdag_usv_add(alpha,(x-points-1),a);
	      gdag_usv_add(alpha,(x-points),a);
	      gdag_add_node(d,x,alpha,0,prec);
	      gdag_usv_free(alpha);
	    }
	}
    }
}

void add_obs(gdag *d,double prec, float a)
{
  gdag_usv * alpha;
   
  alpha=gdag_usv_basis(index(8,10));
  gdag_add_observation(d,alpha,0.0,prec,100);
  gdag_usv_free(alpha);

  alpha=gdag_usv_basis(index(4,5));
  gdag_add_observation(d,alpha,0.0,prec,0);
  gdag_usv_free(alpha);

  alpha=gdag_usv_basis(index(12,5));
  gdag_add_observation(d,alpha,0.0,prec,0);
  gdag_usv_free(alpha);

  alpha=gdag_usv_basis(index(4,15));
  gdag_add_observation(d,alpha,0.0,prec,0);
  gdag_usv_free(alpha);

  alpha=gdag_usv_basis(index(12,15));
  gdag_add_observation(d,alpha,0.0,prec,0);
  gdag_usv_free(alpha);


}

/* eof */
