/* set up the prior distribution without using gdag. n*m spatial grid and T time points, for use when we need a FULL dataset*/

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <math.h>

int main(int argc, char *argv[])
{
  gsl_matrix *state, * newstate, * observations;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_mt19937);
  FILE * file;
  double state_p, obs_p, a;
  int n, m, T, i, j, t;
  if (argc != 4) 
    {
      printf("Usage: %s <n> <m> <T> \n",argv[0]);
      exit(1);
    }

  n=atoi(argv[1]);
  m=atoi(argv[2]);
  T=atoi(argv[3]); /* if needed could also have it so a, state_p & obs_p are inputted*/

  if ((file=fopen("prior.dat","w")) == NULL) 
    {
      printf("Cannot open prior.dat\n");
      exit(1);
    }	
  
  state_p=0.1; /* true precisions */
  obs_p=0.05; 
  a=0.7;       /* true alpha*/
 
  state=gsl_matrix_calloc(n,m);
  newstate=gsl_matrix_calloc(n,m);
  observations=gsl_matrix_calloc(n,m);

  for(t=0;t<T;t++)
    {
      if (t==0)
	{
	  for (j=0;j<m;j++)
	    {
	      for (i=0;i<n;i++)
		{
		  gsl_matrix_set(state,i,j,gsl_ran_gaussian(r,sqrt(1/state_p)));
		}
	    }

	  gsl_matrix_memcpy(observations,state); 
	  
	  for(i=0;i<n;i++)
	    {
	      for (j=0;j<m;j++)
		{
		  gsl_matrix_set(observations,i,j,gsl_matrix_get(observations,i,j)+gsl_ran_gaussian(r,sqrt(1/obs_p)));
		}
	    }
	}  
    
      else
	{
	  for(j=0;j<m;j++)
	    {
	      for(i=0;i<n;i++)
		{
		  /*create the "newstate"*/
		  if ((i==0)&(j==0))
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i,j+1))+(1/16)*gsl_matrix_get(state,i+1,j+1)));
		    }
		  
		  else if ((i==0)&(j==m-1))
		    {
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i,j-1))+(1/16)*gsl_matrix_get(state,i+1,j-1)));
		    }
		  
		  else if ((i==n-1)&(j==0))
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i-1,j)+gsl_matrix_get(state,i,j+1))+(1/16)*gsl_matrix_get(state,i-1,j+1)));
		    }
		  
		  else if((i==n-1)&(j==m-1))
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i-1,j)+gsl_matrix_get(state,i,j-1))+(1/16)*gsl_matrix_get(state,i-1,j-1)));
		    }
		  
		  
		  else if (i==0)
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i,j-1)+gsl_matrix_get(state,i,j+1))+(1/16)*(gsl_matrix_get(state,i+1,j-1)+gsl_matrix_get(state,i+1,j+1))));
		    }
		  
		  else if(j==0)
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i-1,j)+gsl_matrix_get(state,i,j+1))+(1/16)*(gsl_matrix_get(state,i-1,j+1)+gsl_matrix_get(state,i+1,j+1))));
		    }
		  
		  else if(i==n-1)
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i-1,j)+gsl_matrix_get(state,i,j-1)+gsl_matrix_get(state,i,j+1))+(1/16)*(gsl_matrix_get(state,i-1,j-1)+gsl_matrix_get(state,i-1,j+1))));
		    }
		  
		  else if(j==m-1)
		    { 
		  gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i,j-1)+gsl_matrix_get(state,i-1,j))+(1/16)*(gsl_matrix_get(state,i+1,j-1)+gsl_matrix_get(state,i-1,j-1))));
		    }
		  
		  else
		    { 
		      gsl_matrix_set(newstate,i,j,(0.25*gsl_matrix_get(state,i,j)+(1/8)*(gsl_matrix_get(state,i+1,j)+gsl_matrix_get(state,i,j-1)+gsl_matrix_get(state,i-1,j)+gsl_matrix_get(state,i,j+1))+(1/16)*(gsl_matrix_get(state,i+1,j-1)+gsl_matrix_get(state,i-1,j-1)+gsl_matrix_get(state,i+1,j+1)+gsl_matrix_get(state,i-1,j+1))));
		    }
		}
	    }
	  gsl_matrix_scale (newstate,a); 
	  gsl_matrix_memcpy(state,newstate); 
	  gsl_matrix_memcpy(observations,state); 
	  for(j=0;j<m;j++)
	    {
	      for (i=0;i<n;i++)
		{
		  gsl_matrix_set(observations,i,j,gsl_matrix_get(observations,i,j)+gsl_ran_gaussian(r,sqrt(1/obs_p)));
		}
	    }
	}
      /*print observations to a file*/   
      gsl_matrix_fprintf (file, observations, "%f"); 
    }     
  fclose(file);
  gsl_matrix_free(state);  
  gsl_matrix_free(newstate);  
  gsl_matrix_free(observations);
  return (0);
}
/* eof */
  
