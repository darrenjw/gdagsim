/* now interested in generating the structure for independant spatial obs on a lattice (i.e. space is 2-D).  Data inputed from prior.dat.  We have n*m spatial dim and t time! 
Using an index function which is index(i,j,t), i,j=location in space, t is location in time, and the function returns the associated point on the lattice */

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>
#include <gsl_blas.h>
#include <math.h>
#include <gdag.h>

void add_latent(gdag *, int,int ,int,double,double);
void add_obs(gdag *,gsl_vector *, double);
int indx(int,int,int);
int n;
int m; 
int T;

int main(int argc, char *argv[])
{
  gdag * d;
  gsl_rng *r;
  gsl_vector * s ,* temp2, * scratch, * data;
  gsl_vector temp1; 

  FILE * file, *file2;
  int p,i,j,k,points;
  double obs_p,state_p,sumsq,sum,a_post;
  float b,e,x;
  long it,iters;
  
  if (argc != 5) 
    {
      printf("Usage: %s <n> <m> <T> <iters> \n",argv[0]);
      exit(1);
    }
  n=atoi(argv[1]);
  m=atoi(argv[2]);
  T=atoi(argv[3]);
  iters=atoi(argv[4]);

  r=gsl_rng_alloc(gsl_rng_mt19937);
  gdag_rng_init();
  p=n*m*T;
  points=n*m;
  state_p=1.0;      /* sampled precisions */  
  obs_p=1.0;
  a_post=1.0;
  b=0.6;
  e=0.1;

  /*input data from prior.dat*/
  data=gsl_vector_calloc(n*m*T);
  if ((file=fopen("prior.dat","r")) == NULL) 
    {
      printf("Cannot open prior.dat\n");
      exit(1);
    }	
  gsl_vector_fscanf(file,data);
  fclose(file);

  /* MCMC Loop... */
  if ((file2=fopen("spt-output.tab","w")) == NULL) 
    {
      printf("Cannot open spt-output.tab\n");
      exit(1);
    }	
  fprintf(file2,"Iter state_p obs_p alpha\n");
  
  scratch=gsl_vector_alloc(p);
  d=gdag_alloc(p);
  temp2=gsl_vector_alloc(p-n*m);

  for (it=1;it<=iters;it++) 
    {
      fprintf(stderr,"%ld ",it);
      gdag_set_zero(d);
      gsl_vector_set_zero(temp2);

      add_latent(d,m,n,T,a_post,state_p);
      add_obs(d,data,obs_p); 
      gdag_process(d);
      s=gdag_sim(d);

      gsl_vector_memcpy(scratch,s);      /* estimate obs_p */
      gsl_vector_sub(scratch,data);
      gsl_blas_ddot(scratch,scratch,&sumsq);
      obs_p=gdag_ran_gamma(0.1+0.5*p,0.1+0.5*sumsq); 
      
      gsl_vector_memcpy(scratch,s);      /* estimate alpha */
      temp1=gsl_vector_subvector(scratch,n*m,p-n*m);

      for(k=0;k<T-1;k++)                 /* ignore time 0*/
	{
	  for(j=0;j<m;j++)
	    {
	      for(i=0;i<n;i++ )
		{
		  if ((i==0)&(j==0))
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*gsl_vector_get(scratch,indx(i+1,j+1,k))));
		    }
		  else if((i==0)&(j==m-1))
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k)))+(1/16)*gsl_vector_get(scratch,indx(i+1,j-1,k))));
		    }
		  
		  else if ((i==n-1)&(j==0))
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i-1,j,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*gsl_vector_get(scratch,indx(i-1,j+1,k))));
		    }
		  
		  else if ((i==n-1)&(j==m-1))
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i-1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k)))+(1/16)*gsl_vector_get(scratch,indx(i-1,j-1,k))));
		    }
		  /*edges*/
		  else if(j==0)
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i-1,j,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*(gsl_vector_get(scratch,indx(i-1,j+1,k))+gsl_vector_get(scratch,indx(i+1,j+1,k)))));
		    }
		  
		  else if(i==0)
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*(gsl_vector_get(scratch,indx(i+1,j-1,k))+gsl_vector_get(scratch,indx(i+1,j+1,k)))));
		    }
		  
		  else if(i==(n-1))
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i-1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*(gsl_vector_get(scratch,indx(i-1,j-1,k))+gsl_vector_get(scratch,indx(i-1,j+1,k)))));
		    }
		  
		  else if (j==m-1)
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k))+gsl_vector_get(scratch,indx(i-1,j,k)))+(1/16)*(gsl_vector_get(scratch,indx(i+1,j-1,k))+gsl_vector_get(scratch,indx(i-1,j-1,k)))));
		    }
		  /*middle*/		  
		  else
		    { 
		      gsl_vector_set(temp2,indx(i,j,k),(0.25*gsl_vector_get(scratch,indx(i,j,k))+(1/8)*(gsl_vector_get(scratch,indx(i+1,j,k))+gsl_vector_get(scratch,indx(i,j-1,k))+gsl_vector_get(scratch,indx(i-1,j,k))+gsl_vector_get(scratch,indx(i,j+1,k)))+(1/16)*(gsl_vector_get(scratch,indx(i+1,j-1,k))+gsl_vector_get(scratch,indx(i-1,j-1,k))+gsl_vector_get(scratch,indx(i+1,j+1,k))+gsl_vector_get(scratch,indx(i-1,j+1,k)))));
		    }
		}
	    }
	}
      gsl_blas_ddot(temp2,temp2,&sumsq);
      gsl_blas_ddot(&temp1,temp2,&sum);
      x=gsl_ran_gaussian(r,sqrt(1/(state_p*sumsq+1/e)));
      a_post=(e*sum+b/state_p)/(e*sumsq+1/state_p)+x;
      
      gsl_vector_scale (temp2,a_post);         /* estimate state_p */
      gsl_vector_sub(&temp1,temp2); 
      gsl_blas_ddot(&temp1,&temp1,&sumsq); 
      state_p=gdag_ran_gamma(0.1+0.5*(p-n*m),0.1+0.5*sumsq);
      
      fprintf(file2,"%ld %f %f %f\n",it,log(1/state_p),log(1/obs_p),a_post);
    }
   
  gsl_vector_free(temp2);
  gdag_free(d);
  fclose(file2);
  return(0);
}
  
int indx(int i, int j,int t)
{
  int indx;
  indx=(t*n*m)+(j*n)+i;
  return(indx);
}


void add_latent(gdag * d,int m, int n, int T,double a,double prec)
{
  int i,j,k;
  gdag_usv * alpha;
  for(k=0;k<T;k++)
    {
      for(j=0;j<m;j++)
	{
	  for(i=0;i<n;i++)
	    {
	      if(k==0)  /*time 0 */
		{
		  gdag_add_root(d,indx(i,j,k),0.0,prec);
		}
	      /*corners*/
	      else if ((i==0)&(j==0))
		{ 
		  alpha=gdag_usv_alloc(4);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	      else if((i==0)&(j==m-1))
		{ 
		  alpha=gdag_usv_alloc(4);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j-1,k-1),a/16);		  
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	      
	      else if ((i==n-1)&(j==0))
		{ 
		  alpha=gdag_usv_alloc(4);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);		  
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	      
	      else if ((i==n-1)&(j==m-1))
		
		{ 
		  alpha=gdag_usv_alloc(4);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j-1,k-1),a/16);		  
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	      /*edges*/
	      else if(j==0)
		{ 
		  alpha=gdag_usv_alloc(6);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i-1,j+1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	  
	      else if(i==0)
		{ 
		  alpha=gdag_usv_alloc(6);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i+1,j-1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}

	      else if(i==(n-1))
		{ 
		  alpha=gdag_usv_alloc(6);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i-1,j-1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}

	  
	      else if (j==m-1)
		{ 
		  alpha=gdag_usv_alloc(6);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i-1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i-1,j+1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	     
	      else
		{ 
		  alpha=gdag_usv_alloc(9);
		  gdag_usv_add(alpha,indx(i,j,k-1),a/4);
		  gdag_usv_add(alpha,indx(i,j+1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j,k-1),a/8);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i,j-1,k-1),a/8);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i+1,j+1,k-1),a/16);
		  gdag_usv_add(alpha,indx(i+1,j-1,k-1),a/16);
		  gdag_add_node(d,indx(i,j,k),alpha,0,prec);
		  gdag_usv_free(alpha);
		}
	    }
	}
    }
}

void add_obs(gdag * d,gsl_vector * data, double prec)
{
  gdag_usv * alpha;
  int i;
  for (i=0;i < (data->size);i++) 
	{
	  alpha=gdag_usv_basis(i);
	  gdag_add_observation(d,alpha,0.0,prec,gsl_vector_get(data,i));
	  gdag_usv_free(alpha);
	}
}

/* eof */
