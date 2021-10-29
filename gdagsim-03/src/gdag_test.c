/*
GDAGsim 0.3
(C) 2000-2002, Darren J Wilkinson

gdag_test.c
Test code for the GDAGsim library
*/

#include <math.h>

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>

#include <gdag.h>

void _FPU_SETCW(int * cw){}  /* crap for compaq alphas only */

int dtest(char *,double,double,double);

int main(void)
{
  gdag *d;
  gdag_usv *alpha;
  gsl_vector *s, *x;
  gsl_rng *r;
  int error;
  long i,n;
  double sum0,sum1,sum2,sum3,sum4,
    ss0,ss1,ss2,ss3,ss4,
    s0,s1,s2,s3,s4;
  
  error=0;
  n=1000000;


  /* test prior structure */
  d=gdag_calloc(3);
  r=gsl_rng_alloc(gsl_rng_mt19937);

  gdag_add_root(d,0,1.0,2.0);
  gdag_add_root(d,1,3.0,4.0); 
  alpha=gdag_usv_alloc(2);
  gdag_usv_add(alpha,0,1.0);
  gdag_usv_add(alpha,1,2.0);
  gdag_add_node(d,2,alpha,0,1);
  gdag_usv_free(alpha);
  gdag_process(d);

  alpha=gdag_usv_basis(0);
  error+=dtest("v0",0.5,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);
  alpha=gdag_usv_basis(1);
  error+=dtest("v1",0.25,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);
  alpha=gdag_usv_basis(2);
  error+=dtest("v2",2.5,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);

  s=gdag_mean(d);
  error+=dtest("e1",1.0,gsl_vector_get(s,0),0.00001);
  error+=dtest("e2",3.0,gsl_vector_get(s,1),0.00001);
  error+=dtest("e3",7.0,gsl_vector_get(s,2),0.00001);

  sum0=0;sum1=0;sum2=0;
  ss0=0;ss1=0;ss2=0;
  for (i=0;i<n;i++) {
    s=gdag_sim(r,d);
    s0=gsl_vector_get(s,0);
    s1=gsl_vector_get(s,1);
    s2=gsl_vector_get(s,2);
    sum0+=s0;
    sum1+=s1;
    sum2+=s2;
    ss0+=(s0-1)*(s0-1);
    ss1+=(s1-3)*(s1-3);
    ss2+=(s2-7)*(s2-7);
  }
  error+=dtest("sm0",1.0,sum0/n,0.005);
  error+=dtest("sm1",3.0,sum1/n,0.005);
  error+=dtest("sm2",7.0,sum2/n,0.005);
  error+=dtest("sv0",0.5,ss0/n,0.01);
  error+=dtest("sv1",0.25,ss1/n,0.01);
  error+=dtest("sv2",2.5,ss2/n,0.01); 

  gdag_free(d);


  /* test adjusted structure */
  d=gdag_calloc(3);

  gdag_add_root(d,0,1.0,2.0);
  gdag_add_root(d,1,3.0,4.0); 
  alpha=gdag_usv_alloc(2);
  gdag_usv_add(alpha,0,1.0);
  gdag_usv_add(alpha,1,2.0);
  gdag_add_node(d,2,alpha,0,1);
  gdag_usv_free(alpha);

  gdag_prior_process(d);

  alpha=gdag_usv_alloc(2);
  gdag_usv_add(alpha,1,1.0);
  gdag_usv_add(alpha,2,1.0);
  gdag_add_observation(d,alpha,5,0.5,20);
  gdag_usv_free(alpha);

  gdag_process(d);

  error+=dtest("mll",log(gsl_ran_gaussian_pdf(20-15,sqrt(5.75))),
	       gdag_mloglik(d),0.00001);
  s=gdag_sim(r,d);
  x=gsl_vector_calloc(3);
  gsl_vector_memcpy(x,s);
  /* cross-check vll against ll, but should really have a proper test! */
  error+=dtest("vll",gdag_loglik(d),gdag_vloglik(d,x),0.00001);

  s=gdag_mean(d);
  error+=dtest("ae1",33.0/23,gsl_vector_get(s,0),0.00001);
  error+=dtest("ae2",84.0/23,gsl_vector_get(s,1),0.00001);
  error+=dtest("ae3",221.0/23,gsl_vector_get(s,2),0.00001);
  
  alpha=gdag_usv_basis(0);
  error+=dtest("av0",42.0/92,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);
  alpha=gdag_usv_basis(1);
  error+=dtest("av1",14.0/92,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);
  alpha=gdag_usv_basis(2);
  error+=dtest("av2",86.0/92,gdag_var(d,alpha),0.00001);
  gdag_usv_free(alpha);

  sum0=0;sum1=0;sum2=0;
  ss0=0;ss1=0;ss2=0;
  for (i=0;i<n;i++) {
    s=gdag_sim(r,d);
    s0=gsl_vector_get(s,0);
    s1=gsl_vector_get(s,1);
    s2=gsl_vector_get(s,2);
    sum0+=s0;
    sum1+=s1;
    sum2+=s2;
    ss0+=(s0-33.0/23)*(s0-33.0/23);
    ss1+=(s1-84.0/23)*(s1-84.0/23);
    ss2+=(s2-221.0/23)*(s2-221.0/23);
  }
  error+=dtest("sam0",33.0/23,sum0/n,0.005);
  error+=dtest("sam1",84.0/23,sum1/n,0.005);
  error+=dtest("sam2",221.0/23,sum2/n,0.005);
  error+=dtest("sav0",42.0/92,ss0/n,0.01);
  error+=dtest("sav1",14.0/92,ss1/n,0.01);
  error+=dtest("sav2",86.0/92,ss2/n,0.01); 

  gdag_free(d);


  /* misc other tests */
  sum0=0;sum1=0;sum2=0;sum3=0;sum4=0;ss2=0;ss3=0;ss4=0;
  for (i=0;i<n;i++) {
    if (gdag_accept_p(r,0.25)) { sum0++; }
    if (gdag_accept_lp(r,log(0.4))) { sum1++; }
    s2=gdag_ran_gamma(r,2.0,3.0);
    sum2+=s2;
    ss2+=gdag_sqr(s2-(2.0/3.0));
    s3=gdag_ran_gaussian(r,3.0,0.2);
    sum3+=s3;
    ss3+=gdag_sqr(s3-3.0);
    s4=gdag_ran_uniform(r,1.0,3.0);
    sum4+=s4;
    ss4+=gdag_sqr(s4-2.0);
  }
  error+=dtest("ap",0.25,sum0/n,0.005);
  error+=dtest("alp",0.4,sum1/n,0.005);
  error+=dtest("egam",(2.0/3.0),sum2/n,0.005);
  error+=dtest("vgam",(2.0/9.0),ss2/n,0.01);
  error+=dtest("egauss",3.0,sum3/n,0.005);
  error+=dtest("vgauss",5.0,ss3/n,0.01);
  error+=dtest("eunif",2.0,sum4/n,0.005);
  error+=dtest("vunif",1.0/3.0,ss4/n,0.01);




  /* end of testing summary */
  if (error>0) {
    printf("***** %d TEST(S) FAILED. *****\n",error);
    return(1);
  }
  printf("All tests passed.\n");
  return(0);

}


int dtest(char * name,double exp,double act,double tol)
{
  if ((exp-act)*(exp-act) > tol*tol) {
    printf("ERROR TESTING %s: Expected: %f, Actual: %f, Tolerance: %f\n",
	   name,exp,act,tol);
    return(1);
  }
  printf("%s test passed.\n",name);
  return(0);
}


/* eof */

