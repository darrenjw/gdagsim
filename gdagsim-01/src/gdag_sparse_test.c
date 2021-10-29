/*
GDAGsim 0.1
(C) 2000, Darren J Wilkinson

gdag_test.c
Test code for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_linalg.h>

#include <gdag.h>

int dtest(char *,double,double,double);

int main(void)
{
  int i,j,error;
  gdag_usv * usv,* usvv;
  gsl_vector * y;
  gsl_matrix * m;

  error=0;

  /* set up dense vector */
  y=gsl_vector_alloc(6);
  for (i=0;i<6;i++) {
    gsl_vector_set(y,i,i+1);
  }

  /* set up sparse vector */
  usv=gdag_usv_alloc(3);
  gdag_usv_add(usv,1,11.0);
  gdag_usv_add(usv,3,13.0);
  gdag_usv_add(usv,4,14.0);

  /* test sparse update operation */
  gdag_dusaxpy(2.0,usv,y);

  error+=dtest("y0",1.0,gsl_vector_get(y,0),0.000001);
  error+=dtest("y1",24.0,gsl_vector_get(y,1),0.000001);
  error+=dtest("y2",3.0,gsl_vector_get(y,2),0.000001);
  error+=dtest("y3",30.0,gsl_vector_get(y,3),0.000001);
  error+=dtest("y4",33.0,gsl_vector_get(y,4),0.000001);
  error+=dtest("y5",6.0,gsl_vector_get(y,5),0.000001);

  /* set up dense matrix */
  m=gsl_matrix_alloc(3,6);
  for (i=0;i<3;i++) {
    for (j=0;j<6;j++) {
      gsl_matrix_set(m,i,j,i+j);
    }
  }

  /* another sparse vector */
  usvv=gdag_usv_alloc(2);
  gdag_usv_add(usvv,0,5.0);
  gdag_usv_add(usvv,2,6.0);

  /* test matrix update operation */
  gdag_dusger(1.0,usvv,usv,m);

  error+=dtest("m00",0.0,gsl_matrix_get(m,0,0),0.000001);
  error+=dtest("m01",56.0,gsl_matrix_get(m,0,1),0.000001);
  error+=dtest("m02",2.0,gsl_matrix_get(m,0,2),0.000001);
  error+=dtest("m03",68.0,gsl_matrix_get(m,0,3),0.000001);
  error+=dtest("m04",74.0,gsl_matrix_get(m,0,4),0.000001);
  error+=dtest("m05",5.0,gsl_matrix_get(m,0,5),0.000001);

  error+=dtest("m10",1.0,gsl_matrix_get(m,1,0),0.000001);
  error+=dtest("m11",2.0,gsl_matrix_get(m,1,1),0.000001);
  error+=dtest("m12",3.0,gsl_matrix_get(m,1,2),0.000001);
  error+=dtest("m13",4.0,gsl_matrix_get(m,1,3),0.000001);
  error+=dtest("m14",5.0,gsl_matrix_get(m,1,4),0.000001);
  error+=dtest("m15",6.0,gsl_matrix_get(m,1,5),0.000001);

  error+=dtest("m20",2.0,gsl_matrix_get(m,2,0),0.000001);
  error+=dtest("m21",69.0,gsl_matrix_get(m,2,1),0.000001);
  error+=dtest("m22",4.0,gsl_matrix_get(m,2,2),0.000001);
  error+=dtest("m23",83.0,gsl_matrix_get(m,2,3),0.000001);
  error+=dtest("m24",90.0,gsl_matrix_get(m,2,4),0.000001);
  error+=dtest("m25",7.0,gsl_matrix_get(m,2,5),0.000001);

  /* test sparse scatter */
  gdag_dussc(usvv,y);
  
  error+=dtest("yy0",5.0,gsl_vector_get(y,0),0.000001);
  error+=dtest("yy1",24.0,gsl_vector_get(y,1),0.000001);
  error+=dtest("yy2",6.0,gsl_vector_get(y,2),0.000001);
  error+=dtest("yy3",30.0,gsl_vector_get(y,3),0.000001);
  error+=dtest("yy4",33.0,gsl_vector_get(y,4),0.000001);
  error+=dtest("yy5",6.0,gsl_vector_get(y,5),0.000001);

  gdag_usv_free(usv);

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

