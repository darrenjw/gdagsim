/*
GDAGsim 0.3
(C) 2000-2002, Darren J Wilkinson

gdag_sparse_test.c
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
  int i,error;
  gdag_usv * usv,* usvv;
  gsl_vector * y;
  SPMAT * m;

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

  /* set up sparse matrix */
  m=sp_get(3,6,2);
  for (i=0;i<3;i++) {
    sp_set_val(m,i,i,i+5);
  }

  /* another sparse vector */
  usvv=gdag_usv_alloc(2);

  gdag_usv_add(usvv,0,5.0);
  gdag_usv_add(usvv,2,6.0);

  /* test matrix update operation */
  gdag_dusger(1.0,usvv,usv,m);

  error+=dtest("m00",5.0,sp_get_val(m,0,0),0.000001);
  error+=dtest("m01",55.0,sp_get_val(m,0,1),0.000001);
  error+=dtest("m02",0.0,sp_get_val(m,0,2),0.000001);
  error+=dtest("m03",65.0,sp_get_val(m,0,3),0.000001);
  error+=dtest("m04",70.0,sp_get_val(m,0,4),0.000001);
  error+=dtest("m05",0.0,sp_get_val(m,0,5),0.000001);

  error+=dtest("m10",0.0,sp_get_val(m,1,0),0.000001);
  error+=dtest("m11",6.0,sp_get_val(m,1,1),0.000001);
  error+=dtest("m12",0.0,sp_get_val(m,1,2),0.000001);
  error+=dtest("m13",0.0,sp_get_val(m,1,3),0.000001);
  error+=dtest("m14",0.0,sp_get_val(m,1,4),0.000001);
  error+=dtest("m15",0.0,sp_get_val(m,1,5),0.000001);

  error+=dtest("m20",0.0,sp_get_val(m,2,0),0.000001);
  error+=dtest("m21",66.0,sp_get_val(m,2,1),0.000001);
  error+=dtest("m22",7.0,sp_get_val(m,2,2),0.000001);
  error+=dtest("m23",78.0,sp_get_val(m,2,3),0.000001);
  error+=dtest("m24",84.0,sp_get_val(m,2,4),0.000001);
  error+=dtest("m25",0.0,sp_get_val(m,2,5),0.000001);

  /* test row update */

  gdag_row_update(2.0,usv,m,1);

  error+=dtest("mu00",5.0,sp_get_val(m,0,0),0.000001);
  error+=dtest("mu01",55.0,sp_get_val(m,0,1),0.000001);
  error+=dtest("mu02",0.0,sp_get_val(m,0,2),0.000001);
  error+=dtest("mu03",65.0,sp_get_val(m,0,3),0.000001);
  error+=dtest("mu04",70.0,sp_get_val(m,0,4),0.000001);
  error+=dtest("mu05",0.0,sp_get_val(m,0,5),0.000001);

  error+=dtest("mu10",0.0,sp_get_val(m,1,0),0.000001);
  error+=dtest("mu11",28.0,sp_get_val(m,1,1),0.000001);
  error+=dtest("mu12",0.0,sp_get_val(m,1,2),0.000001);
  error+=dtest("mu13",26.0,sp_get_val(m,1,3),0.000001);
  error+=dtest("mu14",28.0,sp_get_val(m,1,4),0.000001);
  error+=dtest("mu15",0.0,sp_get_val(m,1,5),0.000001);

  error+=dtest("mu20",0.0,sp_get_val(m,2,0),0.000001);
  error+=dtest("mu21",66.0,sp_get_val(m,2,1),0.000001);
  error+=dtest("mu22",7.0,sp_get_val(m,2,2),0.000001);
  error+=dtest("mu23",78.0,sp_get_val(m,2,3),0.000001);
  error+=dtest("mu24",84.0,sp_get_val(m,2,4),0.000001);
  error+=dtest("mu25",0.0,sp_get_val(m,2,5),0.000001);


  /* test column update */

  gdag_col_update(3.0,usvv,m,2);

  error+=dtest("muc00",5.0,sp_get_val(m,0,0),0.000001);
  error+=dtest("muc01",55.0,sp_get_val(m,0,1),0.000001);
  error+=dtest("muc02",15.0,sp_get_val(m,0,2),0.000001);
  error+=dtest("muc03",65.0,sp_get_val(m,0,3),0.000001);
  error+=dtest("muc04",70.0,sp_get_val(m,0,4),0.000001);
  error+=dtest("muc05",0.0,sp_get_val(m,0,5),0.000001);

  error+=dtest("muc10",0.0,sp_get_val(m,1,0),0.000001);
  error+=dtest("muc11",28.0,sp_get_val(m,1,1),0.000001);
  error+=dtest("muc12",0.0,sp_get_val(m,1,2),0.000001);
  error+=dtest("muc13",26.0,sp_get_val(m,1,3),0.000001);
  error+=dtest("muc14",28.0,sp_get_val(m,1,4),0.000001);
  error+=dtest("muc15",0.0,sp_get_val(m,1,5),0.000001);

  error+=dtest("muc20",0.0,sp_get_val(m,2,0),0.000001);
  error+=dtest("muc21",66.0,sp_get_val(m,2,1),0.000001);
  error+=dtest("muc22",25.0,sp_get_val(m,2,2),0.000001);
  error+=dtest("muc23",78.0,sp_get_val(m,2,3),0.000001);
  error+=dtest("muc24",84.0,sp_get_val(m,2,4),0.000001);
  error+=dtest("muc25",0.0,sp_get_val(m,2,5),0.000001);

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

