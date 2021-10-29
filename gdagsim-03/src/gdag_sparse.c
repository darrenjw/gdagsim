/*
GDAGsim 0.3
(C) 2000-2002, Darren J Wilkinson

gdag_sparse.c
Source file for the GDAGsim library
*/

#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_rng.h>
#include <gdag.h>

/* sparse vector functions */

/* allocation */

gdag_usv * gdag_usv_alloc(int nz)
{
  gdag_usv * usv;
  usv=malloc(sizeof(gdag_usv));
  usv->nz = nz;
  usv->x = gsl_vector_alloc(nz);
  usv->indx = gsl_vector_int_alloc(nz);
  usv->count = 0;
  return(usv);
}

void gdag_usv_free(gdag_usv * usv)
{
  gsl_vector_free(usv->x);
  gsl_vector_int_free(usv->indx);
  free(usv);
}

gdag_usv * gdag_usv_basis(int i)
{
  gdag_usv * usv;
  usv=gdag_usv_alloc(1);
  gsl_vector_set(usv->x,0,1);
  gsl_vector_int_set(usv->indx,0,i);
  return(usv);
}

/* build */

void gdag_usv_add(gdag_usv * usv,int indx,double x)
{
  if (usv->count >= usv->nz) {
    printf("Error: Out of bounds in gdag_usv_add.\n");
    exit(1);
  }
  gsl_vector_int_set(usv->indx,usv->count,indx);
  gsl_vector_set(usv->x,usv->count,x);
  usv->count++;
}

/* sparse BLAS style operations */

void gdag_dusaxpy(double alpha,gdag_usv * x,gsl_vector * y)
{
  int i;
  for (i=0;i < (x->nz);i++) {
    gsl_vector_set(y,
		   gsl_vector_int_get(x->indx,i),
		   gsl_vector_get(y,gsl_vector_int_get(x->indx,i)) 
		   + alpha*gsl_vector_get(x->x,i)
		   );
  }
}

void gdag_row_update(double alpha,gdag_usv * v,SPMAT * m,size_t idx)
{
  int i;
  for (i=0;i<(v->nz);i++) {
    sp_set_val(m,idx,gsl_vector_int_get(v->indx,i),
	       alpha*gsl_vector_get(v->x,i)
	       + sp_get_val(m,idx,gsl_vector_int_get(v->indx,i))
	       );
  }
}

void gdag_col_update(double alpha,gdag_usv * v,SPMAT * m,size_t idx)
{
  int i;
  for (i=0;i<(v->nz);i++) {
    sp_set_val(m,gsl_vector_int_get(v->indx,i),idx,
	       alpha*gsl_vector_get(v->x,i)
	       + sp_get_val(m,gsl_vector_int_get(v->indx,i),idx)
	       );
  }
}

void gdag_dusger(double alpha,gdag_usv * x,gdag_usv * y,SPMAT * a)
{
  int i,j;
  for (i=0;i < (x->nz);i++) {
    for (j=0;j < (y->nz);j++) {
      sp_set_val(a,
		     gsl_vector_int_get(x->indx,i),
		     gsl_vector_int_get(y->indx,j),
		     sp_get_val(a,gsl_vector_int_get(x->indx,i),
				      gsl_vector_int_get(y->indx,j))
		     + alpha*gsl_vector_get(x->x,i)
		            *gsl_vector_get(y->x,j)
		     );
    }
  }
}

void gdag_dussc(gdag_usv * x,gsl_vector * y)
{
  int i;
  for (i=0;i < (x->nz);i++) {
    gsl_vector_set(y,gsl_vector_int_get(x->indx,i),
		     gsl_vector_get(x->x,i));
  }
}

void gdag_dusdot(gdag_usv * x,gsl_vector * y,double * res)
{
  int i;
  *res=0.0;
  for (i=0;i < (x->nz);i++) {
    *res+=gsl_vector_get(x->x,i)
              *gsl_vector_get(y,gsl_vector_int_get(x->indx,i));
  }
}

/* accessor */

void gdag_usv_dump(gdag_usv * usv)
{
  int i;
  printf("nz: %d\n",usv->nz);
  printf("x: (%f",gsl_vector_get(usv->x,0));
  for (i=1;i < usv->nz;i++) {
    printf(",%f",gsl_vector_get(usv->x,i));
  }
  printf(")\n");
  printf("indx: (%d",gsl_vector_int_get(usv->indx,0));
  for (i=1;i < usv->nz;i++) {
    printf(",%d",gsl_vector_int_get(usv->indx,i));
  }
  printf(")\n");
  printf("count: %d\n",usv->count);
}



/* eof */

