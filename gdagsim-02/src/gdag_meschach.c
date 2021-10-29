/* 

gdag_meschach.c

Just copies of the meschach function spCHsolve with the
   forward (backwards) bits deleted (respectively). 

Also modified copies of the multn functions for lower
   triangular matrices.

Finally, a function creating a Meschach vector view of a GSL vector.

*/

#include "gsl_matrix.h"
#include "sparse2.h"

/* spCHforward -- solve L.out=b where L is a sparse matrix,
	-- out, b dense vectors
	-- returns out; operation may be in-situ */

VEC	*spCHforward(L,b,out)
SPMAT	*L;
VEC	*b, *out;
{
	int	i, j_idx, n;
	SPROW	*row;
	row_elt	*elt;
	Real	sum, *out_ve;

	if ( L == SMNULL || b == VNULL )
		error(E_NULL,"spCHforward");
	if ( L->m != L->n )
		error(E_SQUARE,"spCHforward");
	if ( b->dim != L->m )
		error(E_SIZES,"spCHforward");

	if ( ! L->flag_col )
		sp_col_access(L);
	if ( ! L->flag_diag )
		sp_diag_access(L);

	out = v_copy(b,out);
	out_ve = out->ve;

	/* forward substitution: solve L.x=b for x */
	n = L->n;
	for ( i = 0; i < n; i++ )
	{
		sum = out_ve[i];
		row = &(L->row[i]);
		elt = row->elt;
		for ( j_idx = 0; j_idx < row->len; j_idx++, elt++ )
		{
		    if ( elt->col >= i )
			break;
		    sum -= elt->val*out_ve[elt->col];
		}
		if ( row->diag >= 0 )
		    out_ve[i] = sum/(row->elt[row->diag].val);
		else
		    error(E_SING,"spCHforward");
	}

	return out;
}


/* spCHbackward -- solve L^T.out=b where L is a sparse matrix,
	-- out, b dense vectors
	-- returns out; operation may be in-situ */
VEC	*spCHbackward(L,b,out)
SPMAT	*L;
VEC	*b, *out;
{
	int	i, n, scan_idx, scan_row;
	SPROW	*row;
	row_elt	*elt;
	Real	diag_val, sum, *out_ve;

	if ( L == SMNULL || b == VNULL )
		error(E_NULL,"spCHbackward");
	if ( L->m != L->n )
		error(E_SQUARE,"spCHbackward");
	if ( b->dim != L->m )
		error(E_SIZES,"spCHbackward");

	if ( ! L->flag_col )
		sp_col_access(L);
	if ( ! L->flag_diag )
		sp_diag_access(L);

	out = v_copy(b,out);
	out_ve = out->ve;

	n = L->n;

	/* backward substitution: solve L^T.out = x for out */
	for ( i = n-1; i >= 0; i-- )
	{
		sum = out_ve[i];
		row = &(L->row[i]);
		/* Note that row->diag >= 0 by above loop */
		elt = &(row->elt[row->diag]);
		diag_val = elt->val;

		/* scan down column */
		scan_idx = elt->nxt_idx;
		scan_row = elt->nxt_row;
		while ( scan_row >= 0 /* && scan_idx >= 0 */ )
		{
		    row = &(L->row[scan_row]);
		    elt = &(row->elt[scan_idx]);
		    sum -= elt->val*out_ve[scan_row];
		    scan_idx = elt->nxt_idx;
		    scan_row = elt->nxt_row;
		}
		out_ve[i] = sum/diag_val;
	}

	return out;
}


/* sp_lv_mlt -- sparse matrix/dense vector multiply
   -- result is in out, which is returned unless out==NULL on entry
   --  if out==NULL on entry then the result vector is created */
VEC	*sp_lv_mlt(A,x,out)
SPMAT	*A;
VEC	*x, *out;
{
   int	i, j_idx, m, n, max_idx;
   Real	sum, *x_ve;
   SPROW	*r;
   row_elt	*elts;
   
   if ( ! A || ! x )
     error(E_NULL,"sp_lv_mlt");
   if ( x->dim != A->n )
     error(E_SIZES,"sp_lv_mlt");
   if ( ! out || out->dim < A->m )
     out = v_resize(out,A->m);
   if ( out == x )
     error(E_INSITU,"sp_lv_mlt");
   m = A->m;	n = A->n;
   x_ve = x->ve;
   
   for ( i = 0; i < m; i++ )
   {
      sum = 0.0;
      r = &(A->row[i]);
      max_idx = r->len;
      elts    = r->elt;
      /*      for ( j_idx = 0; j_idx < max_idx; j_idx++, elts++ ) */
      for ( j_idx = 0; ((j_idx < max_idx) && ((elts->col) <= i)) ; j_idx++, elts++ )
	sum += elts->val*x_ve[elts->col];
      out->ve[i] = sum;
   }
   return out;
}


/* sp_vl_mlt -- sparse matrix/dense vector multiply from left
   -- result is in out, which is returned unless out==NULL on entry
   -- if out==NULL on entry then result vector is created & returned */
VEC	*sp_vl_mlt(A,x,out)
SPMAT	*A;
VEC	*x, *out;
{
   int	i, j_idx, m, n, max_idx;
   Real	tmp, *x_ve, *out_ve;
   SPROW	*r;
   row_elt	*elts;
   
   if ( ! A || ! x )
     error(E_NULL,"sp_vl_mlt");
   if ( x->dim != A->m )
     error(E_SIZES,"sp_vl_mlt");
   if ( ! out || out->dim < A->n )
     out = v_resize(out,A->n);
   if ( out == x )
     error(E_INSITU,"sp_vl_mlt");
   
   m = A->m;	n = A->n;
   v_zero(out);
   x_ve = x->ve;	out_ve = out->ve;
   
   for ( i = 0; i < m; i++ )
   {
     /* r = A->row+i; */
      r = &(A->row[i]);
      /* attempted fix */
      max_idx = r->len;
      elts    = r->elt;
      tmp = x_ve[i];
      /*      for ( j_idx = 0; j_idx < max_idx; j_idx++, elts++ ) */
      for ( j_idx = 0; ((j_idx < max_idx) && ((elts->col) <= i)) ; j_idx++, elts++ )
	out_ve[elts->col] += elts->val*tmp;
   }
   
   return out;
}


VEC * gdag_gsltomes(gsl_vector * g)
{
  VEC * m;
  m=malloc(sizeof(VEC));
  m->dim=g->size;
  m->max_dim=g->size;
  m->ve=g->data;
  return(m);
}





/* eof */

