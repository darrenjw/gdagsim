/*

gdag_meschach.h

Headers for additional Meschach functions

*/

VEC * spCHforward(SPMAT *,VEC *,VEC *);

VEC * spCHbackward(SPMAT *,VEC *,VEC *);

VEC * sp_lv_mlt(SPMAT *,VEC *,VEC *);

VEC * sp_vl_mlt(SPMAT *,VEC *,VEC *);

VEC * gdag_gsltomes(gsl_vector *);

/* eof */

