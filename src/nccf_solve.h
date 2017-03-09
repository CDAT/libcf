/**
 * $Id: nccf_solve.h 846 2011-10-06 03:59:42Z pletzer $
 */

int 
nccf_solve_TYPE(int ndims, _TYPE mat[], const _TYPE rhs[], _TYPE sol[]) {
  int info = 1;
  
#ifdef HAVE_LAPACK_LIB
  int lda = ndims;
  int ldb = ndims;
  int ipv[ndims];
  int i;
  _GETRF(&ndims, &ndims, &mat[0], &lda,  &ipv[0], &info); 
  if (info == 0) {
    char meth = 'T'; // transpose
    int nrhs = 1;
    for (i = 0; i < ndims; ++i) {
      sol[i] = rhs[i];
    }
    _GETRS(&meth, &ndims, &nrhs, &mat[0], &lda,  &ipv[0], 
	   &sol[0], &ldb, &info);
  }
#else
  // no LAPACK
  info = 12345;
#endif

  return info;
}
