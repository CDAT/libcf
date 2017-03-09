/**
 * $Id: tst_solve.c 828 2011-09-14 20:05:08Z pletzer $
 */

#ifdef HAVE_CONFIG_H
#include <cf_config.h>
#endif

#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>
#include <math.h>
#include <assert.h>

int main() {

#ifdef HAVE_LAPACK_LIB
  const int n = 3;
  double mat[n*n];
  double matCopy[n*n];
  double rhs[n];
  double sol[n];
  int i, j, status;
  for (i = 0; i < n; ++i) {
    rhs[i] = 2.0*i + 1.0;
    for (j = 0; j < n; ++j) {
      mat[j + n*i] = i - j;
      if (i == j) mat[j + n*i] = 10.0;
      matCopy[j + n*i] = mat[j + n*i];
    }
  }
  status = nccf_solve_double(n, mat, rhs, sol);
  if (status) ERR;

  // check
  double totError = 0.0;
  for (i = 0; i < n; ++i) {
    double sum = - rhs[i];
    for (j = 0; j < n; ++j) {
      sum += matCopy[i*n + j] * sol[j];
    }
    totError += fabs(sum);
  }
  assert(totError < 1.e-10);
#endif

  return 0;
}
