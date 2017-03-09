/**
 * Test flat index calculation
 * $Id: tst_nccf_get_flat_index.c 589 2011-03-23 22:01:39Z dkindig $
 */

#include <assert.h>
#include "nccf_utility_functions.h"

int main() {
  const int dims[] = {10, 11, 12};
  const int idx[] = {9, 8, 7};
  int idx2[3];
  int ndims;
  int fi;

  /* 1d */
  ndims = 1;
  fi = nccf_get_flat_index(ndims, dims, idx);
  nccf_get_multi_index(ndims, dims, fi, idx2);
  assert(fi == idx[0]);
  assert(idx[0] == idx2[0]);

  /* 2d */
  ndims = 2;
  fi = nccf_get_flat_index(ndims, dims, idx);
  nccf_get_multi_index(ndims, dims, fi, idx2);
  assert(fi == idx[1] + dims[1]*idx[0]);
  assert(idx[0] == idx2[0]);
  assert(idx[1] == idx2[1]);

  /* 3d */
  ndims = 3;
  fi = nccf_get_flat_index(ndims, dims, idx);
  nccf_get_multi_index(ndims, dims, fi, idx2);
  assert(fi == idx[2] + dims[2]*(idx[1] + dims[1]*idx[0]));
  assert(idx[0] == idx2[0]);
  assert(idx[1] == idx2[1]);
  assert(idx[2] == idx2[2]);
  
  return 0;
}
