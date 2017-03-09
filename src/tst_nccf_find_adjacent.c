/**
 * Test flat index calculation
 * $Id: tst_nccf_find_adjacent.c 589 2011-03-23 22:01:39Z dkindig $
 */

#include <assert.h>
#include <stdio.h>
#include <nccf_utility_functions.h>

int main() {

  int i;
  /* Did not inside acutally get tripped? */
  int inside_domain = 1, ntot = 1, count = 0;

  /* 4D */
  int dims4d[] = {2, 3, 2, 3};
  int kji1[] = {0, 0, 0, 0};
  int kji2[] = {0, 0, 0, 0};
  int ndims = 4;

  for( i = 0; i < ndims; i++ ) ntot *= dims4d[i];

  while( inside_domain ){
    for( i = 0; i < ndims; i++ ) kji2[i] = kji1[i];
    inside_domain = nccf_find_next_adjacent( ndims, dims4d, kji1 );
    ++count;

    if( count > ntot ) assert( count == ntot );
  }

  assert( count == ntot );

  /* 3D */
  int dims3d[] = {2, 3, 2};
  kji1[0] = 0; kji1[1] = 0; kji1[2] = 0;
  kji2[0] = 0; kji2[1] = 0; kji2[2] = 0;
  ndims = 3;
  inside_domain = 1; ntot = 1; count = 0;

  for( i = 0; i < ndims; i++ ) ntot *= dims3d[i];

  while( inside_domain ){
    for( i = 0; i < ndims; i++ ) kji2[i] = kji1[i];
    inside_domain = nccf_find_next_adjacent( ndims, dims3d, kji1 );
    ++count;

    if( count > ntot ) assert( count == ntot );
  }

  assert( count == ntot );

  /* 2D */
  int dims2d[] = {2, 3};
  kji1[0] = 0; kji1[1] = 0;
  kji2[0] = 0; kji2[1] = 0;
  ndims = 2;
  inside_domain = 1; ntot = 1; count = 0;

  for( i = 0; i < ndims; i++ ) ntot *= dims2d[i];

  while( inside_domain ){
    for( i = 0; i < ndims; i++ ) kji2[i] = kji1[i];
    inside_domain = nccf_find_next_adjacent( ndims, dims2d, kji1 );
    ++count;

    if( count > ntot ) assert( count == ntot );
  }

  assert( count == ntot );

  /* 1D */
  int dims1d[] = {2};
  kji1[0] = 0;
  kji2[0] = 0;
  ndims = 1;
  inside_domain = 1; ntot = 1; count = 0;

  for( i = 0; i < ndims; i++ ) ntot *= dims1d[i];

  while( inside_domain ){
    for( i = 0; i < ndims; i++ ) kji2[i] = kji1[i];
    inside_domain = nccf_find_next_adjacent( ndims, dims1d, kji1 );
    ++count;

    if( count > ntot ) assert( count == ntot );
  }

  assert( count == ntot );

  return 0;

}
