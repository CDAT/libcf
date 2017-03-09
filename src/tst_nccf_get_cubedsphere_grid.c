/**
 * $Id: tst_nccf_get_cubedsphere_grid.c 560 2011-03-14 15:53:54Z pletzer $
 *
 * Test cubed sphere grid generation.
 */
// std includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <netcdf.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

int main() {
  int faceVec[] = {1, 0, 0};
  int faceVec2[] = {0, 1, 0};
  int dims[] = {3, 3};
 
  double lats[9];
  double lons[9];

  nccf_get_cubedsphere_grid(dims, faceVec, lons, lats);
  
  int i,j;
  for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
#ifdef LOGGING
      printf("%f ", lons[i + j*dims[0]]);
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
  }
  for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
#ifdef LOGGING
      printf("%f ", lats[i + j*dims[0]]);
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
  }
#ifdef LOGGING
  printf("\n");
#endif

  nccf_get_cubedsphere_grid(dims, faceVec2, lons, lats);

  for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
#ifdef LOGGING
      printf("%f ", lons[i + j*dims[0]]);
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
  }
  for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
#ifdef LOGGING
      printf("%f ", lats[i + j*dims[0]]);
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
  }

  return 0;
}
