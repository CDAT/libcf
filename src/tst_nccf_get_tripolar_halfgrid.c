/**
 * $Id: tst_nccf_get_tripolar_halfgrid.c 560 2011-03-14 15:53:54Z pletzer $
 *
 * Test tripolar grid generation, another way to generate the tripolar grid.
 */
// std includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <netcdf.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

int main() {
  int dims[] = {5, 10};
  int capIndx = 7;

  double lats[5*10];
  double lons[5*10];

  nccf_get_tripolar_halfgrid(dims, 0, capIndx, lons, lats);
  
  int i,j;
#ifdef LOGGING
  printf("first half:\n");
#endif
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
  printf("\nsecond half:\n");
#endif

  nccf_get_tripolar_halfgrid(dims, 1, capIndx, lons, lats);
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
