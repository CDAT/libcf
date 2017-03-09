/**
 * $Id: tst_nccf_get_bipolar_cap.c 560 2011-03-14 15:53:54Z pletzer $
 *
 * Test bipolar cap grid generation.
 */
// std includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <netcdf.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

int main() {
  double latPerim = 66.;
  int dims[] = {5, 5};
  int nvert = dims[0] * dims[1];
 
  double lats[nvert];
  double lons[nvert];
  double lonSing = 0.0;

  nccf_get_bipolar_cap(dims, latPerim, lonSing, lons, lats);
  
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
  return 0;
}
