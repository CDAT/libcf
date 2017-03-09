/**
 * Create a longitude-latitude grid.
 *
 * "$Id: ex_create_lonlat.c 784 2011-07-14 19:53:33Z pletzer $"
 */

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <netcdf.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"
#include "nccf_mosaic.h"
#include "nccf_host.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

#include "examples.h"

int ex_create_lonlat(const int dims[], 
                     const double coord_mins[], const double coord_maxs[], 
                     int coordids[], int *grididp){
  const int save = 1;
  const int ndims = 2;
  int i, j, ntot, status;
  double dLon, dLat;
  double *lonData;
  double *latData;
  const char *dimnames[] = {"nj", "ni"};
  char scrip_filename[STRING_SIZE];
  dLon = (coord_maxs[0] - coord_mins[0]) / (dims[0] - 1);
  dLat = (coord_maxs[1] - coord_mins[1]) / (dims[1] - 1);
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }
  lonData = malloc(ntot * sizeof(double));
  latData = malloc(ntot * sizeof(double));
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      lonData[i + dims[0] * j] = coord_mins[0] + i * dLon;
      latData[i + dims[0] * j] = coord_mins[1] + j * dLat;
    }
  }
  if ((status = nccf_def_lon_coord(ndims, dims, dimnames, lonData, save, 
                                   &coordids[0]))) ERR;
  if ((status = nccf_def_lat_coord(ndims, dims, dimnames, latData, save, 
                                   &coordids[1]))) ERR;
  free(lonData);
  free(latData);

  /* Create grid */
  if ((status = nccf_def_grid(coordids, "lonlat", grididp))) ERR;
  strcpy(scrip_filename, "gsex2_lonlat_scrip.nc");
  if ((status = nccf_save_grid_scrip(*grididp, scrip_filename))) ERR;

  return 0;
}
