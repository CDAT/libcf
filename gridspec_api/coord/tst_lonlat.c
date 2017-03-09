/**
 * Test copy-save switch in coordinate creation
 *
 * $Id: tst_lonlat.c 719 2011-04-26 17:39:51Z srinath22 $
 */

#include "nccf_coord.h"
#include <stdlib.h>
#include <assert.h>
#include <netcdf.h>
#include <nccf_handle_error.h>

int main() {

  const int ndims = 2;
  const int dims[] = {10, 11};
  const double lonMin = 0.0;
  const double lonMax = 360.0;
  const double latMin = -90.0;
  const double latMax =  90.0;
  const int save = 1;
  double dLon, dLat;
  double *data, *lonData, *latData;
  int i, j;
  int lonId,latId;
  int status;

  const char *dimnames[] = {"ni", "nj"};

  data = (double *) malloc(sizeof(double) * dims[0] * dims[1]);

  /* Longitudes */
  dLon = (lonMax - lonMin) / (dims[0] - 1);
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      data[i + dims[0] * j] = lonMin + i * dLon;
    }
  }
  if ((status = nccf_def_lon_coord( ndims, dims, dimnames, data, save, 
					 &lonId))) ERR;

  /* Latitudes */
  dLat = (latMax - latMin) / (dims[1] - 1);
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      data[i + dims[0] * j] = latMin + j * dLat;
    }
  }
  if ((status = nccf_def_lat_coord( ndims, dims, dimnames, data, save, 
					 &latId))) ERR;

  /* Check that the lon and lat pointers are different */
  if ((status = nccf_get_coord_data_pointer(lonId, &lonData))) ERR;
  if ((status = nccf_get_coord_data_pointer(latId, &latData))) ERR;
  assert(lonData != latData);
  
  /* Free */

  if ((status = nccf_free_coord(lonId))) ERR;
  if ((status = nccf_free_coord(latId))) ERR;
  free(data);
  return status;
}
