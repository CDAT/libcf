/**
 * Test creation of curvilinear coordinates from axes
 *
 * $Id: tst_lonlat_from_axes.c 775 2011-06-13 20:36:52Z pletzer $
 */

#include <nccf_axis.h>
#include <nccf_coord.h>
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
  double dLon, dLat;
  int i;
  int lonId, latId;
  int status;
  // lon and lat axes need not be of the same type
  float *lonData;
  double *latData;
  int axisids[ndims];
  const int positive_up = 1;

  dLon = (lonMax - lonMin) / (dims[0] - 1);
  lonData = malloc(dims[0]*sizeof(float));
  for (i = 0; i < dims[0]; ++i) {
    lonData[i] = lonMin + dLon*i;
  }
  status = nccf_def_axis("longitude", dims[0], NC_FLOAT, 
                         lonData, "longitude", "degrees_east", NCCF_LONGITUDE,
                         "Y", positive_up, NULL, &axisids[0]); 
  if (status) ERR;
  free(lonData);

  dLat = (latMax - latMin) / (dims[1] - 1);
  latData = malloc(dims[1]*sizeof(double));
  for (i = 0; i < dims[1]; ++i) {
    latData[i] = latMin + dLat*i;
  }
  status = nccf_def_axis("latitude", dims[1], NC_DOUBLE, 
                         latData, "latitude", "degrees_north", NCCF_LATITUDE,
                         "X", positive_up, NULL, &axisids[1]);
  if (status) ERR;
  free(latData);

  status = nccf_def_coord_from_axes(ndims, axisids, 0, "lon", 
                                      "longitude", "degrees_east", &lonId);
  if (status) ERR;
  status = nccf_def_coord_from_axes(ndims, axisids, 1, "lat", 
                                      "latitude", "degrees_north", &latId);
  if (status) ERR;

  
  /* Free */

  if ((status = nccf_free_coord(lonId))) ERR;
  if ((status = nccf_free_coord(latId))) ERR;
  for (i = 0; i < ndims; ++i) {
    if ((status = nccf_free_axis(axisids[i]))) ERR;
  }
  return status;
}
