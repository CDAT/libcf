/**
 * Test building a grid with large dimension sizes.
 *
 * "$Id: tst_big_dimsizes.c 767 2011-06-06 23:20:19Z pletzer $"
 *
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "libcf_src.h"
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

int main(){

  const int save = 0;
  const int ndims = 2;

  /* Define a large grid */
  const int dims[] = {1000, 1000};

  const long nvertex = dims[0] * dims[1];
  const double lonMin = 0.0;
  const double lonMax = 360.0;
  const double latMin = -90.0;
  const double latMax = +90.0;
  double dLon, dLat;
  double *lonData;
  double *latData;
  double period[] = {360.0, 0.0};
  int gridid, globalId, mosaicid;
  int coordIds[ndims];
  char coordinates_id[36+1];

  lonData = (double *) malloc(nvertex * sizeof(double));
  latData = (double *) malloc(nvertex * sizeof(double));

  nc_type nc_mode = NC_CLOBBER;
  int i, j;
  int status;
  int ncid;
  const int ngrids = 1;

  char filename[STRING_SIZE];
  char name[STRING_SIZE] = "tst_sg_periodic_latlon";
  const char *dimnames[] = {"ni", "nj"};

  /* Create some coordinates */
  dLon = (lonMax - lonMin) / (dims[1] - 1);
  dLat = (latMax - latMin) / (dims[0] - 1);
  for (j = 0; j < dims[0]; ++j) {
    for (i = 0; i < dims[1]; ++i) {
      lonData[i + dims[1] * j] = lonMin + i * dLon;
      latData[i + dims[1] * j] = latMin + j * dLat;
    }
  }
  if ((status = nccf_def_lon_coord(ndims, dims, dimnames, lonData, save,
                &coordIds[0]))) ERR;
  if ((status = nccf_def_lat_coord(ndims, dims, dimnames, latData, save,
                &coordIds[1]))) ERR;

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Create grid */
  if ((status = 
       nccf_def_grid(coordIds, name, &gridid ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Write to file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_GRID );
  if ((status = nc_create(filename, nc_mode, &ncid))) ERR;
  if ((status = nccf_put_grid(gridid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if ((status = nc_close(ncid))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;
  
  /* Define the contact file */
  if(( status = nccf_def_mosaic( ngrids, &gridid, name, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

  /* Write to file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  if ((status = nc_create(filename, nc_mode, &ncid))) ERR;
  if ((status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if (( status = nccf_put_global(globalId, ncid)));
  if ((status = nc_close(ncid))) ERR;

  /* Free */
  if (( status = nccf_free_global( globalId ))) ERR;
  if (( status = nccf_free_mosaic( mosaicid ))) ERR;
  if (( status = nccf_free_grid( gridid ))) ERR;
  if (( status = nccf_free_coord( coordIds[0] ))) ERR;
  if (( status = nccf_free_coord( coordIds[1] ))) ERR;
  free(lonData);
  free(latData);

  return 0;
}
