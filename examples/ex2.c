/**
 * Create a host file aggregation comprising static and time dependent
 * data stored in different files. The time files contain multiple time
 * steps. Time dependent and static data exist on multiple tiles.
 *
 * The grids are created in ex_create_cube_grids.c
 * The static data are created by ex_create_static_data.c
 * The time data are created by ex_create_time_data.c
 *
 * "$Id: ex2.c 851 2011-11-08 14:37:20Z pletzer $"
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

int main(  ){

  int status = NC_NOERR;
  int i, j, ndims;

  /* create the coordinate id */
  char coordinates_id[36+1];
  char data_id[36+1];
  const int seed = 12345;
  if ((status += nccf_generate_id(seed, coordinates_id))) ERR;
  if ((status += nccf_generate_id(seed, data_id))) ERR;

  /* Dimensions */
  const int nGrids = 3;
  int nvars = 2;
  int nTimes = 2;
  int nTimesPerFile = 3;
  int nTimeFiles = nGrids * nTimes * nvars;
  int nStatFiles = nGrids;

  const int nCells = 10;

  struct time_struct local;

  local.nVars   = nvars;
  local.nCells  = nCells;
  local.nTimes  = nTimes;
  local.nTimesPerFile = nTimesPerFile;
  local.nGrids  = nGrids;
  local.value   = 273.15;
  local.varname = NAME_TIME_DATA_VARN_t;
  local.st_name = "temperature";
  local.units   = "K";
  local.time_stanname = "Time";
  local.time_longname = "Time";
  local.time_units    = "days since 2000-01-01";

  /* create the grids */
  char **gridFiles   = calloc(nGrids, sizeof(char *));
  char **staticFiles = calloc(nStatFiles, sizeof(char *));
  char **timeFiles   = calloc(nTimeFiles, sizeof(char *));
  for (i = 0; i < nGrids; ++i) {
    gridFiles[i]   = calloc(STRING_SIZE, sizeof(char));
  }
  for( i = 0; i < nStatFiles; i++ ) staticFiles[i] = calloc(STRING_SIZE, sizeof(char));
  for( i = 0; i < nTimeFiles; i++ ) timeFiles[i]   = calloc(STRING_SIZE, sizeof(char));
  int gridids[nGrids];
  int staticids[nStatFiles];
  int timeids[nTimeFiles];
  ex_create_cube_grids(coordinates_id, nCells, nGrids, gridFiles,
                       gridids );

  /* create the mosaic */
  int mosaicid, globalid;
  const double periods[] = {0., 0., 0.};
  status += nccf_def_mosaic(nGrids, gridids, NAME_MOSAIC, &mosaicid);
  status += nccf_compute_mosaic_contacts( mosaicid, periods );
  if (status != NC_NOERR) ERR;
  status += nccf_def_global( &globalid );
  status += nccf_add_global_att(globalid, CF_COORDINATES_ID, coordinates_id, 0);
  status += nccf_add_global_att(globalid, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0);
  if (status != NC_NOERR) ERR;
  int ncid;
  status += nc_create(NAME_MOSAIC_FILE, NC_CLOBBER, &ncid);
  if (status != NC_NOERR) ERR;
  status += nccf_put_mosaic(mosaicid, ncid);
  status += nccf_put_global(globalid, ncid);
  status += nc_close(ncid);
  if (status != NC_NOERR) ERR;
  status += nccf_free_global(globalid);

  /* create static data */
  status += ex_create_static_data( data_id, coordinates_id, nCells, nStatFiles,
                                   nGrids, gridids, staticFiles, staticids );

 /* create time dependent data */

  status = ex_create_time_data( data_id, coordinates_id, 0, &local, gridids, timeFiles, timeids);
  if( status ) ERR;
  local.value   = 10.0;
  local.varname = NAME_TIME_DATA_VARN_v;
  local.st_name = "velocity";
  local.units = "m/s";
  status = ex_create_time_data( data_id, coordinates_id, 1, &local, gridids, timeFiles, timeids);
  if( status ) ERR;

  /* create the host file */
  int hostid;
  if (( status = nccf_def_host(coordinates_id, data_id, nTimes, &hostid))) ERR;
  if (( status += nc_create( NAME_HOST_FILE, NC_CLOBBER, &ncid ))) ERR;
  for( i = 0; i < nGrids; i++ ){
    if (( status += nccf_add_host_file( hostid, gridFiles[i], 0 ))) ERR;
  }
  for( i = 0; i < nStatFiles; i++ ){
    if (( status += nccf_add_host_file( hostid, staticFiles[i], 0 ))) ERR;
  }
  for( i = 0; i < nTimeFiles; i++ ){
    if (( status += nccf_add_host_file( hostid, timeFiles[i], 0 ))) ERR;
  }
  if (( status += nccf_add_host_file( hostid, NAME_MOSAIC_FILE, 0 ))) ERR;

  /* Define the global attributes */
  if(( status = nccf_def_global( &globalid ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_FILETYPE,
                                     CF_GLATT_FILETYPE_HOST, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_DATA_ID,
                                     data_id , 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID,
                                     coordinates_id, 0 ))) ERR;

  if (( status = nccf_put_host(hostid, ncid))) ERR;
  if (( status = nccf_put_global(globalid, ncid))) ERR;
  if (( status = nc_close( ncid ))) ERR;
  if (( status = nccf_free_global( globalid ))) ERR;

  /* Build data from host */
  int gfindx = 0; // tile  index
  int read_data = 1; 
  int dataidT;
  if (( status = nccf_def_host_data(hostid, "T", gfindx, 
                                    read_data, &dataidT) )) ERR;

  /* Check that we can access the data */
  nc_type xtype;
  void *data_ptrT;
  const void *fill_value;
  if (( status = nccf_get_data_pointer(dataidT, &xtype, 
                                       &data_ptrT, &fill_value) )) ERR;
  if (( status = nccf_inq_data_ndims(dataidT, &ndims) )) ERR;
  int dimsT[ndims];
  if (( status = nccf_inq_data_dims(dataidT, dimsT) )) ERR;
  int ntot = 1;
  for (i = 0; i < ndims; ++i) ntot *= dimsT[i];
#if 0
  if (xtype == NC_FLOAT) {
    float *data = (float *) data_ptrT;
    for (i = 0; i < ntot; ++i) {
      printf("i = %d float data = %f\n", i, data[i]);
    }
  }
  else if (xtype == NC_DOUBLE) {
    double *data = (double *) data_ptrT;
    for (i = 0; i < ntot; ++i) {
      printf("i = %d double data = %f\n", i, data[i]);
    }
  }
  else {
    printf("unsupported data type.\n");
    return 1;
  }
#endif

  /* Clean up data and unlerlying grid and coordinates */
  int grididT;
  if (( status = nccf_inq_data_gridid(dataidT, &grididT) )) ERR;
  if (( status = nccf_inq_grid_ndims(grididT, &ndims) )) ERR;
  int coordidsT[ndims];
  if (( status = nccf_inq_grid_coordids(grididT, coordidsT) )) ERR;
  if (( status = nccf_free_data(dataidT) )) ERR;
  if (( status = nccf_free_grid(grididT) )) ERR;
  for (i = 0; i < ndims; ++i) {
    if (( status = nccf_free_coord(coordidsT[i]) )) ERR;
  }
  
  /* Generate target grid for SCRIP interpolation */
  int lonlat_dims[] = {90, 180};
  int lonlat_coordids[2];
  int lonlat_gridid;
  const double coord_mins[] = {  0.0, -90.0};
  const double coord_maxs[] = {360.0, +90.0};
  ex_create_lonlat(lonlat_dims, coord_mins, coord_maxs,
                   lonlat_coordids, &lonlat_gridid);

  /* clean up */
  if ((status = nccf_free_grid(lonlat_gridid))) ERR;
  if ((status = nccf_free_coord(lonlat_coordids[0]))) ERR;
  if ((status = nccf_free_coord(lonlat_coordids[1]))) ERR;

  if ((status = nccf_free_host(hostid))) ERR;
  if ((status = nccf_free_mosaic(mosaicid))) ERR;

  if ((status = nccf_inq_grid_ndims(gridids[0], &ndims))) ERR;
  int cid[ndims];
  for (i = 0; i < nGrids; ++i) {
    status += nccf_inq_grid_coordids( gridids[i], cid );

    for( j = 0; j < ndims; j++ )
      if ((status += nccf_free_coord(cid[j]))) ERR;
    if ((status += nccf_free_grid(gridids[i]))) ERR;
  }

  for( i = 0; i < nStatFiles; i++ ){
    if ((status += nccf_free_data(staticids[i]))) ERR;
    free( staticFiles[i] );
  }
  for( i = 0; i < nTimeFiles; i++ ){
    if ((status += nccf_free_data(timeids[i]))) ERR;
    free( timeFiles[i] );
  }
  for (i = 0; i < 2*nGrids; ++i) {
  }

  for (i = 0; i < nGrids; ++i) {
    free(gridFiles[i]);
  }

  free(gridFiles);
  free(staticFiles);
  free(timeFiles);

  return status;
}
