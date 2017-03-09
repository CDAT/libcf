/**
 * Test generation of a simple lon-lat structured grid.
 *
 * "$Id: tst_struct_grid_latlon.c 893 2011-12-21 22:15:08Z pletzer $"
 *
 */

#include "nccf_grid.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include <libcf_src.h>

#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(){

  const int save = 0;
  const int ndims = 2;
  const int dims[] = {10, 11};
  const double lonMin = 0.0;
  const double lonMax = 360.0;
  const double latMin = -90.0;
  const double latMax = +90.0;
  double dLon, dLat;
  int nvertex = dims[0] * dims[1];
  double lonData[nvertex], latData[nvertex];
  int imask[nvertex]; 
  char coordinates_id[36+1];
  char coordinates_id2[36+1];

  nc_type ncmode = NC_CLOBBER;
  int i, j;
  int coordids[ndims];
  char **coordNames;
  int gridid, gridid2, globalId;
  int status;
  int ncid;

  char name[STRING_SIZE];
  sprintf( name, "%s", "lonlatgrid");
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Populate coordinates and create lon/lat coordinate objects */
  dLon = (lonMax - lonMin) / (dims[0] - 1);
  dLat = (latMax - latMin) / (dims[1] - 1);
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      lonData[i + dims[0] * j] = lonMin + i * dLon;
      latData[i + dims[0] * j] = latMin + j * dLat;
      imask[i + dims[0] * j] = 1; // 1=valid
    }
  }
  if ((status = nccf_def_lon_coord(ndims, dims, dimnames, lonData, save, 
            &coordids[0]))) ERR;
  if ((status = nccf_def_lat_coord(ndims, dims, dimnames, latData, save, 
            &coordids[1]))) ERR;

  /* Create grid */
  sprintf( filename, "%s_grid.nc", name );
  if ((status = nccf_def_grid( coordids, name, &gridid))) ERR;
  if ((status = nccf_set_grid_validmask(gridid, imask))) ERR;
  if ((status = nccf_def_global( &globalId ))) ERR;
  if ((status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                     coordinates_id, 0 ))) ERR;
  if ((status = nccf_add_global_att( globalId, CF_FILETYPE,   
                                     CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if ((status = nccf_add_global_att( globalId, CF_GRIDNAME, 
                                     name, 0 ))) ERR;

  double coord_periodicity[ndims];
  if ((status = nccf_inq_grid_periodicity(gridid, coord_periodicity))) ERR;
  assert(coord_periodicity[0] < CF_HUGE_DOUBLE);  // periodic
  assert(coord_periodicity[1] == CF_HUGE_DOUBLE); // non-periodic
 
  /* Write to file */
  if(( status = nc_create( filename, ncmode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;

  /* Free */
  if (( status = nccf_free_grid(gridid))) ERR;
  if (( status = nccf_free_global( globalId ))) ERR;
  if (( status = nccf_free_coord( coordids[0] ))) ERR;
  if (( status = nccf_free_coord( coordids[1] ))) ERR;

  /* Read from file */
  coordNames = (char **) malloc(sizeof(char *) * ndims);
  coordNames[0] = (char *) calloc(strlen("lon") + 1, sizeof(char));
  strcpy(coordNames[0], "lon");
  coordNames[1] = (char *) calloc(strlen("lat") + 1, sizeof(char));
  strcpy(coordNames[1], "lat");

  if(( status = nccf_def_grid_from_file(filename, 
               ndims, (const char**)coordNames, name, 
               &gridid2))) ERR;

  if ((status = nccf_save_grid_scrip(gridid2, 
					    "tst_struct_grid_latlon_scrip.nc"))) ERR;

  /* Test the global attribute writer */
  if((status = nccf_def_global_from_file( filename, &globalId ))) ERR;
  nccf_inq_global_att( globalId, CF_COORDINATES_ID, coordinates_id2 );
  assert( strcmp( coordinates_id, coordinates_id2 ) == 0 );
  nccf_free_global( globalId );

  for (i = 0; i < ndims; ++i) {
    free(coordNames[i]);
  }
  free(coordNames);

  int cid[ndims];
  status = nccf_inq_grid_coordids( gridid2, cid );
  if (( status = nccf_free_grid(gridid2))) ERR;
  for( i = 0; i < ndims; i++ ) status += nccf_free_coord( cid[i] );

  return 0;
}
