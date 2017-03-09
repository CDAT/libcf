/**
 * Test generation of a structured grid.
 *
 * "$Id: tst_1d_grid.c 767 2011-06-06 23:20:19Z pletzer $"
 *
 */

#include "nccf_grid.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <netcdf.h>

#include "nccf_global.h"
#include "libcf_src.h"
#include "nccf_coord.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(){

  const int save = 1;
  const int ndims = 1;
  const int dims[] = {10};
  const double lonMin[] = {0.0, 180.};
  const double lonMax[] = {180.0, 360.};
  double dLon;
  int nvertex = dims[0];
  double lonData[nvertex];
  int ngrids = 2;
  char coordinates_id[36+1];

  nc_type ncmode = NC_CLOBBER;
  int i, j;
  int coordids[ngrids][ndims];
  char **coordNames;
  int gridid[ngrids], gridid2, globalId;
  int status;
  int ncid;

  char name[STRING_SIZE];
  char gridname[STRING_SIZE];
  sprintf( name, "%s", "grid_1d");
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni"};

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Populate coordinates and create lon/lat coordinate objects */
  for( j = 0; j < ngrids; j++ ){
    dLon = (lonMax[j] - lonMin[j]) / (dims[0] - 1);
    for (i = 0; i < dims[0]; ++i) {
      lonData[i] = lonMin[j] + i * dLon;
    }
    if ((status = nccf_def_lon_coord(ndims, dims, dimnames, lonData, 
                              save, &coordids[j][0]))) ERR;

  /* Create grid */
    sprintf( gridname, "%s_grid%d", name, j );
    sprintf( filename, "%s.nc", gridname );
    if ((status = nccf_def_grid( coordids[j], gridname, &gridid[j]))) ERR;
    if ((status = nccf_def_global( &globalId ))) ERR;
    if ((status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                       coordinates_id, 0 ))) ERR;
    if ((status = nccf_add_global_att( globalId, CF_FILETYPE, 
                                       CF_GLATT_FILETYPE_GRID, 0))) ERR;
    if ((status = nccf_add_global_att( globalId, CF_GRIDNAME, 
                                       gridname, 0 ))) ERR;

  /* Write to file */
    if(( status = nc_create( filename, ncmode, &ncid ))) ERR;
    if(( status = nccf_put_grid(gridid[j], ncid))) ERR;
    if(( status = nccf_put_global(globalId, ncid))) ERR;
    if(( status = nc_close( ncid ))) ERR;
    if(( status = nccf_free_global( globalId ))) ERR;

  }

  /* Free */
  for ( j = 0; j < ngrids; j++ ){
    if (( status = nccf_free_grid(gridid[j]))) ERR;
    if (( status = nccf_free_coord( coordids[j][0] ))) ERR;
  }

  /* Read from file */
  coordNames = (char **) malloc(sizeof(char *) * ndims);
  coordNames[0] = (char *) calloc(strlen("lon") + 1, sizeof(char));
  strcpy(coordNames[0], "lon");

  int cid[ndims];
  for( j = 0; j < ngrids; j++ ){
    sprintf( filename, "%s_grid%d.nc", name, j );
    if(( status = nccf_def_grid_from_file(filename,
               ndims, (const char**)coordNames, name, 
               &gridid2))) ERR;
    status = nccf_inq_grid_coordids( gridid2, cid );
    if (( status = nccf_free_grid(gridid2))) ERR;
    for (i = 0; i < ndims; ++i) {
      status += nccf_free_coord( cid[i] );
    }
  }

  free(coordNames[0]);
  free(coordNames);
  return 0;
}
