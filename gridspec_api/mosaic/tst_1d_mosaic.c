/**
 * Test generation of a structured grid and write a mosaic file.
 *
 * "$Id: tst_1d_mosaic.c 767 2011-06-06 23:20:19Z pletzer $"
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

  const int save = 1;
  const int ndims = 1;
  const int dims[] = {10};
  const double lonMin[] = {0.0, 180.};
  const double lonMax[] = {180.0, 360.};
  double dLon;
  int nvertex = dims[0];
  double lonData[nvertex];
  int ngrids = 2;
  double period[] = {0.};
  char coordinates_id[36+1];

  nc_type nc_mode = NC_CLOBBER;
  int i, j;
  int coordIds[ngrids][ndims], gridid[ngrids], mosaicid, globalId;
  int status, ncid;

  char name[STRING_SIZE];
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
                    save, &coordIds[j][0]))) ERR;

    /* Create grid */
    sprintf( name, "%s_grid%d", "grid_1d", j);
    sprintf( filename, "%s.nc", name );
    if ((status = nccf_def_grid(coordIds[j], name, &gridid[j]))) ERR;
    if(( status = nccf_def_global( &globalId ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

    /* Write to file */
    if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
    if(( status = nccf_put_grid(gridid[j], ncid))) ERR;
    if(( status = nccf_put_global(globalId, ncid))) ERR;
    if(( status = nc_close( ncid ))) ERR;

    if(( status = nccf_free_global( globalId ))) ERR;

  }

  /* Write the mosaic file */
  ncid = -1;
  if(( status = nccf_def_mosaic( ngrids, gridid, name, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  sprintf( filename, "%s%s.nc", "test_", CF_FILENAME_MOSAIC );
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;

  /* write mosaic as polytopes */
//  nccf_print_mosaic_as_polytopes(mosaicid,"1d_mosaic.pt");

  /* Free */
  for ( j = 0; j < ngrids; j++ ){
    if (( status = nccf_free_grid(gridid[j]))) ERR;
    if (( status = nccf_free_coord( coordIds[j][0] ))) ERR;
  }

  if(( status = nccf_free_global( globalId ))) ERR;
  printf( "Befor free\n" );
  if(( status = nccf_free_mosaic( mosaicid ))) ERR;
  printf( "After free\n" );

  return 0;

}
