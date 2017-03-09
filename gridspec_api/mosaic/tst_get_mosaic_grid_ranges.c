/*
 * Test grid contact maps.
 *
 * "$Id: tst_get_mosaic_grid_ranges.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>

#include <libcf_src.h>
#include <netcdf.h>
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_varObj.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

#include <assert.h>

int main(  ){

  /* Initialize the the cflists */

  /* Variables for the grid */
  // Each Id goes with a ncid
  int status, ncid;
  int iGrid; // Some counters
  const int ndims = 2;         // Number of x, y dims
  const int ngrids = 2;
  const int nx   = 20;         // number of x faces, gives nlon + 1 edges
  const int ny   = 10;         // number of x faces, gives nlon + 1 edges
  const int dims[] = {ny, nx};
  const int nvertex = ny * nx;
  /* Return contact indices for i, j etc. and gridids */
  int ibeg0[ndims], iend0[ndims], ibeg1[ndims], iend1[ndims], g0, g1;
  const double period[] = {0., 0.};

  int gridid[ngrids], mosaicid, i, j, globalId;
  int coordId[ngrids][ndims];
  int save = 1;
  char coordinates_id[36+1];

  nc_type nc_mode = NC_CLOBBER;

  char *name = "tile";
  char grid_filename[STRING_SIZE] = "tst_get_mosaic_grid_ranges.nc";
  char gridname[STRING_SIZE];
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  double x1[nvertex], y1[nvertex], x2[nvertex], y2[nvertex], dx, dy;
  const double xMin =   0.0;
  const double xMax = 360.0;
  const double yMin = -90.0;
  const double yMax =  90.0;

  dy = ( yMax - yMin )/( ny - 1 );
  dx = ( xMax - xMin )/( nx - 1 )/2;

  for( j = 0; j < nx; j++){
    for( i = 0; i < ny; i++ ){
      x1[j + nx * i] = xMin + j * dx;
      x2[j + nx * i] = xMin + j * dx + 180;
      y1[j + nx * i] = yMin + i * dy;
      y2[j + nx * i] = yMin + i * dy;
    }
  }

  /* Define tile 0 */
  iGrid = 0;
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, x1, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, y1, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%d", name, iGrid + 1 );
  sprintf( filename, "%s.nc", gridname );

  if(( status = nccf_def_grid( coordId[iGrid],
            gridname, 
            &gridid[iGrid] ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Write tile 0 */
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

  /* Define tile 1 */
  iGrid = 1;
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, x2, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, y2, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%d", name, iGrid + 1 );
  sprintf( filename, "%s.nc", gridname );
  if(( status = nccf_def_grid( coordId[iGrid],
            gridname, 
            &gridid[iGrid] ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Write tile 1 */
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

  /* Define the contacts */
  if(( status = nccf_def_mosaic( ngrids, gridid, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

  /* Write the file */
  if(( status = nc_create( grid_filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

  /* Get the contact map and tile contacts */
  nccf_inq_mosaic_gridranges( mosaicid, 0, &g0, &g1, 
			       ibeg0, iend0, 
			       ibeg1, iend1 );

  /* free the grid */
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    if(( status = nccf_inq_grid_coordids( gridid[iGrid],
        coordId[iGrid] ))) ERR;
    for( i = 0; i < ndims; i++ ){
      if(( status = nccf_free_coord( coordId[iGrid][i] ))) ERR;
    }
    if(( status = nccf_free_grid( gridid[iGrid] ))) ERR;
  }
  if(( status = nccf_free_mosaic( mosaicid ))) ERR;

  assert( g0 == 0 );
  assert( g1 == 1 );
  assert( ibeg0[0] ==  0 );
  assert( ibeg0[1] == 19 );
  assert( iend0[0] ==  9 );
  assert( iend0[1] == 19 );
  assert( ibeg1[0] == 0  );
  assert( ibeg1[1] == 0  );
  assert( iend1[0] ==  9 );
  assert( iend1[1] ==  0 );

/* Close the file*/
  return 0;
}
