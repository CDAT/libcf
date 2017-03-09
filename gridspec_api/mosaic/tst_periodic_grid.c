/**
 * Test the definition of a structured grid.
 *
 * "$Id: tst_periodic_grid.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <netcdf.h>
#include <libcf_src.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_varObj.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(  ){

  /* whether or not coordinate data are copy-saved in coordinate object */
  const int save = 1;
  /* periodicities used in determining contacts */
  const double period[] = {360., 0.};
  /* number of space dimensions */
  const int ndims = 2;
  /* number of grids */
  const int ngrids = 2;
  /* unique grid id */
  char coordinates_id[36+1];

  int i, j, iGrid;
  int status, ncid;
  int nx = 10, ny = 10; // number of x faces, gives nlon + 1 edges
  int dims[] = {ny, nx};
  int nvertex = dims[0] * dims[1];

  /* Coordinate ids */
  int coordId[ngrids][ndims];
  int gridid[ngrids], globalId, mosaicid;


  nc_type nc_mode = NC_CLOBBER;

  char *name = "tst_periodic_latlon";
  char gridname[STRING_SIZE];
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  double x1[nvertex], y1[nvertex], x2[nvertex], y2[nvertex], dx, dy;
  const double xMin =   0.0;
  const double xMax = 360.0;
  const double yMin = -90.0;
  const double yMax =  90.0;
  
  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;  

  dy = ( yMax - yMin )/( ny - 1 );
  dx = ( xMax - xMin )/( nx - 1 )/2;

  for( j = 0; j < nx; j++){
    for( i = 0; i < ny; i++ ){
      x1[j + nx * i] = xMin + j * dx;
      x2[j + nx * i] = xMin + j * dx + 180;
      y1[j + ny * i] = yMin + i * dy;
      y2[j + ny * i] = yMin + i * dy;
    }
  }

/* Calculate the coordinates for the current tile */

/* Define tile 0 */
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, x1, save,
                                         &coordId[0][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, y1, save,
                                         &coordId[0][1] ))) ERR;

  iGrid = 0;
  sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
  if(( status = nccf_def_grid( coordId[0],
					  gridname, 
					  &gridid[iGrid] ))) ERR;

  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  sprintf( filename, "%s.nc", gridname );
  if(( status = nc_create( filename, nc_mode, &ncid )));
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define tile 1 */
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, x2, save,
                                         &coordId[1][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, y2, save,
                                         &coordId[1][1] ))) ERR;

  iGrid = 1;
  sprintf( gridname, "%s%s%d.nc", name, CF_FILENAME_GRID, iGrid );
  if(( status = nccf_def_grid( coordId[1],
					  gridname, 
					  &gridid[iGrid] ))) ERR;

  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  sprintf( filename, "%s.nc", gridname );
  if(( status = nc_create( filename, nc_mode, &ncid )));
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define the contacts */
  if(( status = nccf_def_mosaic( ngrids, gridid, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

/* Write Polytopes */
//  if(( status = nccf_print_mosaic_as_polytopes( mosaicid, "periodic.pt" ))) ERR;

/* Write the contacts file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Checks */
  int ncontacts;
  if(( status = nccf_inq_mosaic_ncontacts(mosaicid, &ncontacts))) ERR;
  assert(ncontacts == 2);

  //if( status=nccf_print_mosaic_as_polytopes(mosaicid,"periodic_grid.pt")) ERR;

/* Free */
  for (j = 0; j < ngrids; ++j) {
    for (i = 0; i < ndims; ++i) {
      if(( status = nccf_free_coord(coordId[j][i]))) ERR;
    }
  }
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    if(( status = nccf_inq_grid_coordids( gridid[iGrid],
        coordId[iGrid] ))) ERR;
    if(( status = nccf_free_grid( gridid[iGrid] ))) ERR;
  }

  if(( status = nccf_free_mosaic( mosaicid ))) ERR;

  return 0;
}
