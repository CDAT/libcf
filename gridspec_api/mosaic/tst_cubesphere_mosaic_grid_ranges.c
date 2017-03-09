/**
 * Test reading the grid contacts for a cube sphere from a mosaic file.
 *
 * "$Id: tst_cubesphere_mosaic_grid_ranges.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libcf_src.h>
#include <netcdf.h>
#include "nccf_coord.h"
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(){

  /* to save the data in the coordinate object */
  const int save = 1;

  /* Variables for the grid */
  int ncid, mosaicid, globalId;
  int status;
  int iGrid; // Some counters

  const int ndims = 2;         // Number of x, y dims
  const int ngrids = 6;
  const int nlon = 10;         // number of x faces, gives nlon + 1 edges
  const int elon = nlon + 1;   // # of edges
  const int dims[] = {elon, elon};
  const int nvertex = elon*elon;
  int faceVec[3] = {0, 0, 0}, pos, sign;
  int i;
  char coordinates_id[36+1];

  nc_type nc_mode = NC_CLOBBER;

  double clon[nvertex], clat[nvertex];
  double period[] = {0., 0.};

  char *name = "tst_cubesphere_periodic";
  char gridname[STRING_SIZE], filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  int coordIds[ngrids][ndims];
  int gridids[ngrids];

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Define a grid */
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    pos = iGrid % 3;
    sign = ( iGrid/3 )  == 0 ? 1 : -1;
    faceVec[pos] = sign;

    nccf_get_cubedsphere_grid( dims, faceVec, clon, clat );
    faceVec[pos] = 0; // reset to zero.

    if(( status = nccf_def_lon_coord( ndims, dims, dimnames, clon, save,
                                            &coordIds[iGrid][0] ))) ERR;
    if(( status = nccf_def_lat_coord( ndims, dims, dimnames, clat, save,
                                            &coordIds[iGrid][1] ))) ERR;

    /* Define the grid for the current tile. */
    sprintf( gridname, "%s%d", name, iGrid );
    if(( status = nccf_def_grid( coordIds[iGrid],
					    gridname, 
					    &gridids[iGrid] ))) ERR;
    if(( status = nccf_def_global( &globalId ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

    sprintf( filename, "%s%s%d.nc", name, CF_FILENAME_GRID, iGrid );
    if(( status = nc_create( filename, nc_mode, &ncid )));
    if(( status = nccf_put_grid(gridids[iGrid], ncid)));
    if(( status = nccf_put_global(globalId, ncid)));
    if(( status = nc_close( ncid )));
    if(( status = nccf_free_global( globalId )));

  }

  /* Write the mosaic file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  if(( status = nccf_def_mosaic( ngrids, gridids, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if(( status = nc_close( ncid ))) ERR;

  /* Find the contacts */
  int g0, g1, ij0_min[ndims], ij0_max[ndims], ij1_min[ndims], ij1_max[ndims];
  int ncont, iContact;
  nccf_inq_mosaic_ncontacts( mosaicid, &ncont );

  for( iContact = 0; iContact < ncont; iContact++ ){
    status = nccf_inq_mosaic_gridranges( 
          mosaicid, iContact, &g0, &g1, ij0_min, ij0_max, ij1_min, ij1_max );
  }

  /* free the contacts object */
  if(( status = nccf_free_global( globalId ))) ERR;
  if(( status = nccf_free_mosaic( mosaicid ))) ERR;

  /* free the grid and the coordinates */
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    if(( status = nccf_free_grid( gridids[iGrid] ))) ERR;
    for (i = 0; i < ndims; ++i) {
      if(( status = nccf_free_coord( coordIds[iGrid][i] ))) ERR;
    }
  }

/* Close the file*/
  return 0;
}
