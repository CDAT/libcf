/*
 * Test the definition of a structured grid.
 *
 * "$Id: tst_two_tiles.c 767 2011-06-06 23:20:19Z pletzer $"
 * */

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
#include "nccf_constants.h"

#include <assert.h>

int main(int argc, char *argv[]){

  printf("Running %s\n", argv[0]);

  /* Variables for the grid */
  // Each Id goes with a ncid
  int status, ncid;
  int iGrid; // Some counters
  int ndims = 2;         // Number of x, y dims
  int ngrids = 2;
  int nlon = 10;         // number of x faces, gives nlon + 1 edges
  int dims[] = {nlon, nlon};
  int nvertex = nlon * nlon;
  double period[] = {0., 0.};

  int gridid[ngrids], mosaicid, mosaicid2, i, j, globalId;
  int coordId[ngrids][ndims];
  const int save = 1;
  char coordinates_id[36+1];

  nc_type nc_mode = NC_CLOBBER;

  char *name = "tst_two_tiles";
//  char grid_filename[STRING_SIZE] = "tst_two_tiles_mosaic.nc";
  char gridname[STRING_SIZE];
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  double lon1[nvertex], lat1[nvertex], lon2[nvertex], lat2[nvertex], dlon, dlat;
  const double lonMin =   0.0;
  const double lonMax = 360.0;
  const double latMin = -90.0;
  const double latMax =  90.0;

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;  

  dlat = ( latMax - latMin )/( nlon - 1 );
  dlon = ( lonMax - lonMin )/( nlon - 1 )/2;

  for( j = 0; j < nlon; j++){
    for( i = 0; i < nlon; i++ ){
      lon1[j + nlon * i] = lonMin + j * dlon;
      lon2[j + nlon * i] = lonMin + j * dlon + 180;
      lat1[j + nlon * i] = latMin + i * dlat;
      lat2[j + nlon * i] = latMin + i * dlat;
    }
  }

/* Calculate the coordinates for the current tile */

/* Define tile 0 */
  iGrid = 0;
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, lon1, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, lat1, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
  sprintf( filename, "%s.nc", gridname );
  if(( status = nccf_def_grid( coordId[iGrid],
					  gridname, 
					  &gridid[iGrid] ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, gridname, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Write tile 0 */
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define tile 1 */
  iGrid = 1;
  if(( status = nccf_def_lon_coord( ndims, dims, dimnames, lon2, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims, dimnames, lat2, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
  sprintf( filename, "%s.nc", gridname );
  if(( status = nccf_def_grid( coordId[iGrid],
					  gridname, 
					  &gridid[iGrid] ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, gridname, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Write tile 1 */
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define the contacts */
  if(( status = nccf_def_mosaic( ngrids, gridid, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;

/* Write the file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;

/* write mosaic as polytopes */
//  nccf_print_mosaic_as_polytopes(mosaicid,"two_tiles.pt");

  status = nccf_free_mosaic( mosaicid );
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    status = nccf_free_grid( gridid[iGrid] );
    status = nccf_free_coord( coordId[iGrid][0] );
    status = nccf_free_coord( coordId[iGrid][1] );
  }  

/* Load the contacts, grids, and coordinates from file */
  char *tst_name  = "tst_two_tiles_mosaic";
  if ((status = nccf_def_mosaic_from_file(filename, 
					  tst_name,
   					  &mosaicid2))) ERR;

/* Check */
  int ncontacts;
  if ((status = nccf_inq_mosaic_ncontacts(mosaicid2, &ncontacts))) ERR;
  assert(ncontacts == 1);

/* free the grid */
  if(( status = nccf_free_mosaic( mosaicid ))) ERR;

  return 0;
}
