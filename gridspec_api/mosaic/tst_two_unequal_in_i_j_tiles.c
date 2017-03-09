/**
 * Test mosaic construction of two tiles with unequal i and j dimensions.
 *
 * "$Id: tst_two_unequal_in_i_j_tiles.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>

#include <libcf_src.h>
#include <netcdf.h>
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

#include <assert.h>

int main(  ){

  /* Variables for the grid */
  // Each Id goes with a ncid
  int status, ncid;
  int iGrid; // Some counters
  const int ndims = 2;         // Number of x, y dims
  const int ngrids = 2;
  const int dims0[] = { 5,10};
  const int dims1[] = { 11, 7};
  const double period[] = {360., 0.};
  char coordinates_id[36+1];

  int gridid[ngrids], mosaicid, mosaicid2, i, j, globalId;
  int coordId[ngrids][ndims], save = 1;

  nc_type nc_mode = NC_CLOBBER;

  char *name = "tst_two_unequal_in_i_j";
  char gridname[STRING_SIZE];
  char filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};
  char **files;
  files = ( char** )calloc( 3, sizeof( char* ));
  for( iGrid = 0; iGrid < 3; iGrid++ )
      files[iGrid] = ( char* )calloc( STRING_SIZE, sizeof( char ));

  int nv0 = dims0[0] * dims0[1];
  int nv1 = dims1[0] * dims1[1];
  double lon1[nv0], lat1[nv0], lon2[nv1], lat2[nv1], dlon, dlat;
  const double lonMin[2] = {  0., 270.};
  const double lonMax[2] = {270., 360.};
  const double latMin    = -90.0;
  const double latMax    =  90.0;

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Grid 0 */
  dlat = ( latMax    - latMin )    / ( dims0[0] - 1 );
  dlon = ( lonMax[0] - lonMin[0] ) / ( dims0[1] - 1 );

  for( j = 0; j < dims0[0]; j++ ){
    for( i = 0; i < dims0[1]; i++){
      lon1[i + dims0[1] * j] = lonMin[0] + i * dlon;
      lat1[i + dims0[1] * j] = latMin    + j * dlat;
    }
  }

  /* Grid 1 */
  dlat = ( latMax    - latMin )    / ( dims1[0] - 1 );
  dlon = ( lonMax[1] - lonMin[1] ) / ( dims1[1] - 1 );

  for( j = 0; j < dims1[0]; j++ ){
    for( i = 0; i < dims1[1]; i++){
      lon2[i + dims1[1] * j] = lonMin[1] + i * dlon;
      lat2[i + dims1[1] * j] = latMin    + j * dlat;
    }
  }

/* Calculate the coordinates for the current tile */

/* Define tile 0 */
  iGrid = 0;
  if(( status = nccf_def_lon_coord( ndims, dims0, dimnames, lon1, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims0, dimnames, lat1, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
  sprintf( filename, "%s.nc", gridname );
  if(( status = nccf_def_grid( coordId[iGrid],
					  gridname, 
					  &gridid[iGrid] ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;
  sprintf( files[0], "%s", filename );

/* Write tile 0 */
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;  
  if(( status = nccf_put_grid(gridid[iGrid], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define tile 1 */
  iGrid = 1;
  if(( status = nccf_def_lon_coord( ndims, dims1, dimnames, lon2, save,
                                         &coordId[iGrid][0] ))) ERR;
  if(( status = nccf_def_lat_coord( ndims, dims1, dimnames, lat2, save,
                                         &coordId[iGrid][1] ))) ERR;

  sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
  sprintf( filename, "%s.nc", gridname );
  sprintf( files[1], "%s", filename );
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
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* Define the contacts */
  if(( status = nccf_def_mosaic( ngrids, gridid, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

/* Write the file */
  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  sprintf( files[2], "%s", filename );
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* write mosaic as polytopes */
//  nccf_print_mosaic_as_polytopes(mosaicid,"two_unequal_in_i_j.pt");

  status = nccf_free_mosaic( mosaicid );
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    status = nccf_free_grid( gridid[iGrid] );
    status = nccf_free_coord( coordId[iGrid][0] );
    status = nccf_free_coord( coordId[iGrid][1] );
  }  

/* Load the contacts, grids, and coordinates from file */
  if ((status = nccf_def_mosaic_from_file( filename,
					    "tst_two_tiles_mosaic",
					    &mosaicid2))) ERR;

/* Check */
  int ncontacts;
  if ((status = nccf_inq_mosaic_ncontacts(mosaicid2, &ncontacts))) ERR;
  assert(ncontacts == 2);

  char contact_map[STRING_SIZE], tile_contact[STRING_SIZE], sep[10];
  nccf_inq_mosaic_contactmap(  mosaicid2, 0, contact_map );
  nccf_inq_mosaic_tilecontact( mosaicid2, 0, tile_contact );
  nccf_inq_mosaic_tileseparator( sep );
  printf( "%s, %s, %s\n", contact_map, tile_contact, sep );


/* free the grid */
  if(( status = nccf_free_mosaic( mosaicid2 ))) ERR;

/* Add the signature to these files */
  char *id;
  char *names[] = {"grid0", "grid1", "mosaic"};
  char *types[] = {"grid", "grid", "mosaic"};
  id = ( char* )calloc( STRING_SIZE, sizeof( char ));

  nccf_generate_id( 77, id );
  if(( status = nccf_add_id_to_files( id, 3, (const char**)names, 
                                             (const char**)types, 
                                             (const char**)files ))) ERR;

  free( id );
  for (i = 0; i < 3; ++i) {
    free( files[i] );
  }
  free(files);

/* Close the file*/
  return 0;
}
