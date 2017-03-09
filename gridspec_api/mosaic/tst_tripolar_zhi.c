/*
 * Test the creation of a tripolar grid
 *
 * "$Id: tst_tripolar_zhi.c 767 2011-06-06 23:20:19Z pletzer $"
 * */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libcf_src.h>
#include <netcdf.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(int arg, char *argv[]){

/* Grid variables */
  int ndims = 2;
  int ngrids = 2;
  int gldims[] = {10, 10};
  int glnvertex = gldims[0] * gldims[1];
  const char *dimnames[] = {"ni", "nj"};
  int capIndx = 7;
  double gllons[glnvertex], gllats[glnvertex];
  int j, k, status, ncid, mosaicid;
  const int save = 1;
  char coordinates_id[36+1];

  printf("Running %s\n", argv[0]);

  double period[] = {0, 0};
// 4 for the lower grid, and 2 for the bipolar grid
  int gridid[ngrids], coordIds[2][2], globalId;

  nc_type nc_mode = NC_CLOBBER;

  char *tname = "2tile_tripolar_zhi";
  char *name, filename[STRING_SIZE];
  name = ( char * )calloc( STRING_SIZE, sizeof( char ) );
  sprintf( name, "%s", tname );

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  for (k=0; k<2; ++k) {
  /* Get the first half */
    status = nccf_get_tripolar_halfgrid( gldims, k, capIndx, gllons, gllats );

    if ((status =
      nccf_def_lon_coord(ndims, gldims, dimnames, gllons, save, &coordIds[k][0]))) ERR;
    if ((status =
      nccf_def_lat_coord(ndims, gldims, dimnames, gllats, save, &coordIds[k][1]))) ERR;

  /* Write the grid */
    sprintf( filename, "%s_grid%d.nc", name, k );
    status = nccf_def_grid( coordIds[k], name, &gridid[k] );
    if(( status = nccf_def_global( &globalId ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

    if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
    if(( status = nccf_put_grid(gridid[k], ncid))) ERR;
    if(( status = nccf_put_global(globalId, ncid))) ERR;
    if(( status = nc_close( ncid ))) ERR;
    if(( status = nccf_free_global( globalId ))) ERR;
  }
/* contacts */
  if(( status = nccf_def_mosaic( ngrids, gridid, name, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;
  sprintf( filename, "%s_contacts.nc", name );
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

/* write mosaic as polytopes */
//  nccf_print_mosaic_as_polytopes(mosaicid,"tripolar_zhi.pt");

  if(( status = nccf_free_mosaic( mosaicid ))) ERR;
  free(name);

  for (k=0; k<ngrids; ++k) {
    nccf_free_grid(gridid[k]);
    for (j=0; j<ndims; ++j)
        if(( status = nccf_free_coord(coordIds[k][j]))) ERR;
  }

  return NC_NOERR;
}
