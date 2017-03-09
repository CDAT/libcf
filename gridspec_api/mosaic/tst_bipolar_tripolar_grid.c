/*
 * Test the creation of a tripolar grid and write a mosaic file for it.
 * \author Dave Kinding and Alexander Pletzer, Tech-X Corp.
 *
 * "$Id: tst_bipolar_tripolar_grid.c 767 2011-06-06 23:20:19Z pletzer $"
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

int main(){

/* Grid variables */
  const int ndims = 2;
  const int ngrids = 3;
  const int bpdims[] = {6, 6}, gldims[] = {6, 6};
  const int bpnvertex = bpdims[0] * bpdims[1];
  const int glnvertex = gldims[0] * gldims[1];
  const double latPerim = 60;
  const double lonSing = 0.0;
  double bplons[bpnvertex], bplats[bpnvertex];
  double gllons[glnvertex], gllats[glnvertex];
  double dLon, dLat, lonMin[] = {-180, 0}, lonMax[] = {0, 180};
  double latMin = -90, latMax = latPerim;
  int i, j, k, status, ncid;
  const int save = 1;
  char coordinates_id[36+1];  
  
  double period[] = {0., 0.};

// 4 for the lower grid, and 2 for the bipolar grid
  int gridid[ngrids], coordIds[ngrids][ndims], mosaicid, globalId;

  nc_type nc_mode = NC_CLOBBER;

  char *tname = "3tile_tripolar";
  char *name, filename[STRING_SIZE];
  name = ( char * )calloc( STRING_SIZE, sizeof( char ) );
  sprintf( name, "%s", tname );
  const char *dimnames[] = {"ni", "nj"};
  const char *bpdimnames[] = {"bp_ni", "bp_nj"};
  
  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* Get the bipolar cap 1 grid */
  status = nccf_get_bipolar_cap(bpdims, latPerim, lonSing, bplons, bplats);

  /* A two grids below the cap */
  for( k = 0; k < 2; k++ ){
    dLon = (lonMax[k] - lonMin[k]) / (gldims[0] - 1);
    dLat = (latMax    - latMin) / (gldims[1] - 1);
    for (j = 0; j < gldims[1]; ++j) {
      for (i = 0; i < gldims[0]; ++i) {
        gllons[i + gldims[0] * j] = lonMin[k] + i * dLon;
        gllats[i + gldims[0] * j] = latMin    + j * dLat;
      }
    }
    if ((status = 
      nccf_def_lon_coord(ndims, gldims, dimnames, gllons, save, &coordIds[k][0]))) ERR;
    if ((status = 
      nccf_def_lat_coord(ndims, gldims, dimnames, gllats, save, &coordIds[k][1]))) ERR;

/* Write the grid */
    sprintf( filename, "%s_grid%d.nc", name, k );
    if(( status = nccf_def_grid( coordIds[k], name, &gridid[k] )));
    if(( status = nccf_def_global( &globalId ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;
    if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
    if(( status = nccf_put_grid(gridid[k], ncid))) ERR;
    if(( status = nccf_put_global(globalId, ncid)));
    if(( status = nc_close( ncid ))) ERR;
    if(( status = nccf_free_global( globalId ))) ERR;

  }

  if ((status = 
    nccf_def_lon_coord(ndims, bpdims, bpdimnames, bplons, save, &coordIds[2][0]))) ERR;
  if ((status = 
    nccf_def_lat_coord(ndims, bpdims, bpdimnames, bplats, save, &coordIds[2][1]))) ERR;

/* Write the mosaics */
  sprintf( filename, "%s_grid%d.nc", name, 2 );
  status = nccf_def_grid( coordIds[2], name, &gridid[k] ); //bipolar cap grid
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_grid(gridid[k], ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

  if(( status = nccf_def_mosaic( ngrids, gridid, name, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;
  sprintf( filename, "%s_mosaics.nc", name );
  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid)));
  if(( status = nc_close( ncid ))) ERR;

/* write mosaic as polytopes */
//  nccf_print_mosaic_as_polytopes(mosaicid,"tripolar.pt");

  if(( status = nccf_free_global( globalId ))) ERR;
  if(( status = nccf_free_mosaic( mosaicid ))) ERR;
  free(name);
 
  for (k=0; k<ngrids; ++k) {
    nccf_free_grid(gridid[k]);
    for (j=0; j<ndims; ++j)
        if( (status = nccf_free_coord(coordIds[k][j])) ) ERR;
  }
    
  return NC_NOERR;
}
