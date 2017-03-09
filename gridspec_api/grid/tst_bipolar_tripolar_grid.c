/**
 * Test the creation of a tripolar grid
 * \author Dave Kinding and Lexander Pletzer, Tech-X Corp.
 *
 * "$Id: tst_bipolar_tripolar_grid.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libcf_src.h>
#include <netcdf.h>
#include "nccf_global.h"
#include <nccf_coord.h>
#include <nccf_handle_error.h>
#include <nccf_utility_functions.h>

int main(){

  const int ndims = 2;
  const int ngrids = 3;
  const int bpdims[] = {6, 6}, gldims[] = {6, 6};
  const int bpnvertex = bpdims[0] * bpdims[1];
  const int glnvertex = gldims[0] * gldims[1];
  const double latPerim = 60;
  const double lonSing = 0.0;
  double bplons[bpnvertex], bplats[bpnvertex];
  double gllons[glnvertex], gllats[glnvertex];
  double dLon, dLat, lonMin[] = {0, 180}, lonMax[] = {180, 360};
  double latMin = -90, latMax = latPerim;
  int i, j, iGrid, status, ncid;
  const int save = 1;
  char coordinates_id[36+1];

  int gridid[ngrids], coordids[ngrids][ndims];  // 3 for the lower grid, and 2 for the bipolar grid
  int globalId;

  nc_type nc_mode = NC_CLOBBER;

  const char *dimnames[] = {"ni", "nj"};
  char *tname = "3tile_tripolar";
  char *name, filename[STRING_SIZE];
  name = ( char * )calloc( STRING_SIZE, sizeof( char ) );
  sprintf( name, "%s", tname );

  /* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;
  
  /* Get the bipolar cap 1 grid */
  status = nccf_get_bipolar_cap(bpdims, latPerim, lonSing, bplons, bplats);

  /* A two grids below the cap */
  for( iGrid = 0; iGrid < ngrids-1; iGrid++ ){
    dLon = (lonMax[iGrid] - lonMin[iGrid]) / (gldims[0] - 1);
    dLat = (latMax    - latMin) / (gldims[1] - 1);
    for (j = 0; j < gldims[1]; ++j) {
      for (i = 0; i < gldims[0]; ++i) {
        gllons[i + gldims[0] * j] = lonMin[iGrid] + i * dLon;
        gllats[i + gldims[0] * j] = latMin    + j * dLat;
      }
    }
    if ((status = nccf_def_lon_coord(ndims, gldims, dimnames, 
                    gllons, save, &coordids[iGrid][0]))) ERR;
    if ((status = nccf_def_lat_coord(ndims, gldims, dimnames, 
                    gllats, save, &coordids[iGrid][1]))) ERR;

    /* Write the grid */
    sprintf( filename, "%s_grid%d.nc", name, iGrid );
    status = nccf_def_grid( coordids[iGrid], name, &gridid[iGrid] );
    if (status) ERR;
    if ((status = nccf_def_global( &globalId ))) ERR;
    if ((status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                       coordinates_id, 0 ))) ERR;
    status = nc_create( filename, nc_mode, &ncid );
    if (status) ERR;
    status = nccf_put_grid(gridid[iGrid], ncid);
    status = nccf_put_global(globalId, ncid);
    if (status) ERR;
    status = nc_close( ncid );
    if (status) ERR;
    if((status = nccf_free_global( globalId ))) ERR;
  }

  if ((status = nccf_def_lon_coord(ndims, bpdims, dimnames, bplons, 
                save, &coordids[iGrid][0]))) ERR;
  if ((status = nccf_def_lat_coord(ndims, bpdims, dimnames, bplats, 
                save, &coordids[iGrid][1]))) ERR;

  /* Write the polar cap grid */
  sprintf( filename, "%s_grid%d.nc", name, 2 );
  status = nccf_def_grid( coordids[iGrid], name, &gridid[iGrid] );
  if ((status = nccf_def_global( &globalId ))) ERR;
  if ((status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                     coordinates_id, 0 ))) ERR;
  nc_create( filename, nc_mode, &ncid );
  nccf_put_grid(gridid[iGrid], ncid);
  status = nccf_put_global(globalId, ncid);
  nc_close( ncid );

  /* Free all of the grids */
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    status = nccf_inq_grid_coordids( gridid[iGrid], coordids[iGrid] );
    status += nccf_free_grid( gridid[iGrid] );
    status += nccf_free_coord( coordids[iGrid][0] );
    status += nccf_free_coord( coordids[iGrid][1] );
  }
  nccf_free_global( globalId );

  free( name );

  return NC_NOERR;
}
