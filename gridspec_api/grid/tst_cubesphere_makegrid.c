/**
 * Test the definition of a cube-sphere structured grid.
 *
 * "$Id: tst_cubesphere_makegrid.c 767 2011-06-06 23:20:19Z pletzer $"
 */



#include "nccf_grid.h"
#include <stdio.h>
#include <math.h>
#include <netcdf.h>

#include "libcf_src.h"
#include "nccf_global.h"
#include "nccf_coord.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

int main(){

  /* to save the data in the coordinate object */
  const int save = 1;

  /* Variables for the grid */
  int ncid;
  int status;
  int iGrid; // Some counters

  int ndims = 2;         // Number of x, y dims
  int ngrids = 6;
  int nlon = 10;         // number of x faces, gives nlon + 1 edges
  int elon = nlon + 1;   // # of edges
  int dims[] = {elon, elon};
  int nvertex = elon*elon;
  int faceVec[3] = {0, 0, 0}, pos, sign;
  char coordinates_id[36+1];
  int i;

  nc_type nc_mode = NC_CLOBBER;

  double clon[nvertex], clat[nvertex];

  char *name = "cubesphere_tile";
  char gridname[STRING_SIZE], filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  int coordids[ngrids][ndims];
  int gridid[ngrids];  // Each Id goes with a ncid
  int globalId;

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
                                            &coordids[iGrid][0] ))) ERR;
    if(( status = nccf_def_lat_coord( ndims, dims, dimnames, clat, save,
                                           &coordids[iGrid][1] ))) ERR;
    /* Define the grid for the current tile. */
    sprintf( gridname, "%s%d", name, iGrid );
    if(( status = nccf_def_grid( coordids[iGrid], 
					    gridname, &gridid[iGrid] ))) ERR;
    if ((status = nccf_def_global( &globalId ))) ERR;
    if ((status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                       coordinates_id, 0 ))) ERR;

    /* Put the grid into netcdf*/
    sprintf( filename, "%s%d.nc", name, iGrid );
    if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
    if(( status = nccf_put_grid(gridid[iGrid], ncid)));
    if(( status = nccf_put_global(globalId, ncid))) ERR;
    if(( status = nc_close( ncid ))) ERR;

    /* Free the global grid */
    if(( status = nccf_free_global( globalId ))) ERR;

  }

  /* free the grid */
  for( i = 0; i < ngrids; i++ ){
    if(( status = nccf_free_grid( gridid[i] )));
    if(( status = nccf_free_coord( coordids[i][0] ))) ERR;
    if(( status = nccf_free_coord( coordids[i][1] ))) ERR;
  }

  /* Close the file*/
  return NC_NOERR;
}
