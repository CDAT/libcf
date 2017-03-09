/**
 * Test the creation of a cubed-sphere structured grid.
 *
 * "$Id: tst_cubesphere_mosaic.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_mosaic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>

#include "libcf_src.h"
#include "nccf_coord.h"
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"
#include "nccf_handle_error.h"
#include "nccf_utility_functions.h"

void fill_func(double *in,double *out){
   out[0]=in[0]+in[1];
}
int main(){

/* to save the data in the coordinate object */
  const int save = 1;

/* Variables for the grid */
  int ncid, mosaicid, mosaicid2;
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

  double clon[nvertex], clat[nvertex], data_values[nvertex];
  double period[] = {0., 0.};

  char *name = "tst_cubesphere";
  char gridname[STRING_SIZE], filename[STRING_SIZE];
  const char *dimnames[] = {"ni", "nj"};

  int coordIds[ngrids][ndims], globalId, gridids[ngrids];

/* Generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;  


/* Define a grid */
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    pos = iGrid % 3;
    sign = ( iGrid/3 )  == 0 ? 1 : -1;
    faceVec[pos] = sign;
    
    nccf_get_cubedsphere_grid( dims, faceVec, clon, clat );
    faceVec[pos] = 0; // reset to zero.

    // set data
    for (i = 0; i < nvertex; ++i) {
      data_values[i] = sin(clat[i]*M_PI/180.0) * cos(clon[i]*2.0*M_PI/180.0);
    }

    if(( status = nccf_def_lon_coord( ndims, dims, dimnames, clon, save,
                                            &coordIds[iGrid][0] ))) ERR;
    if(( status = nccf_def_lat_coord( ndims, dims, dimnames, clat, save,
                                            &coordIds[iGrid][1] ))) ERR;

/* Define the grid for the current tile. */
    sprintf( gridname, "%s%s%d", name, CF_FILENAME_GRID, iGrid );
    if(( status = nccf_def_grid( coordIds[iGrid], 
					    gridname, 
					    &gridids[iGrid] ))) ERR;
    if(( status = nccf_def_global( &globalId ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_GRIDNAME, name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalId, CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

    sprintf( filename, "%s.nc", gridname );
    if(( status = nc_create(filename, nc_mode, &ncid) )) ERR;
    if(( status = nccf_put_grid(gridids[iGrid], ncid) )) ERR;
    int data_id;
    const char *timeDim = NULL;
    if(( status = nccf_def_data(gridids[iGrid], 
					   "p", "pressure", "Pa", 
					   timeDim, &data_id) )) ERR;
    if(( status = nccf_add_global_att( globalId, CF_FILETYPE,   
                                       CF_GLATT_FILETYPE_STATIC_DATA, 1 ))) ERR;
    if(( status = nccf_set_data_double(data_id, data_values, 
					      save, NC_FILL_DOUBLE) )) ERR;
    if(( status = nccf_put_data(data_id, ncid) )) ERR;
    if(( status = nccf_put_global(globalId, ncid)));
    if(( status = nc_close( ncid )));
    if(( status = nccf_free_global( globalId ))) ERR;

  }

  sprintf( filename, "%s%s.nc", name, CF_FILENAME_MOSAIC );
  if(( status = nccf_def_mosaic( ngrids, gridids, gridname, &mosaicid ))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, period ))) ERR;

  if(( status = nc_create( filename, nc_mode, &ncid ))) ERR;
  if(( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nc_close( ncid ))) ERR;

//  if(( status = nccf_print_mosaic_as_polytopes( mosaicid, "cubesphere.pt" ))) ERR;

  if(( status = nccf_free_mosaic( mosaicid ))) ERR;
  
  /* free the contacts object */
  /* free the grid and the coordinates */
  int coordids[ngrids];
  for( iGrid = 0; iGrid < ngrids; iGrid++ ){
    if(( status = nccf_inq_grid_coordids( gridids[iGrid], coordids ))) ERR;
    if(( status = nccf_free_grid( gridids[iGrid] ))) ERR;
    for (i = 0; i < ndims; ++i) {
      if(( status = nccf_free_coord( coordids[i] ))) ERR;
    }
    if(( status = nccf_free_data(gridids[iGrid]) )) ERR;
  }

  /* Load the contacts, grids, and coordinates from file */
  if ((status = nccf_def_mosaic_from_file( filename,
					    "tst_two_tiles_mosaic",
					    &mosaicid2))) ERR;

  /* free the grid and the coordinates */
  status = nccf_free_mosaic( mosaicid2 );

  /* Close the file*/
  return 0;
}
