/**
 * Unit test for nccf_def_data_from_file
 *
 * "$Id: tst_nccf_def_SD_from_grid.c 851 2011-11-08 14:37:20Z pletzer $"
 */

#include "nccf_data.h"
#include <assert.h>
#include <stdio.h>
#include "nccf_grid.h"
#include "nccf_coord.h"
#include "nccf_constants.h"

int
main(){

  const char* filename = "tst_timedata.nc";
  const char* varname  = "p";
  int status = NC_NOERR, read_data = 1, i, ndims, ncid;
  int nGrids = 1;
  int dataid, *gridids, nDatadims, ntot;
  nc_type type;
  float *data;
  char *coordNames[] = {"x", "y", "z"};
  char tilename[STRING_SIZE];
  float *fv = NULL;

  /* Get the dimensions */
  status = nc_open( filename, NC_NOWRITE, &ncid );
  if ( status > 0 ) return status;
  if ( status ) ERR;
  if ((status = nc_inq_ndims( ncid, &ndims ))) ERR;
  if ((status = nc_close( ncid ))) ERR;

  /* Get the grids 
   * NOTE: ndims -1 is because we only have x, y, z as defined coordinates 
   *       time is a dimension and used in the data, but no values are given 
   *       for the time dimension */
  gridids = (int *) malloc(nGrids * sizeof(int));
  if ((status = nccf_def_grid_from_file( filename,
						ndims-1, (const char**)coordNames, tilename, &gridids[0]))) ERR;

  if ((status = nccf_def_data_from_file( filename, gridids[0], varname,
						read_data, &dataid))) ERR;

  /* Inquire about the dimensionality of the data */
  if ((status = nccf_inq_data_ndims(dataid, &nDatadims))) ERR;
  int dataDims[nDatadims];
  if ((status = nccf_inq_data_dims(dataid, dataDims))) ERR;

  ntot = 1;
  for( i = 0; i < nDatadims; i++ ) {
    ntot *= dataDims[i];
  }

  /* Get the data */
  if ((status = nccf_inq_data_type(dataid, &type))) ERR;
  nc_type xtype;
  if ((status = nccf_get_data_pointer(dataid, &xtype, 
					     (void **) &data, 
					     (const void **) &fv ))) ERR;

  /* Make some assertions */
  assert(dataid == 0);
  assert(type   == 5);
  assert(ntot   == 21483);
  assert(abs( *fv - NC_FILL_FLOAT ) < 1.e-7);

  /* This fails in valgrind for some reason */
  int ind = ntot-2;
  assert( data[ind] >= 1.99 );

  /* Clean up */
  int *coordIds;
  coordIds = ( int* )malloc((ndims-1) * sizeof(int));
  if ((status = nccf_inq_grid_coordids( gridids[0], coordIds))) ERR;
  for (i = 0; i < ndims-1; i++) {
    if ((status = nccf_free_coord(coordIds[i]))) ERR;
  }
  for (i = 0; i < nGrids; i++) {
    if ((status = nccf_free_grid(gridids[i]))) ERR;
  }
  if ((status = nccf_free_data( dataid ))) ERR;

  free(coordIds);
  free(gridids);

  return status;
}
