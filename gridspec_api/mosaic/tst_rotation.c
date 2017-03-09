/**
 * Two tile mosaic with contact as a rotation in index space
 *
 *
 * "$Id: tst_rotation.c 767 2011-06-06 23:20:19Z pletzer $"
 */

// libcf/gridspec includes
#include "nccf_mosaic.h"
#include "nccf_global.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

// standard includes
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

int main() {

  const int ngrids = 2;
  const int ndims = 2;
  char coordinates_id[36+1];

  // first tile
  const int dims1[] = {3 + 1, 2 + 1};
  const double xyMins1[] = {0., 0.};
  const double xyMaxs1[] = {1., 1.};

  // second tile
  const int dims2[] = {4 + 1, 3 + 1};
  const double xyMins2[] = {2., 0.};
  const double xyMaxs2[] = {1., 1.};

  assert(dims1[0] == dims2[1]); // because of rotation

  // compute coordinates
  int i, j;
  double *x1s, *y1s;
  double *x2s, *y2s;

  int n1 = 1, n2 = 1;
  for (i = 0; i < ndims; ++i) {
    n1 *= dims1[i];
    n2 *= dims2[i];
  }
  x1s = (double *) malloc(n1 * sizeof(double));
  y1s = (double *) malloc(n1 * sizeof(double));
  x2s = (double *) malloc(n2 * sizeof(double));
  y2s = (double *) malloc(n2 * sizeof(double));

  double dx1 = (xyMaxs1[0] - xyMins1[0]) / (dims1[1] - 1);
  double dy1 = (xyMaxs1[1] - xyMins1[1]) / (dims1[0] - 1);
  for(j = 0; j < dims1[0]; j++){
    for(i = 0; i < dims1[1]; i++){
      x1s[i + dims1[1]*j] = xyMins1[0] + i * dx1;
      y1s[i + dims1[1]*j] = xyMins1[1] + j * dy1;
    }
  }

  double dx2 = (xyMaxs2[0] - xyMins2[0]) / (dims2[1] - 1);
  double dy2 = (xyMaxs2[1] - xyMins2[1]) / (dims2[0] - 1);
  for(j = 0; j < dims2[0]; j++){
    for(i = 0; i < dims2[1]; i++){
      x2s[i + dims2[1]*j] = xyMins2[0] + i * dx2;
      y2s[i + dims2[1]*j] = xyMins2[1] + j * dy2;
    }
  }

  // create grid objects
  int coordIds[ngrids][ndims];
  const char *coordnames[] = {"x", "y"};
  const char *coordDimsNames[] = {"ny", "nx"};
  const int save = 0;
  int status;

  /* generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  if(( status = nccf_def_coord(ndims, dims1, coordDimsNames,
				    x1s, save,
				    coordnames[0], NULL, NULL,
				    &coordIds[0][0] ))) ERR;
  if(( status = nccf_def_coord(ndims, dims1, coordDimsNames,
				    y1s, save,
				    coordnames[1], NULL, NULL,
				    &coordIds[0][1] ))) ERR;
  if(( status = nccf_def_coord(ndims, dims2, coordDimsNames,
				    x2s, save,
				    coordnames[0], NULL, NULL,
				    &coordIds[1][0] ))) ERR;
  if(( status = nccf_def_coord(ndims, dims2, coordDimsNames,
				    y2s, save,
				    coordnames[1], NULL, NULL,
				    &coordIds[1][1] ))) ERR;

  // create grids
  int gridids[ngrids], globalIds[ngrids];
  char *gridNames[] = {"left", "right"};
  for (i = 0; i < ngrids; ++i) {
    if (( status = nccf_def_grid(coordIds[i],
					    gridNames[i],
					    &gridids[i] ))) ERR;
    if(( status = nccf_def_global( &globalIds[i] ))) ERR;
    if(( status = nccf_add_global_att( globalIds[i], CF_FILETYPE,
                                       CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalIds[i], CF_GRIDNAME,
                                       gridNames[i], 0 ))) ERR;
    if(( status = nccf_add_global_att( globalIds[i], CF_COORDINATES_ID,
                                       coordinates_id, 0 ))) ERR;
  }

  // write grid files
  const nc_type nc_mode = NC_CLOBBER;
  char filename[STRING_SIZE];
  int ncid;
  for (i = 0; i < ngrids; ++i) {
    sprintf(filename, "%s_grid%d.nc", "tst_rotation", i);
    nc_create(filename, nc_mode, &ncid);
    nccf_put_grid(gridids[i], ncid);
    nccf_put_global(globalIds[i], ncid);
    nc_close(ncid);
  }

  // define the mosaic, do this after writing the grid files to disk
  int mosaicid, globalId;
  const double periods[] = {0.0, 0.0};
  if ((status = nccf_def_mosaic(ngrids, gridids, "tst_rotation", &mosaicid))) ERR;
  if(( status = nccf_compute_mosaic_contacts( mosaicid, periods ))) ERR;
  if(( status = nccf_def_global( &globalId ))) ERR;
  if(( status = nccf_add_global_att( globalId, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

  // write mosaic
  sprintf(filename, "%s_mosaic.nc", "tst_rotation");
  if(( status = nc_create(filename, nc_mode, &ncid))) ERR;
  if((status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if(( status = nccf_put_global(globalId, ncid))) ERR;
  if(( status = nc_close(ncid))) ERR;
  if(( status = nccf_free_global( globalId ))) ERR;

  // clean up
  if (( status = nccf_free_mosaic(mosaicid))) ERR;
  for (i = 0; i < ngrids; ++i) {
    if ((status ==  nccf_free_global( globalIds[i]))) ERR;
    if ((status ==  nccf_free_grid(gridids[i]))) ERR;
    for (j = 0; j < ndims; ++j) {
      if ((status ==  nccf_free_coord(coordIds[i][j]) )) ERR;
    }
  }
  free(x1s);
  free(y1s);
  free(x2s);
  free(y2s);

  return 0;
}
