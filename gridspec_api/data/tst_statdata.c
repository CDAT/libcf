/**
 * Unit test exercising static data creation
 *
 * $Id: tst_statdata.c 767 2011-06-06 23:20:19Z pletzer $
 */

#include "nccf_data.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "nccf_utility_functions.h"
#include "nccf_coord.h"
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_handle_error.h"

int main(int argc, char *argv[]) {
  
  int status;
  const int ndims = 1;
  int coordIds[ndims];
  int globalId;
  const int save = 0;
  int i;
  char coordinates_id[36+1];
  char data_id[36+1];

  printf("Running %s\n", argv[0]);

  /* generate unique ids */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;
  if ((status = nccf_generate_id(456, data_id))) ERR;

  /* create coordinates */
  const int nx1 = 11;
  const double dx = 1.0;
  double *xs = (double *) malloc(nx1 * sizeof(double));
  float *ps = (float *) malloc(nx1 * sizeof(float));
  for (i = 0; i < nx1; ++i) {
    xs[i] = dx * i; 
    ps[i] = 0.0; //(float)(xs[nx1-1] - xs[i])/(float)(xs[nx1-1] - xs[i]); 
  }
  int dims[] = {nx1};
  const char *dimNames[] = {"nx1"};
  status = nccf_def_coord(ndims, dims, dimNames, xs, save,
			       "x", "distance", "m", &coordIds[0]);
  if (status) ERR;

  /* create grid */
  int gridid;
  status = nccf_def_grid(coordIds, "x_grid", &gridid);
  status = nccf_def_global( &globalId );
  status = nccf_add_global_att( globalId, CF_FILETYPE, 
                                CF_GLATT_FILETYPE_GRID, 0 );

  /* create data */
  int pId;
  status = nccf_def_data(gridid, "p", "air_pressure", "Pa", NULL, &pId);
  if (status) ERR;
  status = nccf_set_data_float(pId, ps, save, NC_FILL_FLOAT);
  if (status) ERR;
  
  /* Add global attributes */
  status = nccf_add_global_att( globalId, CF_FILETYPE, 
                                CF_GLATT_FILETYPE_GRID, 1 );
  status = nccf_add_global_att( globalId, CF_FILETYPE, 
                                CF_GLATT_FILETYPE_STATIC_DATA, 0 );
  status = nccf_add_global_att( globalId, CF_COORDINATES_ID, 
                                coordinates_id, 0 );
  status = nccf_add_global_att( globalId, CF_DATA_ID, 
                                data_id, 0 );

  /* write everything to disk */
  int ncid;
  if ((status = nc_create("tst_statdata.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_grid(gridid, ncid))) ERR;
  if ((status = nccf_put_data(pId, ncid))) ERR;
  if ((status = nccf_put_global(globalId, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  /* clean up */
  if ((status = nccf_free_data(pId))) ERR;
  if ((status = nccf_free_grid(gridid))) ERR;
  if ((status = nccf_free_global(globalId))) ERR;
  if ((status = nccf_free_coord(coordIds[0]))) ERR;
  free(xs);
  free(ps);

  return 0;
}
