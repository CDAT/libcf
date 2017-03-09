/**
 * Test interpolation on a periodic grid with the source grid undoing a 
 * periodicity jump somewhere in the middle
 *
 * $Id: tst_periodic2.c 904 2011-12-28 21:58:46Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

#include "nccf_regrid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <netcdf.h>
#include <libcf_src.h>

#include <nccf_coord.h>
#include <nccf_grid.h>
#include <nccf_global.h>
#include <nccf_data.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

#include <cf_config.h>

void setValues(int data_id) {
  int lonid, gridid, status, i, n;
  double *clon, *data;
  const double *fillval;
  nc_type xtype;
  if ((status = nccf_inq_data_gridid(data_id, &gridid))) ERR;
  if ((status = nccf_inq_grid_coordids(gridid, &lonid))) ERR;
  if ((status = nccf_get_coord_data_pointer(lonid, &clon))) ERR;
  if ((status = nccf_inq_coord_dims(lonid, &n))) ERR;
  if ((status = nccf_get_data_pointer(data_id, &xtype, (void **)&data, 
                                      (const void **)&fillval))) ERR;
  for (i = 0; i < n; ++i) {
    data[i] = sin(M_PI*clon[i]/180.0);
  }
}

void createLon1D(int n, double xmin, int *lon_id, int *grid_id, int *data_id) {
  const int ndims = 1;
  const int save = 1;
  int status;
  int i;
  double *clon;
  double *data;
  double dx;

  const char *dimnames[] = {"ni"};

  clon = (double *) malloc(n * sizeof(double));
  data = (double *) malloc(n * sizeof(double));

  dx = 360.0 / ( (double) (n - 1) );
  for (i = 0; i < n; ++i) {
    clon[i] = xmin + i * dx;
  }

  if ((status = nccf_def_lon_coord(ndims, &n, dimnames, clon, save, 
                                   lon_id))) ERR;
  if ((status = nccf_def_grid(lon_id, "lon_grid", grid_id))) ERR;
  if ((status = nccf_def_data(*grid_id, "lon_data", 
                              NULL, NULL, NULL, data_id))) ERR;

  /* Initialize the data */
  for (i = 0; i < n; ++i) {
    data[i] = NC_FILL_DOUBLE;
  }
  if ((status = nccf_set_data_double(*data_id, data, save,
					    NC_FILL_DOUBLE))) ERR;

  /* Set the data */
  setValues(*data_id);

  free(clon);
  free(data);
}

//////////////////////////////////////////////////////////////////////

double 
checkInterp(int dataId1, int dataId2) {
  int status, k, ndims, ntot, i;
  double err;
  int gridId;
  if ((status = nccf_inq_data_gridid(dataId1, &gridId))) ERR;
  if ((status = nccf_inq_grid_ndims(gridId, &ndims))) ERR;
  int coordIds[ndims];
  if ((status = nccf_inq_grid_coordids(gridId, coordIds))) ERR;
  int dims[ndims];
  if ((status = nccf_inq_coord_dims(coordIds[0], dims))) ERR;
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }
  
  /* compare the two sets of data */
  double *data1;
  double *data2;
  nc_type xtype;
  const void *fill_value;
  if ((status = nccf_get_data_pointer(dataId1, &xtype, (void **) &data1,
					     &fill_value))) ERR;
  if ((status = nccf_get_data_pointer(dataId2, &xtype, (void **) &data2,
					     &fill_value))) ERR;

  err = 0.0;
  for (k = 0; k < ntot; ++k) {
    err += fabs(data1[k] - data2[k]);
  }
  err /= ntot;

#ifdef LOGGING
    printf("Average interpolation error: %lf\n", err);
#endif /* LOGGING */

  return err;
}

//////////////////////////////////////////////////////////////////////

int main(){

  int status, i;
  int ori_lon_id, ori_grid_id, ori_data_id; 
  int tgt_lon_id, tgt_grid_id, tgt_data_id; 
  int tgt_data_ref_id;
  double *tgt_data_ref;
  const int save = 1;

  /* Create grids */
  const int ori_n = 11;
  const int tgt_n = 11;

  // original nodal grid
  createLon1D(ori_n, -180.0, &ori_lon_id, &ori_grid_id, &ori_data_id);

  // target nodal grid
  createLon1D(tgt_n, 0.0, &tgt_lon_id, &tgt_grid_id, &tgt_data_id);

  // reference data on target grid
  if ((status = nccf_def_data(tgt_grid_id, "data", 
                              NULL, NULL, NULL, &tgt_data_ref_id))) ERR;
  tgt_data_ref = (double *) malloc(tgt_n * sizeof(double));
  for (i = 0; i < tgt_n; ++i) {
    tgt_data_ref[i] = NC_FILL_DOUBLE;
  }
  if ((status = nccf_set_data_double(tgt_data_ref_id, 
                                     tgt_data_ref, save,
                                     NC_FILL_DOUBLE))) ERR;
  setValues(tgt_data_ref_id);       
  
#ifdef HAVE_LAPACK_LIB
  /* Create regrid object */
  int regrid_id;
  const int nitermax = 4;
  const double tolpos = 1.e-3;
  if ((status = nccf_def_regrid(ori_grid_id, tgt_grid_id, 
                                &regrid_id))) ERR;

  if ((status = nccf_compute_regrid_weights(regrid_id, nitermax, 
					    tolpos))) ERR;

  /* Interpolate */
  if ((status = nccf_apply_regrid(regrid_id, ori_data_id,
                                  tgt_data_id))) ERR;

  /* Check */
  int nvalid, ntargets;
  double ratio;
  if ((status = nccf_inq_regrid_ntargets(regrid_id, &ntargets))) ERR;
  if ((status = nccf_inq_regrid_nvalid(regrid_id, &nvalid))) ERR;
  ratio = (double)(nvalid) / (double)(ntargets);
#ifdef LOGGING
  printf("ratio of valid to num target points (nodal) = %f\n", ratio);
#endif /* LOGGING */
  assert(ratio == 1.0);

  checkInterp(tgt_data_id, tgt_data_ref_id);

  /* Clean up */
  if ((status = nccf_free_regrid(regrid_id))) ERR;
#endif

  if ((status = nccf_free_data(ori_data_id))) ERR;
  if ((status = nccf_free_data(tgt_data_id))) ERR;
  if ((status = nccf_free_data(tgt_data_ref_id))) ERR;
    
  if ((status = nccf_free_grid(ori_grid_id))) ERR;
  if ((status = nccf_free_grid(tgt_grid_id))) ERR;

  if ((status = nccf_free_coord(ori_lon_id))) ERR;
  if ((status = nccf_free_coord(tgt_lon_id))) ERR;

  free(tgt_data_ref);

  return 0;
}
