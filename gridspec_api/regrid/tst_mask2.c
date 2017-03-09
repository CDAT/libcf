/**
 * Test masked interpolation 
 * $Id: tst_mask2.c 923 2012-03-22 18:44:09Z dkindig $
 */

#include <nccf_regrid.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <netcdf.h>
#include <nccf_utility_functions.h>
#include <nccf_coord.h>
#include <nccf_grid.h>
#include <nccf_data.h>
#include <nccf_handle_error.h>

#include <cf_config.h> // to access HAVE_LAPACK_LIB

const int SAVE = 1;

/** return total number of nodes */
int get_ntot(int ndims, const int dims[]) {
  int ntot = 1;
  int i;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }
  return ntot;
}

/** create a uniform grid */
void create_unigrid(const int dims[], 
		    const double lo_corner[], const double hi_corner[],
		    int coordids[], int *gridid) {
  const int ndims = 2;
  int status;
  int i, j;
  int indx[ndims];
  const char *dimnames[] = {"ny", "nx"};
  const char *coordnames[] = {"y", "x"};
  int ntot = get_ntot(ndims, dims);
  double **coords;
  coords = (double **) malloc(ndims * sizeof(double *));
  double deltas[ndims];
  for (i = 0; i < ndims; ++i) {
    coords[i] = (double *) malloc(ntot * sizeof(double));
    deltas[i] = (hi_corner[i] - lo_corner[i])/(double)(dims[i]-1);
  }
  for (j = 0; j < ntot; ++j) {
    nccf_get_multi_index(ndims, dims, j, indx);
    for (i = 0; i < ndims; ++i) {
      coords[i][j] = lo_corner[i] + indx[i]*deltas[i];    
    }
  }
  // create the coordinates
  for (i = 0; i < ndims; ++i) {
    status = nccf_def_coord(ndims, dims, dimnames, coords[i], SAVE, 
			    coordnames[i], NULL, NULL, &coordids[i]);
    if (status) ERR;
  }
  // create the grid
  status = nccf_def_grid(coordids, "yx", gridid);
  if (status) ERR;
  // clean up, all the data are copied
  for (i = 0; i < ndims; ++i) {
    free(coords[i]);
  }
  free(coords);
}

int main() {

  int i, j, status;

  // global settings
  const double fill_value = 0;
  const int ndims = 2;
  const double lo_corner[] = {0, 0};
  const double hi_corner[] = {1, 1};
  int indx[ndims];

  // target grid
  const int dims_tgt[] = {5, 11};
  int coordids_tgt[ndims];
  int gridid_tgt;
  create_unigrid(dims_tgt, lo_corner, hi_corner,
		 coordids_tgt, &gridid_tgt);

  // put data on the target grid
  int ntot_tgt = get_ntot(ndims, dims_tgt);
  int dataid_tgt;
  double *data_tgt = (double *) malloc(ntot_tgt * sizeof(double));
  for (j = 0; j < ntot_tgt; ++j) {
    // set data to bad value
    data_tgt[j] = fill_value;
  }
  status = nccf_def_data(gridid_tgt, "data_tgt", NULL, NULL, NULL, 
				 &dataid_tgt);
  if (status) ERR;
  status = nccf_set_data_double(dataid_tgt, data_tgt, 
				       SAVE, fill_value);
  if (status) ERR;

  // original grid
  const int dims_ori[] = {11, 21};
  int coordids_ori[ndims];
  int gridid_ori;
  create_unigrid(dims_ori, lo_corner, hi_corner,
		 coordids_ori, &gridid_ori);

  // put data on the original grid
  double **coords_ori;
  coords_ori = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    status = nccf_get_coord_data_pointer(coordids_ori[i], 
                                         (double **)&coords_ori[i]);
    if (status) ERR;
  }
  
  int dataid_ori;
  int ntot_ori = get_ntot(ndims, dims_ori);
  int *imask = (int *) malloc(ntot_ori * sizeof(int));
  double *data_ori = (double *) malloc(ntot_ori * sizeof(double));
  for (j = 0; j < ntot_ori; ++j) {
    nccf_get_multi_index(ndims, dims_ori, j, indx);
    double y = coords_ori[0][j];
    double x = coords_ori[1][j];
    // set mask
    imask[j] = (y >= x? 1: 0);
    data_ori[j] = 1.0;
  }
  status = nccf_def_data (gridid_ori, "data_ori", NULL, NULL, NULL, 
				 &dataid_ori);
  if (status) ERR;
  status = nccf_set_data_double (dataid_ori, data_ori, 
					SAVE, NC_FILL_DOUBLE);
  if (status) ERR;

  // put a mask on the original grid
  status = nccf_set_grid_validmask(gridid_ori, imask);
  if (status) ERR;

  // interpolate
#ifdef HAVE_LAPACK_LIB

  int regridid;
  status = nccf_def_regrid(gridid_ori, gridid_tgt, &regridid);
  if (status) ERR;
  const int nitermax = 20;
  double tolpos = 1.e-4;
  status = nccf_compute_regrid_weights(regridid, nitermax, tolpos);
  if (status) ERR;
  status = nccf_apply_regrid(regridid, dataid_ori, dataid_tgt);
  if (status) ERR;

  // check
  double **coords_tgt;
  coords_tgt = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    status = nccf_get_coord_data_pointer(coordids_tgt[i], 
                                         (double **)&coords_tgt[i]);
    if (status) ERR;
  }
  double *data_tgt_interp;
  nc_type xtype;
  const void *fv;
  status = nccf_get_data_pointer(dataid_tgt, &xtype, 
                                 (void **)&data_tgt_interp, &fv);
  if (status) ERR;
  for (j = 0; j < ntot_tgt; ++j) {
    nccf_get_multi_index(ndims, dims_tgt, j, indx);
    double y = coords_tgt[0][j];
    double x = coords_tgt[1][j];
    if (y >= x && data_tgt_interp[j] != 1) {
      printf("node = %d x = %lf y = %lf interp value = %lf (expected 1)\n", 
             j, x, y, data_tgt_interp[j]);
      assert(data_tgt_interp[j] == 1);
    }
    if (y < x && data_tgt_interp[j] != NC_FILL_DOUBLE) {
      printf("node = %d x = %lf y = %lf interp value = %lf (expected 0)\n", 
             j, x, y, data_tgt_interp[j]);
      assert(data_tgt_interp[j] == NC_FILL_DOUBLE);
    }
  }
  
  status = nccf_free_regrid(regridid);
  if (status) ERR;
  free(coords_tgt);
#endif

  // save
  int ncid;
  status = nc_create("tst_mask2_tgt.nc", NC_CLOBBER, &ncid);
  if (status) ERR;
  status = nccf_put_grid(gridid_tgt, ncid);
  if (status) ERR;
  status = nccf_put_data(dataid_tgt, ncid);
  if (status) ERR;
  status = nc_close(ncid);
  if (status) ERR;
   status = nc_create("tst_mask2_ori.nc", NC_CLOBBER, &ncid);
  if (status) ERR;
  status = nccf_put_grid(gridid_ori, ncid);
  if (status) ERR;
  status = nccf_put_data(dataid_ori, ncid);
  if (status) ERR;
  status = nc_close(ncid);
  if (status) ERR;

  // clean up
  free(coords_ori);
  free(data_tgt);
  free(data_ori);
  status = nccf_free_data(dataid_tgt);
  if (status) ERR;
  free(imask);
  status = nccf_free_data(dataid_ori);
  if (status) ERR;
  status = nccf_free_grid(gridid_tgt);
  if (status) ERR;
  status = nccf_free_grid(gridid_ori);
  if (status) ERR;
  for (i = 0; i < ndims; ++i) {
    status = nccf_free_coord(coordids_tgt[i]);
    if (status) ERR;
    status = nccf_free_coord(coordids_ori[i]);
    if (status) ERR;
  }
 
  return 0;
}
