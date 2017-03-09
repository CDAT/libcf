/**
 * Test masked interpolation 
 * $Id: tst_mask.c 881 2011-12-17 21:53:14Z pletzer $
 */

#include <nccf_regrid.h>
#include <math.h>
#include <stdlib.h>
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
void create_unigrid(int ndims, const int dims[], 
		    const double lo_corner[], const double hi_corner[],
		    int coordids[], int *gridid) {
  int status;
  int i, j;
  int indx[ndims];
  const char *dimnames[] = {"nx", "ny", "nz"};
  const char *coordnames[] = {"x", "y", "z"};
  int ntot = get_ntot(ndims, dims);
  double **coords;
  coords = (double **) malloc(ndims * sizeof(double *));
  double deltas[ndims];
  for (i = 0; i < ndims; ++i) {
    coords[i] = (double *) malloc(ntot * sizeof(double));
    deltas[i] = (hi_corner[i]-lo_corner[i])/(double)(dims[i]-1);
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
  status = nccf_def_grid (coordids, "xyz", gridid);
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
  const double fill_value = -1;
  const int ndims = 3;
  const double lo_corner[] = {0, 0, 0};
  const double hi_corner[] = {1, 1, 1};
  int indx[ndims];

  // target grid
  const int dims_tgt[] = {50, 51, 52};
  int coordids_tgt[ndims];
  int gridid_tgt;
  create_unigrid(ndims, dims_tgt, lo_corner, hi_corner,
		 coordids_tgt, &gridid_tgt);
  status = nccf_save_grid_scrip(gridid_tgt, "tst_mask_tgt_scrip.nc");
  if (status) ERR;

  // put data on the target grid
  int ntot_tgt = get_ntot(ndims, dims_tgt);
  int dataid_tgt;
  double *data_tgt = (double *) malloc(ntot_tgt * sizeof(double));
  for (j = 0; j < ntot_tgt; ++j) {
    // set data to bad value
    data_tgt[j] = fill_value;
  }
  status = nccf_def_data (gridid_tgt, "data_tgt", NULL, NULL, NULL, 
				 &dataid_tgt);
  if (status) ERR;
  status = nccf_set_data_double(dataid_tgt, data_tgt, 
				       SAVE, fill_value);
  if (status) ERR;

  // original grid
  const int dims_ori[] = {100, 101, 102};
  int coordids_ori[ndims];
  int gridid_ori;
  create_unigrid(ndims, dims_ori, lo_corner, hi_corner,
		 coordids_ori, &gridid_ori);
  status = nccf_save_grid_scrip(gridid_ori, "tst_mask_ori_scrip.nc");
  if (status) ERR;

  // put data on the original grid
  int dataid_ori;
  int ntot_ori = get_ntot(ndims, dims_ori);
  int *imask = (int *) malloc(ntot_ori * sizeof(int));
  double *data_ori = (double *) malloc(ntot_ori * sizeof(double));
  for (j = 0; j < ntot_ori; ++j) {
    nccf_get_multi_index(ndims, dims_ori, j, indx);
    data_ori[j] = 0;
    for (i = 0; i < ndims; ++i) {
      double xyz = indx[i]/(double)(dims_ori[i] - 1);
      data_ori[j] += xyz*xyz;
      // only data_ori[j] < 1.0 is valid
      imask[j] = (data_ori[j] < 1.0? 1: 0);
    }
    data_ori[j] = sqrt(data_ori[j]);

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
  const int nitermax = 100;
  double tolpos = 1.e-4;
  status = nccf_compute_regrid_weights(regridid, nitermax, tolpos);
  if (status) ERR;
  status = nccf_apply_regrid(regridid, dataid_ori, dataid_tgt);
  if (status) ERR;
  status = nccf_free_regrid(regridid);
  if (status) ERR;
#endif

  // save
  int ncid;
  status = nc_create("tst_mask_tgt.nc", NC_CLOBBER, &ncid);
  if (status) ERR;
  status = nccf_put_grid(gridid_tgt, ncid);
  if (status) ERR;
  status = nccf_put_data(dataid_tgt, ncid);
  if (status) ERR;
  status = nc_close(ncid);
  if (status) ERR;
   status = nc_create("tst_mask_ori.nc", NC_CLOBBER, &ncid);
  if (status) ERR;
  status = nccf_put_grid(gridid_ori, ncid);
  if (status) ERR;
  status = nccf_put_data(dataid_ori, ncid);
  if (status) ERR;
  status = nc_close(ncid);
  if (status) ERR;

  // clean up
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
