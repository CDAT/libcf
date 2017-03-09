/**
 * Test interpolation from one lon-lat grid to another
 *
 * $Id: tst_put_and_def_from_file.c 881 2011-12-17 21:53:14Z pletzer $
 *
 * \author Dave Kindig, Tech-X Corp.
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
#include <nccf_data.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

#include <cf_config.h>

void createLonLat(const double xymin[], const double xymax[], 
		  const int dims[],
		  int coordIds[], int *gridId, int *dataId){
  const int ndims = 2;
  int nvertex = dims[0]*dims[1];
  const int save = 1;
  int status;
  int i, j, k;
  double *clon, *clat;
  double *data;
  double dxs[ndims];

  const char *dimnames[] = {"nj", "ni"};

  clat = ( double* )malloc( nvertex * sizeof( double ));
  clon = ( double* )malloc( nvertex * sizeof( double ));
  data = ( double* )malloc( nvertex * sizeof( double ));

  for (i = 0; i < ndims; ++i) {
    dxs[i] = (xymax[i] - xymin[i]) / (dims[i] - 1);
  }

  /* Populate coordinates and create lon/lat coordinate objects */
  for (j = 0; j < dims[0]; ++j) {
    for (i = 0; i < dims[1]; ++i) {
      k = i + dims[1]*j;
      clat[k] = xymin[0] + j*dxs[0];
      clon[k] = xymin[1] + i*dxs[1];
    }
  }  

  if ((status = nccf_def_lat_coord(ndims, dims, dimnames, clat, save, &coordIds[0]))) ERR;
  if ((status = nccf_def_lon_coord(ndims, dims, dimnames, clon, save, &coordIds[1]))) ERR;
  if ((status = nccf_def_grid(coordIds, "tripolar_grid", gridId))) ERR;
  if ((status = nccf_def_data(*gridId, "data_lonlat", 
				     NULL, NULL, NULL, dataId))) ERR;

  /* Set the data */
  for (k = 0; k < nvertex; ++k) {
    data[k] = cos(M_PI*clon[k]/180.0) * sin(M_PI*clat[k]/180.0);
  }
  if ((status = nccf_set_data_double(*dataId, data, save, 
					    NC_FILL_FLOAT))) ERR;

  free(clat);
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

  const int ndims = 2;
  int status;
  int ori_grid_id, oriDataId, ori_coord_ids[ndims]; 
  int targetGridId, targetDataId, targetCoordIds[ndims];
  int targetGridIdRef, targetDataIdRef, targetCoordIdsRef[ndims];

  /* Create grids */
  const double oriXymin[] = {-90.0, 0.0};
  const double oriXymax[] = {+90.0, 360.0};
  const int oriDims[] = {60, 61};
  const double targetXymin[] = {-89.0, 0.0};
  const double targetXymax[] = {+89.0, 360.0};
  const int targetDims[] = {12, 13};
  int i;

  createLonLat(oriXymin, oriXymax, oriDims, ori_coord_ids, 
	       &ori_grid_id, &oriDataId);
  createLonLat(targetXymin, targetXymax, targetDims, targetCoordIdsRef, 
		   &targetGridIdRef, &targetDataIdRef);
  createLonLat(targetXymin, targetXymax, targetDims, targetCoordIds, 
		   &targetGridId, &targetDataId);
  
#ifdef HAVE_LAPACK_LIB
  int ncid;
  const char* put_filename = "tst_put_and_def_from_file_test.nc";
  int nnodes, nvalid, ntargets;
  int nts, nns, rgd_id;
  /* Create regrid object */
  int regrid_id;
  const int nitermax = 100;
  const double tolpos = 1.e-3;
  if ((status = nccf_def_regrid(ori_grid_id, targetGridId, &regrid_id))) ERR;
  if ((status = nccf_compute_regrid_weights(regrid_id, nitermax, 
					    tolpos))) ERR;

  /* Write the weights, indices and inside_domain vars */
  if (( status = nc_create( put_filename, NC_CLOBBER, &ncid ))) ERR;
  if (( status = nccf_put_regrid(regrid_id, ncid)));
  if (( status = nc_close( ncid ))) ERR; 

  /* Write the data and coordinate files */
  if ((status = nc_create("tst_lonlat_regrid_target_ref.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_grid(targetGridId, ncid))) ERR;
  if ((status = nccf_put_data(targetDataIdRef, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  /* Interpolate */
  if ((status = nccf_apply_regrid(regrid_id, oriDataId, targetDataId))) ERR;

  /* Write the data and coordinate files */
  if ((status = nc_create("tst_lonlat_regrid_target.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_grid(targetGridId, ncid))) ERR;
  if ((status = nccf_put_data(targetDataId, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  /* Check */
  if ((status = nccf_inq_regrid_ntargets(regrid_id, &ntargets))) ERR;
  if ((status = nccf_inq_regrid_nvalid(regrid_id, &nvalid))) ERR;
  if ((status = nccf_inq_regrid_nnodes( regrid_id, &nnodes))) ERR;
  double ratio = (double)(nvalid) / (double)(ntargets);
#ifdef LOGGING
  printf("ratio of valid to num target points = %f\n", ratio);
  printf( "\n%f\n", ratio );
#endif /* LOGGING */
  assert(ratio > 0.995);
  assert( checkInterp(targetDataId, targetDataIdRef) < 2.e-4 );

  /* Clean up */
  if ((status = nccf_free_regrid(regrid_id))) ERR;
#endif

  if ((status = nccf_free_data(oriDataId))) ERR;
  if ((status = nccf_free_data(targetDataId))) ERR;
  if ((status = nccf_free_data(targetDataIdRef))) ERR;
  if ((status = nccf_free_grid(ori_grid_id))) ERR;
  if ((status = nccf_free_grid(targetGridId))) ERR;
  if ((status = nccf_free_grid(targetGridIdRef))) ERR;
  for (i = 0; i < ndims; ++i) {
    if ((status = nccf_free_coord(ori_coord_ids[i]))) ERR;
    if ((status = nccf_free_coord(targetCoordIds[i]))) ERR;
    if ((status = nccf_free_coord(targetCoordIdsRef[i]))) ERR;
  }

#ifdef HAVE_LAPACK_LIB
  /* Retrieve the weights etc from file */
  if (( status = nccf_def_regrid_from_file( put_filename, &rgd_id ))) ERR;
  if (( status = nccf_inq_regrid_ntargets( rgd_id, &nts ))) ERR;
  if (( status = nccf_inq_regrid_nnodes( rgd_id, &nns ))) ERR;
  assert( ntargets == nts );
  assert( nnodes == nns );
  if (( status = nccf_free_regrid( rgd_id ))) ERR;
#endif


  return 0;
}
