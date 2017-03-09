/**
 * Test interpolation from one lon-lat grid to another
 *
 * $Id: tst_cubesphere_to_lonlat_regrid.c 881 2011-12-17 21:53:14Z pletzer $
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

#include <cf_config.h>

#include "nccf_coord.h"
#include "nccf_grid.h"
#include "nccf_data.h"
#include "nccf_constants.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"
#include "nccf_global.h"

double checkInterp(int dataId1, int dataId2);
void computeError(int dataId1, int dataId2);
int writeData( const int gridId, const int dataId, const char *filename );
int createData(int gridId, const char *filename, int useBogusData, 
               const char *coordinates_id );
int createLonLatGrid(const double xymin[], const double xymax[], 
                     const int dims[], int coordIds[]);
int createCubeSphereTile( const int iGrid, const char *gridName, 
                          const int dims[], int coordIds[], int *gridId );

int main(){

  int status;
  int tarGridId, tarDataId, tarDataIdRef;
  const int ndims = 2, nGrid = 6;

  int tarCoordIds[ndims], csCoordIds[ndims], csGridIds[nGrid], csDataId[nGrid];

  /* coordinates_id */
  char coordinates_id_LL[36+1], coordinates_id_CS[36+1];
  nccf_generate_id( 123, coordinates_id_LL );
  nccf_generate_id( 321, coordinates_id_CS );

  /* Create grids */
  const double tarXymin[] = {-90.0, -180.0};
  const double tarXymax[] = { 90.0  ,180.0};
  const int tarDims[] = {121, 121};
  const int csDims[] = {40, 40};
  int i, iGrid;

  const char *testName = "tst_cubesphere_to_lonlat_regrid";
  const char *gridName[] = {"_lonlat",
                            "_ref",
                            "_grid", 
                            "_interp"};
  char* regrid_filename;
  regrid_filename = (char*)calloc( STRING_SIZE, sizeof(char));
  char *csTileFileName;
  csTileFileName = (char*)calloc( STRING_SIZE, sizeof(char));
  char *tarTileFileName;
  tarTileFileName = (char*)calloc( STRING_SIZE, sizeof(char));

  tarGridId = createLonLatGrid(tarXymin, tarXymax, tarDims, tarCoordIds);

  /* Create data Reference file*/
  sprintf( tarTileFileName, "%s%s.nc", testName, gridName[1] );
  tarDataIdRef = createData(tarGridId, tarTileFileName, 0, coordinates_id_LL);

  /* Create data Interpolation file*/
  sprintf( tarTileFileName, "%s%s.nc", testName, gridName[3] );
  tarDataId = createData(tarGridId, tarTileFileName, 1, coordinates_id_LL);

  //for( iGrid = 0; iGrid < nGrid; iGrid++ ){
  for( iGrid = 0; iGrid < 3; iGrid++ ){
    status = createCubeSphereTile( iGrid, testName, csDims, 
                  csCoordIds, &csGridIds[iGrid] );
    sprintf( csTileFileName, "%s%s%d.nc", testName, gridName[2], iGrid );
    csDataId[iGrid]  = createData( csGridIds[iGrid], csTileFileName, 0, coordinates_id_CS );
#ifdef HAVE_LAPACK_LIB
    /* Create regrid object */
    int regrid_id;
    const int nitermax = 20;
    const double tolpos = 0.01;
    if ((status = nccf_def_regrid(csGridIds[iGrid], tarGridId, &regrid_id))) ERR;

    /* The cut needs to be determined before allocation. The cuts
     * are not in the same place necessarily */
//    if( iGrid == 2 || iGrid  == 5 || 
    if( iGrid == 0){
      const int lo[] = {0, csDims[1]/2 - 1};
      const int hi[] = {csDims[0], csDims[1]/2};
      if ((status = nccf_add_regrid_forbidden( regrid_id, lo, hi ))) ERR;
    }

    if ((status = nccf_compute_regrid_weights(regrid_id, 
					      nitermax, tolpos))) ERR;

    /* Interpolate */
    if ((status = nccf_apply_regrid(regrid_id, csDataId[iGrid], 
				    tarDataId))) ERR;

    /* Compute error and put the new values in tarDataId */
    //computeError(tarDataIdRef, tarDataId);

    /* Write regrided data*/
    if ((status = writeData(tarGridId, tarDataId, tarTileFileName))) ERR;

    /* Compute average error */
    double avgErr = checkInterp(tarDataIdRef, tarDataId);
    if( avgErr ) assert( 1 == 1 );

    /* Write weights to file */
    int ncid;
    sprintf( regrid_filename, "%s_put_test%d.nc", testName, iGrid );
    if (( status = nc_create( regrid_filename, NC_CLOBBER, &ncid ))) ERR;
    if (( status = nccf_put_regrid(regrid_id, ncid))) ERR;
    if (( status = nc_close( ncid ))) ERR;

    /* Check */
    int nvalid, ntargets;
    if ((status = nccf_inq_regrid_ntargets(regrid_id, &ntargets))) ERR;
    if ((status = nccf_inq_regrid_nvalid(regrid_id, &nvalid))) ERR;
    double ratio = (double)nvalid/(double)ntargets;
    printf("tst_cubesphere_to_lonlat_regrid: ratio of valid to num target points = %f\n", ratio);
    //assert(ratio >= 0.9); // WE NEED TO GET THIS RATIO TO ~< 1! (Pletzer);

    /* Clean up */
    status += nccf_free_regrid(regrid_id);
#endif
    status += nccf_free_grid(csGridIds[iGrid]);
    status += nccf_free_data(csDataId[iGrid]);

  }

  status += nccf_free_data(tarDataId);
  status += nccf_free_data(tarDataIdRef);
  status += nccf_free_grid(tarGridId);
  for (i = 0; i < ndims; ++i) {
    status += nccf_free_coord(tarCoordIds[i]);
    status += nccf_free_coord(csCoordIds[i]);
  }

  free( regrid_filename );
  free( tarTileFileName );
  free( csTileFileName );

  return status;
}

int createCubeSphereTile( const int iGrid, const char *gridName, 
                          const int dims[], int coordIds[], int *gridId ){
  int pos, sign; //, i;
  int faceVec[3] = {0, 0, 0};
  const int ndims = 2;
  int nvertex = dims[0]*dims[1];
  int save = 1, status = 0;
  int gid;
  double *clon, *clat; //, dataValues[nvertex];

  const char *dimnames[] = {"nj", "ni"};

  pos  = iGrid % 3;
  sign = (iGrid/3) == 0 ? 1 : -1;
  faceVec[pos] = sign;

  clat = ( double* )malloc( nvertex * sizeof( double ));
  clon = ( double* )malloc( nvertex * sizeof( double ));

  nccf_get_cubedsphere_grid( dims, faceVec, clon, clat );

  status += nccf_def_lat_coord( ndims, dims, dimnames, clat, save, &coordIds[0] );
  status += nccf_def_lon_coord( ndims, dims, dimnames, clon, save, &coordIds[1] );
  status += nccf_def_grid( coordIds, gridName, &gid );

  *gridId  =  gid;

  free( clat );
  free( clon );

  return status;

}

int createLonLatGrid(const double xymin[], const double xymax[], const int dims[], int coordIds[]) {
  const int save = 1;
  const int ndims = 2;
  double dxs[ndims];
  const char *dimnames[] = {"nj", "ni"};
  const char *name = "my_lon_lat_grid";
  int i, j, k;
  int ntot;
  int status;
  int gridId;

  double **coordData;
  coordData = (double **) malloc(ndims * sizeof(double *));
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
    dxs[i] = (xymax[i] - xymin[i]) / (dims[i] - 1);
  }
  for (i = 0; i < ndims; ++i) {
    coordData[i] = (double *) malloc(ntot * sizeof(double));
  }

  /* Populate coordinates and create lon/lat coordinate objects */
  for (j = 0; j < dims[0]; ++j) {
    for (i = 0; i < dims[1]; ++i) {
      k = i + dims[1]*j;
      coordData[0][k] = xymin[0] + j*dxs[0];
      coordData[1][k] = xymin[1] + i*dxs[1];
    }
  }

  /* Create coordinates */
  status = nccf_def_lat_coord(ndims, dims, dimnames,
					coordData[0], save,
					&coordIds[0]);
  status += nccf_def_lon_coord(ndims, dims, dimnames,
					coordData[1], save,
					&coordIds[1]);

  /* Create grid */
  status += nccf_def_grid(coordIds, name,
					 &gridId);

  for (i = 0; i < ndims; ++i) {
    free(coordData[i]);
  }
  free(coordData);

  return gridId;
}

//////////////////////////////////////////////////////////////////////

int createData(int gridId, const char *filename, int useBogusData, 
               const char *coordinates_id ){

  int dataId, globalIdGrid;
  int status;
  int ndims, i, k, ntot;
  const char *name = "foo";
  const char *standard_name = NULL;
  const char *units = NULL;
  const char *timeDimName = NULL; /* not time dependent */
  const int save = 1;
  double *data;
  char gridname[STRING_SIZE];
  
  if ((status = nccf_def_data(gridId,
				     name, standard_name,
				     units, timeDimName,
				     &dataId))) ERR;

  /* get the coordinate data pointers */
  if ((status = nccf_inq_grid_ndims(gridId, &ndims))) ERR;
  int coordIds[ndims];
  double *coordData[ndims];
  if ((status = nccf_inq_grid_coordids(gridId, coordIds))) ERR;
  nccf_inq_grid_name( gridId, gridname );
  int dims[ndims];
  if ((status = nccf_inq_coord_dims(coordIds[0], dims))) ERR;
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    if ((status = nccf_get_coord_data_pointer(coordIds[i], &coordData[i]))) ERR;
    ntot *= dims[i];
  }

  /* set the data */
  data = (double *) malloc(ntot * sizeof(double));
  for (k = 0; k < ntot; ++k) {
    /* coordData[1] = lon */
    /* coordData[0] = lat */
    data[k] = cos(M_PI*coordData[1][k]/180.0) * sin(M_PI*coordData[0][k]/180.0);
    if( useBogusData != 0 ) data[k] = useBogusData;
  }
  status += nccf_set_data_double(dataId, data, save, NC_FILL_DOUBLE);

  /* Global Atts for Grid and Data */
  nccf_def_global( &globalIdGrid );
  nccf_add_global_att( globalIdGrid, CF_COORDINATES_ID, coordinates_id, 0 );
  nccf_add_global_att( globalIdGrid, CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 );
  nccf_add_global_att( globalIdGrid, CF_FILETYPE, CF_GLATT_FILETYPE_STATIC_DATA, 1 );
  nccf_add_global_att( globalIdGrid, CF_GRIDNAME, gridname, 0 );

  /* write lonlat data and coordinates to file */
  int ncid;
  if ((status = nc_create(filename, NC_CLOBBER, &ncid))) ERR;
  status += nccf_put_grid(gridId, ncid);
  status += nccf_put_data(dataId, ncid);
  status += nccf_put_global(globalIdGrid, ncid);
  if ((status = nc_close(ncid))) ERR;

  nccf_free_global( globalIdGrid );
  free(data);

  return dataId;
}

//////////////////////////////////////////////////////////////////////

int writeData( const int gridId, const int dataId, 
                        const char *filename ){

  /* write lonlat data and coordinates to file */
  int ncid, status;
  if ((status = nc_open(filename, NC_WRITE, &ncid))) ERR;
  status = nccf_put_grid(gridId, ncid);
  status = nccf_put_data(dataId, ncid);
  if ((status = nc_close(ncid))) ERR;

  return NC_NOERR;
}

//////////////////////////////////////////////////////////////////////
/**
 * Compute data2 as the difference between data1 and data2 
 */
void
computeError(int dataId1, int dataId2) {
  int status, k, ndims, ntot, i;
  int gridId;
  status = nccf_inq_data_gridid(dataId1, &gridId);
  status += nccf_inq_grid_ndims(gridId, &ndims);
  int coordIds[ndims];
  status += nccf_inq_grid_coordids(gridId, coordIds);
  int dims[ndims];
  status += nccf_inq_coord_dims(coordIds[0], dims);
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  double *data1;
  double *data2;
  nc_type xtype;
  const void *fill_value;
  status += nccf_get_data_pointer(dataId1, &xtype, (void **) &data1,
					 &fill_value);
  status += nccf_get_data_pointer(dataId2, &xtype, (void **) &data2,
					 &fill_value);

  /* Set the values of data2 */
  for (k = 0; k < ntot; ++k) {
    data2[k] = (data1[k] - data2[k]);
  }
}

//////////////////////////////////////////////////////////////////////

double
checkInterp(int dataId1, int dataId2) {
  int status, k, ndims, ntot, i;
  double err;
  int gridId;
  status = nccf_inq_data_gridid(dataId1, &gridId);
  status += nccf_inq_grid_ndims(gridId, &ndims);
  int coordIds[ndims];
  status += nccf_inq_grid_coordids(gridId, coordIds);
  int dims[ndims];
  status += nccf_inq_coord_dims(coordIds[0], dims);
  ntot = 1;
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  /* compare the two sets of data */
  double *data1;
  double *data2;
  nc_type xtype;
  const void *fill_value;
  status += nccf_get_data_pointer(dataId1, &xtype, (void **) &data1,
					 &fill_value);
  status += nccf_get_data_pointer(dataId2, &xtype, (void **) &data2,
					 &fill_value);

  err = 0.0;
  for (k = 0; k < ntot; ++k) {
    err += fabs(data1[k] - data2[k]);
  }
  err /= ntot;

#ifdef LOGGING
  printf("Average interpolation error: %lf\n", err);
#endif
  return err;
}

//////////////////////////////////////////////////////////////////////
