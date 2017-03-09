/**
 * Test creation of a longitude curvilinear coordinate
 *
 * $Id: tst_lon.c 767 2011-06-06 23:20:19Z pletzer $
 */

#include "nccf_coord.h"
#include <stdlib.h>
#include <assert.h>
#include <netcdf.h>
#include <nccf_handle_error.h>

int main() {

  const int save = 0;
  const int ndims = 2;
  const int dims[] = {10, 11};
  const double lonMin = 0.0;
  const double lonMax = 360.0;
  double dLon;
  double *data;
  int i, j;
  int coordid, coordid2;
  int status;
  int ncid;

  int ndims2;
  int *dims2;
  double *data2;
  
  const char *dimnames[] = {"ni", "nj"};

  data = (double *) malloc(sizeof(double) * dims[0] * dims[1]);
  dLon = (lonMax - lonMin) / (dims[0] - 1);
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      data[i + dims[0] * j] = lonMin + i * dLon;
    }
  }
  
  /* Define */

  if ((status = nccf_def_lon_coord( ndims, dims, dimnames, 
                                         data, save, &coordid))) ERR;

  /* Add an attribute */
  if ((status = nccf_add_coord_att(coordid, "to", 
					"be or not to be"))) ERR;

  /* Write to file */

  if ((status = nc_create("tst_lon.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_coord(coordid, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  /* Read from file */
  
  if ((status = nccf_def_coord_from_file("tst_lon.nc", 
					      "lon", &coordid2))) ERR;

  /* Check */

  if ((status = nccf_inq_coord_ndims(coordid, &ndims2))) ERR;
  assert(ndims2 == ndims);
  if ((status = nccf_inq_coord_ndims(coordid2, &ndims2))) ERR;
  assert(ndims2 == ndims);  

  dims2 = (int *) malloc(ndims2 * sizeof(int));
  if ((status = nccf_inq_coord_dims(coordid, dims2))) ERR;
  for (i = 0; i < ndims; ++i) {
    assert(dims2[i] == dims[i]);
  }
  if ((status = nccf_inq_coord_dims(coordid2, dims2))) ERR;
  for (i = 0; i < ndims; ++i) {
    assert(dims2[i] == dims[i]);
  }
  free(dims2);

  if ((status = nccf_get_coord_data_pointer(coordid, &data2))) ERR;
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      assert(data[i + dims[0] * j] == data2[i + dims[0] * j]);
    }
  }
  if ((status = nccf_get_coord_data_pointer(coordid2, &data2))) ERR;
  for (j = 0; j < dims[1]; ++j) {
    for (i = 0; i < dims[0]; ++i) {
      assert(data[i + dims[0] * j] == data2[i + dims[0] * j]);
    }
  }
  
  /* Free */

  if ((status = nccf_free_coord(coordid))) ERR;
  if ((status = nccf_free_coord(coordid2))) ERR;
  free(data);
  return status;
}
