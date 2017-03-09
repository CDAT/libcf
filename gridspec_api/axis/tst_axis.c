/*
 * Test axis creation
 * $Id: tst_axis.c 768 2011-06-07 03:32:37Z pletzer $
 */

#include <nccf_axis.h>
#include <stdlib.h>
#include <assert.h>
#include <netcdf.h>
#include <nccf_handle_error.h>

int main() {
  int status, axisid, i, ncid, axisid2, ncid2;

  const char standard_name[] = "longitude";
  const char units[] = "degrees_east";
  const int cdm_axis_type = NCCF_LONGITUDE;
  const char axis[] = "X";
  const int positive_up = 1;
  const char *formula_terms = NULL;
  
  const int n = 11;
  float *data = (float *) malloc(n * sizeof(float));
  for (i = 0; i < n; ++i) {
    data[i] = 0.0f + (360.0f/(float)(n - 1)) * (float)(i);
  }

  if ((status = nccf_def_axis("lon", n, NC_FLOAT, data, standard_name, 
                              units, cdm_axis_type, axis, positive_up, 
                              formula_terms, &axisid))) ERR;

  if ((status = nccf_add_axis_att(axisid, "axis", "X"))) ERR;
  

  if ((status = nc_create("tst_axis.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_axis(axisid, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  if ((status = nccf_def_axis_from_file("tst_axis.nc", "lon", &axisid2))) ERR;

  if ((status = nc_create("tst_axis2.nc", NC_CLOBBER, &ncid2))) ERR;
  if ((status = nccf_put_axis(axisid2, ncid2))) ERR;
  if ((status = nc_close(ncid2))) ERR;

  /* checks */
  int len, len2;
  if ((status = nccf_inq_axis_len(axisid, &len))) ERR;
  if ((status = nccf_inq_axis_len(axisid, &len2))) ERR;
  assert(len == len2);

  if ((status = nccf_free_axis(axisid))) ERR;
  if ((status = nccf_free_axis(axisid2))) ERR;
  free(data);
  return 0;
}
