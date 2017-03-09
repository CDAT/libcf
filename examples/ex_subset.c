/*
 * Example of data subsetting
 * $Id: ex_subset.c 881 2011-12-17 21:53:14Z pletzer $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <netcdf.h>

#include <nccf_utility_functions.h>
#include <nccf_global.h>
#include <nccf_coord.h>
#include <nccf_grid.h>
#include <nccf_data.h>
#include <nccf_regrid.h>

/**
 * Subsetting a field
 * \param srcdata_id data id on source grid
 * \param dstgrid_id destination
 * \param subsampled_data will set the pointer to the sub-sampled data
 *        on return. Caller responsible for free'ing the data
 * \return 0 on success, or a positive number equial to the numeber points
 *        of points that could not be interpolated (likely because these
 *        fell out of the domain, or perhaps due to singularities in the 
 *        coordinate system). A negative number indicates other failures.
 */
int subset(int srcdata_id, int dstgrid_id, float* subsampled_data) {

  int i;
  int status;

  // get the grid sizes of dstgrid_id
  int ndims;
  status = nccf_inq_grid_ndims(dstgrid_id, &ndims);
  if (status) ERR;
  if (ndims <= 0) return -2;

  // get the grid id attached to the source data
  int srcgrid_id;
  status = nccf_inq_data_gridid(srcdata_id, &srcgrid_id);
  if (status) ERR;
  
  // define the regrid object 
  int regrid_id;
  status = nccf_def_regrid(srcgrid_id, dstgrid_id, &regrid_id);
  if (status) ERR;

  // compute the interpolation weights
  // caller may want to pass parameters down this routine....
  const int nitermax = 100;
  const double tolpos = 0.01;
  status = nccf_compute_regrid_weights(regrid_id, nitermax, tolpos);
  if (status) ERR;

  int dstcoordids[ndims];
  status = nccf_inq_grid_coordids(dstgrid_id, dstcoordids);
  if (status) ERR;
  int dstdims[ndims];
  status = nccf_inq_coord_dims(dstcoordids[0], dstdims);
  if (status) ERR;
  int dstntot = 1;
  for (i = 0; i < ndims; ++i) {
    dstntot *= dstdims[i];
  }

  // create container for the interpolated data, no standard_name, 
  // no units, and no time dependence
  int dstdata_id;
  status = nccf_def_data(dstgrid_id, "temporary", "", "", NULL, &dstdata_id);
  if (status) ERR;

  // initialize the interpolated data
  float *dstdata = malloc(dstntot * sizeof(float));
  for (i = 0; i < dstntot; ++i) {
      dstdata[i] = 0.0f;
  }
  const int not_save = 0;
  status = nccf_set_data_float(dstdata_id, dstdata, not_save, NC_FILL_FLOAT);
  if (status) ERR;
 
  // interpolate
  status = nccf_apply_regrid(regrid_id, srcdata_id, dstdata_id);
  if (status) ERR;

  // copy
  for (i = 0; i < dstntot; ++i) {
    subsampled_data[i] = dstdata[i];
  }

  // cleanup
  status = nccf_free_regrid(regrid_id);
  if (status) ERR;
  status = nccf_free_data(dstdata_id);
  if (status) ERR;
  free(dstdata);

  return 0;
}

/**
 * Function to set data 
 * @param lat latitude
 * @param lon longitude
 * @return value
 */
float funct(double lat, double lon) {
  return 292.0 + cosf(2*lon*M_PI/180.0)*sinf(3*lat*M_PI/180.0);
}

/* =================================================================== */

int main() {

  const int ndims = 2;
  int i, j, k, status;
  int srclon_id, srclat_id;
  
  // create source coordinates and grid
  const int nlon = 181;
  const int nlat = 91;
  const double lonmin = 0.0; 
  const double lonmax = 360.0; 
  const double latmin = -90.0; 
  const double latmax = +90.0; 
  const double dlon = (lonmax - lonmin)/(nlon-1);
  const double dlat = (latmax - latmin)/(nlat-1);
  double *lon = malloc(nlat * nlon * sizeof(double));
  double *lat = malloc(nlat * nlon * sizeof(double));
  for (j = 0; j < nlat; ++j) {
    for (i = 0; i < nlon; ++i) {
      k = i + j*nlon;
      lon[k] = lonmin + i*dlon;
      lat[k] = latmin + j*dlat;
    }
  }
  const int srcdims[] = {nlat, nlon};
  const char *srcdimnames[] = {"nlat", "nlon"};
  const int save = 1;
  status = nccf_def_lon_coord(ndims, srcdims, srcdimnames,
                              lon, save, &srclon_id);
  if (status) ERR;
  status = nccf_def_lat_coord(ndims, srcdims, srcdimnames,
                              lat, save, &srclat_id);
  if (status) ERR;
  const int srccoord_ids[] = {srclat_id, srclon_id};
  int srcgrid_id;
  status = nccf_def_grid(srccoord_ids, "lat_lon", &srcgrid_id);
  if (status) ERR;

  // create static data
  const char* time_dimname = NULL;
  int srcta_id;
  status = nccf_def_data(srcgrid_id, "ta", "air_temperature", 
                         "Kelvin", time_dimname, &srcta_id);
  if (status) ERR;

  // set the data to some values, floats
  float *ta = malloc(nlat * nlon * sizeof(float));
  for (j = 0; j < nlat; ++j) {
    for (i = 0; i < nlon; ++i) {
      k = i + j*nlon;
      ta[k] = funct(lat[k], lon[k]);
    }
  }
  const int not_save = 0;
  status = nccf_set_data_float(srcta_id, ta, not_save, NC_FILL_FLOAT);
  if (status) ERR;
  
  // create destination coordinates and grid
  const int nlon_dest = 11;
  const int nlat_dest = 6;
  const double dstlonmin = 20.0;
  const double dstlonmax = 40.0;
  const double dstlatmin = 10.0;
  const double dstlatmax = 50.0;
  const double dstdlon = (dstlonmax - dstlonmin)/(nlon_dest - 1);
  const double dstdlat = (dstlatmax - dstlatmin)/(nlat_dest - 1);
  double *dstlon = malloc(nlon_dest * nlat_dest * sizeof(double));
  double *dstlat = malloc(nlon_dest * nlat_dest * sizeof(double));
  for (j = 0; j < nlat_dest; ++j) {
    for (i = 0; i < nlon_dest; ++i) {
      k = i + j*nlon_dest;
      dstlon[k] = dstlonmin + i*dstdlon;
      dstlat[k] = dstlatmin + j*dstdlat;
    }
  }
  const int dstdims[] = {nlat_dest, nlon_dest};
  const char *dstdimnames[] = {"nlat_dest", "nlon_dest"};
  int dstlon_id, dstlat_id;
  status = nccf_def_lon_coord(ndims, dstdims, dstdimnames,
                              dstlon, save, &dstlon_id);
  if (status) ERR;
  status = nccf_def_lat_coord(ndims, dstdims, dstdimnames,
                              dstlat, save, &dstlat_id);
  if (status) ERR;
  const int dstcoord_ids[] = {dstlat_id, dstlon_id};
  int dstgrid_id;
  status = nccf_def_grid(dstcoord_ids, "lat_lon_dest", &dstgrid_id);
  if (status) ERR;

  //
  // sub-sample
  //

  float *dstta = malloc(nlat_dest * nlon_dest * sizeof(float));
  status = subset(srcta_id, dstgrid_id, dstta);
  if (status != 0) {
    printf("Subsetting error: %d\n", status);
    return 1;
  }

  // check
  float error = 0.0f;
  for (j = 0; j < nlat_dest; ++j) {
    for (i = 0; i < nlon_dest; ++i) {
      k = i + j*nlon_dest;
      error += fabs(dstta[k] - funct(dstlat[k], dstlon[k]));
    }
  }
  error /= (nlon_dest * nlat_dest);
  printf("average interpolation error = %g\n", error);
  

  // clean up
  status = nccf_free_data(srcta_id); if(status) ERR;
  free(ta);
  status = nccf_free_grid(srcgrid_id); if(status) ERR;
  status = nccf_free_coord(srclat_id); if(status) ERR;
  status = nccf_free_coord(srclon_id); if(status) ERR;
  free(lat);
  free(lon);
  status = nccf_free_grid(dstgrid_id); if(status) ERR;
  free(dstlat);
  free(dstlon);
  status = nccf_free_coord(dstlat_id); if(status) ERR;
  status = nccf_free_coord(dstlon_id); if(status) ERR;
  free(dstta);

  // done
  return 0;
}
