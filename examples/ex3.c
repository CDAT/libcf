/**
 * Example showing how to generate cell centered data of two cubed sphere
 * tiles. 
 *
 * $Id: ex3.c 831 2011-09-14 20:24:21Z dkindig $
 */

// std includes
#include <stdio.h>
#include <math.h>

#include <netcdf.h>

#include <nccf_utility_functions.h>
#include <nccf_global.h>
#include <nccf_coord.h>
#include <nccf_grid.h>
#include <nccf_data.h>
#include <nccf_regrid.h>
#include <nccf_mosaic.h>
#include <nccf_host.h>

#include <nccf_handle_error.h>

// function for setting data
double funct(double lon, double lat) {
  // some test field
  return sin(5*M_PI*(lon - 30.0)/180.0 - 0.3) * cos(4*M_PI*lat/180.0);
}

double correct_for_dateline(double dlon) {
  const double maxjump = 200.0;
  const double period = 360.0;
  double res = dlon;
  if (dlon < -maxjump) {
    res += period;
  }
  if (dlon > +maxjump) {
    res -= period;
  }
  return res;
}

int main(int argc, char *argv[]) {

  int status;
  int i, j;
  int k, k0;

  const char *coordinates_id = "db9b6740-25a2-11e0-a724-5c260a1834c1";
  const char *data_id = "e906f520-25a2-11e0-beeb-5c260a1834c1";
  const char *case_name = "gsex3";

  // space dimensionality
  const int ndims = 2;

  // set the tile resolutions
  const int n1 = 20;
  const int dims[] = {n1, n1};
  const int dims_cell[] = {n1-1, n1-1};

  // number of tiles
  const int ntiles = 3;

  // set the tile orientation on the sphere (normal vector pointing out)
  int norm_vect[ntiles][3];
  if (ntiles >= 1) {
    norm_vect[0][0] = 1; norm_vect[0][1] = 0; norm_vect[0][2] = 0;
  }
  if (ntiles >= 2) {
  norm_vect[1][0] = 0; norm_vect[1][1] = 1; norm_vect[1][2] = 0;
  }
  if (ntiles >= 3) {
    norm_vect[2][0] = 0; norm_vect[2][1] = 0; norm_vect[2][2] = 1;
  }

  // generate nodal coordinate data
  double **lon;
  double **lat;
  lon =  (double **) malloc(ntiles * sizeof(double *));
  lat =  (double **) malloc(ntiles * sizeof(double *));
  for (i = 0; i < ntiles; ++i) {
    lon[i] = (double *) malloc( n1*n1 * sizeof(double));
    lat[i] = (double *) malloc( n1*n1 * sizeof(double));
  }
  for (i = 0; i < ntiles; ++i) {
    if ((status = nccf_get_cubedsphere_grid(dims, norm_vect[i],
  					    lon[i], lat[i]))) ERR;
  }

  // generate cell centered coordinate data and scalar field
  double **lonCell;
  double **latCell;
  double **field;
  lonCell =  (double **) malloc(ntiles * sizeof(double *));
  latCell =  (double **) malloc(ntiles * sizeof(double *));
  field =  (double **) malloc(ntiles * sizeof(double *));
  for (i = 0; i < ntiles; ++i) {
    lonCell[i] = (double *) malloc( (n1-1)*(n1-1) * sizeof(double));
    latCell[i] = (double *) malloc( (n1-1)*(n1-1) * sizeof(double));
    field[i] = (double *) malloc( (n1-1)*(n1-1) * sizeof(double));
    int dims_cell[] = {n1 - 1, n1 - 1};
    int dims[] = {n1, n1};
    int ij[ndims];
    int ntot_cell = dims_cell[0] * dims_cell[1];
    double lon0, dlon;
    double lat0;
    for (k = 0; k < ntot_cell; ++k) {
      // compute ij
      nccf_get_multi_index(ndims, dims_cell, k, ij);
      // nodal flat index
      k0 = nccf_get_flat_index(ndims, dims, ij);
      lon0 = lon[i][k0];
      lat0 = lat[i][k0];
      // average nodal values, correct for dateline if necessary
      lonCell[i][k] = lon0;
      latCell[i][k] = lat0;

      ij[0] += 1;
      k0 = nccf_get_flat_index(ndims, dims, ij);
      dlon = correct_for_dateline(lon[i][k0] - lon0);
      lonCell[i][k] += (lon0 + dlon);
      latCell[i][k] += lat[i][k0];

      ij[1] += 1;
      k0 = nccf_get_flat_index(ndims, dims, ij);
      dlon = correct_for_dateline(lon[i][k0] - lon0);
      lonCell[i][k] += (lon0 + dlon);
      latCell[i][k] += lat[i][k0];
      
      ij[0] -= 1;
      k0 = nccf_get_flat_index(ndims, dims, ij);
      dlon = correct_for_dateline(lon[i][k0] - lon0);
      lonCell[i][k] += (lon0 + dlon);
      latCell[i][k] += lat[i][k0];

      // now average
      lonCell[i][k] /= 4.0;
      latCell[i][k] /= 4.0;
      // set the field
      field[i][k] = funct(lonCell[i][k], latCell[i][k]);
    }
  }

  // generate grid and data objects
  int gridids[ntiles];
  int gridids_cell[ntiles];
  int dataids_cell[ntiles];
  int coordids[ntiles][ndims];
  int coordids_cell[ntiles][ndims];
  int cellglobalid[ntiles];
  for (i = 0; i < ntiles; ++i) {
    int const save = 0;
    const char *dimnames[] = {"n_x", "n_y"};
    const char *dimnames_cell[] = {"ncell_x", "ncell_y"};
    if ((status = nccf_def_lon_coord(ndims, dims,
  				     dimnames, lon[i], save,
  				     &coordids[i][0]))) ERR;
    if ((status = nccf_def_lat_coord(ndims, dims,
  				     dimnames, lat[i], save,
  				     &coordids[i][1]))) ERR;
    char tilename[STRING_SIZE];
    sprintf(tilename, "tile%d", i);
    if ((status = nccf_def_grid(coordids[i], tilename,
  				       &gridids[i]))) ERR;
    
    // cell center coordinates
    if ((status = nccf_def_coord(ndims, dims_cell, dimnames_cell, 
				 lonCell[i], save, "lon_cell", 
				 CF_COORD_LON_STNAME, CF_COORD_LON_UNITS,
				 &coordids_cell[i][0]))) ERR;
    if ((status = nccf_def_coord(ndims, dims_cell, dimnames_cell, 
				 latCell[i], save, "lat_cell", 
				 CF_COORD_LAT_STNAME, CF_COORD_LAT_UNITS,
				 &coordids_cell[i][1]))) ERR;
    if ((status = nccf_def_grid(coordids_cell[i], tilename,
  				       &gridids_cell[i]))) ERR;
    const char *time_axis = ""; // static data
    if ((status = nccf_def_data(gridids_cell[i], "some_data",
  				       "put_standard_name_here", 
				       "put_units_here", 
				         time_axis,
  				       &dataids_cell[i]))) ERR;
    if ((status = nccf_def_global( &cellglobalid[i] )));
    if ((status = nccf_add_global_att(cellglobalid[i], CF_DATA_ID, data_id, 0 ))) ERR;
    if ((status = nccf_add_global_att(cellglobalid[i], CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;
    if ((status = nccf_add_global_att(cellglobalid[i], CF_FILETYPE, CF_GLATT_FILETYPE_STATIC_DATA, 0 ))) ERR;
    if ((status = nccf_add_global_att(cellglobalid[i], CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 1 ))) ERR;
    if ((status = nccf_add_global_att(cellglobalid[i], CF_GRIDNAME, tilename, 0 ))) ERR;
    if ((status = nccf_set_data_double(dataids_cell[i], 
					      field[i], save, 
					      NC_FILL_DOUBLE))) ERR;
  }

  // generate the mosaic
  int mosaicid, mosaicglobalid;
  const double periods[] = {360.0, 0.0};
  if ((status = nccf_def_mosaic(ntiles, gridids, "mosaic", &mosaicid))) ERR;
  if ((status = nccf_compute_mosaic_contacts( mosaicid, periods ))) ERR;
  if ((status = nccf_def_global( &mosaicglobalid ))) ERR;
  if ((status = nccf_add_global_att(mosaicglobalid,
					   CF_COORDINATES_ID, 
					   coordinates_id, 0 ))) ERR;
  if ((status = nccf_add_global_att(mosaicglobalid,
					   CF_FILETYPE, 
					   CF_GLATT_FILETYPE_MOSAIC, 0 ))) ERR;

  // host
  const int ntime_slices = 1;
  int hostid, hostglobalid;
  if ((status = nccf_def_host(coordinates_id, data_id, ntime_slices, 
			      &hostid))) ERR;
  const int force = 0;
  if ((status = nccf_def_global( &hostglobalid ))) ERR;
  if ((status = nccf_add_global_att(hostglobalid, CF_FILETYPE, 
                                    CF_GLATT_FILETYPE_HOST, 0 ))) ERR;
  if ((status = nccf_add_global_att(hostglobalid, CF_COORDINATES_ID, 
                                    coordinates_id, 0 ))) ERR;
  if ((status = nccf_add_global_att(hostglobalid, CF_DATA_ID, 
                                    data_id, 0 ))) ERR;

  // write everything to file
  int ncid;
  char filename[STRING_SIZE];
  for (i = 0; i < ntiles; ++i) {
    // data and grid in the same file
    sprintf(filename, "%s_tile%d.nc", case_name, i);
    if ((status = nc_create(filename, NC_CLOBBER, &ncid))) ERR;
    if ((status = nccf_put_grid(gridids[i], ncid))) ERR;
    if ((status = nccf_put_grid(gridids_cell[i], ncid))) ERR;
    if ((status = nccf_put_data(dataids_cell[i], ncid))) ERR;
    if ((status = nccf_put_global(cellglobalid[i], ncid)));
    if ((status = nc_close(ncid))) ERR;
    if ((status = nccf_add_host_file(hostid, filename, force))) ERR;
  }
  // mosaic
  sprintf(filename, "%s_mosaic.nc", case_name);
  if ((status = nc_create(filename, NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if ((status = nccf_put_global(mosaicglobalid, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;
  if ((status = nccf_add_host_file(hostid, filename, force))) ERR;

  // host
  sprintf(filename, "%s_host.nc", case_name);
  if ((status = nc_create(filename, NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_host(hostid, ncid))) ERR;
  if ((status = nccf_put_global(hostglobalid, ncid))) ERR;
  if ((status = nc_close(ncid))) ERR;

  // clean up
  for (i = 0; i < ntiles; ++i) {
    if ((status = nccf_free_global(cellglobalid[i]))) ERR;
    if ((status = nccf_free_data(dataids_cell[i]))) ERR;
    if ((status = nccf_free_grid(gridids_cell[i]))) ERR;
    if ((status = nccf_free_grid(gridids[i]))) ERR;
    for (j = 0; j < ndims; ++j) {
      if ((status = nccf_free_coord(coordids_cell[i][j]))) ERR;
      if ((status = nccf_free_coord(coordids[i][j]))) ERR;
    }
    free(lonCell[i]);
    free(latCell[i]);
    free(lon[i]);
    free(lat[i]);
    free(field[i]);
  }
  if ((status = nccf_free_global(mosaicglobalid))) ERR;
  if ((status = nccf_free_global(hostglobalid))) ERR;
  if ((status = nccf_free_mosaic(mosaicid))) ERR;
  if ((status = nccf_free_host(hostid))) ERR;
  free(lonCell);
  free(latCell);
  free(lon);
  free(lat);
  free(field);
  
  return 0;
}
