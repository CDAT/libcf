/**
 * Unit test exercising time dependent data creation
 *
 * $Id: tst_timedata.c 767 2011-06-06 23:20:19Z pletzer $
 */

#include <nccf_data.h>
#include <assert.h>
#include <math.h>

#include <nccf_utility_functions.h>
#include <nccf_coord.h>
#include <nccf_grid.h>
#include <nccf_handle_error.h>

int main() {
  
  int status;
  const int ndims = 3;
  int coordIds[ndims];
  const int save = 0;
  int i, j, k, bigI;
  char coordinates_id[36+1];

  /* generate unique coordinate id */
  if ((status = nccf_generate_id(123, coordinates_id))) ERR;

  /* create coordinates */
  const int nx1 = 11;
  const int ny1 = 21;
  const int nz1 = 31;
  const int ntot = nx1 * ny1 * nz1;
  const double xmin = 0.0, xmax = 360.0;
  const double ymin = -90.0, ymax = 90.0;
  const double zmin = 0.0, zmax = 10000.0;
  double dx = (xmax - xmin) / (double)(nx1 - 1);
  double dy = (ymax - ymin) / (double)(ny1 - 1);
  double dz = (zmax - zmin) / (double)(nz1 - 1);
  
  double *xs = (double *) malloc(ntot * sizeof(double));
  double *ys = (double *) malloc(ntot * sizeof(double));
  double *zs = (double *) malloc(ntot * sizeof(double));
  float *ps = (float *) malloc(ntot * sizeof(float));

  for (k = 0; k < nz1; ++k) {
    double z = zmin + dz * k;
    for (j = 0; j < ny1; ++j) {
      double y = ymin + dy * j;
      for (i = 0; i < nx1; ++i) {
	double x = xmin + dx * i;
	bigI = i + nx1*( j + ny1*k );
	xs[bigI] = x;
	ys[bigI] = y;
	zs[bigI] = z;
	ps[bigI] = 0.0; 
      }
    }
  }
  const int dims[] = {nx1, ny1, nz1};
  const char *dimNames[] = {"nx1", "ny1", "nz1"};
  status = nccf_def_coord(ndims, dims, dimNames, xs, save,
			       "x", "grid_longitude", "degree", &coordIds[0]);
  if (status) ERR;
  status = nccf_def_coord(ndims, dims, dimNames, ys, save,
			       "y", "grid_latitude", "degree", &coordIds[1]);
  if (status) ERR;
  status = nccf_def_coord(ndims, dims, dimNames, zs, save,
			       "z", "sea_surface_elevation", "m", &coordIds[2]);
  if (status) ERR;

  /* create grid */
  int gridid;
  status = nccf_def_grid(coordIds, "myGrid", &gridid);

  /* create data */
  int pId;
  status = nccf_def_data(gridid, "p", "air_pressure", "Pa", "time", &pId);
  if (status) ERR;
  status = nccf_set_data_float(pId, ps, save, NC_FILL_FLOAT);
  if (status) ERR;

  /* write everything to disk */
  int ncid;
  if ((status = nc_create("tst_timedata.nc", NC_CLOBBER, &ncid))) ERR;
  if ((status = nccf_put_grid(gridid, ncid))) ERR;

  const int nt = 3;
  int it;
  for (it = 0; it < nt; ++it) {
    for (k = 0; k < nz1; ++k) {
      double z = zmin + dz * k;
      for (j = 0; j < ny1; ++j) {
	double y = ymin + dy * j;
	for (i = 0; i < nx1; ++i) {
	  double x = xmin + dx * i;
	  bigI = i + nx1*( j + ny1*k );
	  xs[bigI] = x;
	  ys[bigI] = y;
	  zs[bigI] = z;
	  ps[bigI] = it; 
	}
      }
    }
    if ((status = nccf_put_data(pId, ncid))) ERR;
  }
  if ((status = nc_close(ncid))) ERR;
  
  /* clean up */
  if ((status = nccf_free_data(pId))) ERR;
  if ((status = nccf_free_grid(gridid))) ERR;
  for (i = 0; i < ndims; ++i) {
    if ((status = nccf_free_coord(coordIds[i]))) ERR;
  }
  free(xs);
  free(ys);
  free(zs);
  free(ps);

  return 0;
}
