/**
 * $Id: exparinterp.c 901 2011-12-23 20:33:15Z pletzer $
 *
 * Example of code to perform interpolation in parallel
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <examples.h>

#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

double cost_function(int ny, int nx) {
  // ~ ratio of surface over volume
  return (double)(ny + nx)/(double)(ny * nx);
}

void create_domain_decomp(int nprocs, int *ny, int *nx) {
  // create an optimal 2D domain decomposition, one 
  // one which minimizes the ratio of surface over volume
  int i = 1;
  int j = nprocs;
  double newCost;
  double cost = cost_function(i, j);
  *ny = i; 
  *nx = j;
  for (i = 2; i < nprocs/2; ++i) {
    j = nprocs / i;
    if ( (double)(nprocs)/(double)(i) == j) {
      newCost = cost_function(i, j);
      if (newCost < cost) {
        *ny = i;
        *nx = j;
        cost = newCost;
      }
    }
  }
}

double data_fct_tgt(double lat, double lon) {
  return NC_FILL_DOUBLE;
}

double data_fct(double lat, double lon) {
  return sin(3*M_PI*lat/180.0)*cos(5*M_PI*lon/180.0);
}

int check_data(int dataid, double (*data_fct)(double, double), 
               double tolfct) {
  int fails = 0;
  int dims[2];
  int status;
  int i, j;
  status = nccf_inq_data_dims(dataid, dims);
  int gridid;
  status = nccf_inq_data_gridid(dataid, &gridid);
  int coordids[2];
  status = nccf_inq_grid_coordids(gridid, coordids);
  double *lats;
  double *lons;
  status = nccf_get_coord_data_pointer(coordids[0], &lats);
  status = nccf_get_coord_data_pointer(coordids[1], &lons);
  nc_type type;
  double *data;
  double *missing_value;
  status = nccf_get_data_pointer(dataid, &type, (void **) &data, 
                                 (const void **) &missing_value);
  int index;
  for (j = 0; j < dims[0]; ++j) {
    for (i = 0; i < dims[1]; ++i) {
      index = i + dims[1]*j;
      double fexact = data_fct(lats[index], lons[index]);
      double fintrp = data[index];
      if (fabs(fexact - fintrp) > tolfct) {
        fails++;
        printf("interpolation error %f for lat, lon = %f, %f\n", 
               fabs(fexact - fintrp), 
               lats[index], lons[index]);
      }
    }
  }
  return fails;
}

int create_latlon(const int dims[], 
                  const double coord_mins[], const double coord_maxs[], 
                  int coordids[], int *gridid) {
  const char *dimnames[NC_MAX_NAME+1] = {"nlat", "nlon"};
  int j, i, status;
  double *lats = (double *) malloc(dims[0]*dims[1]*sizeof(double));
  double *lons = (double *) malloc(dims[0]*dims[1]*sizeof(double));
  double dlat = (coord_maxs[0] - coord_mins[0]) / ((double) dims[0]-1);
  double dlon = (coord_maxs[1] - coord_mins[1]) / ((double) dims[1]-1);
  for (j = 0; j < dims[0]; ++j) {
    for (i = 0; i < dims[1]; ++i) {
      int index = i + dims[1]*j;
      lats[index] = coord_mins[0] + j*dlat;
      lons[index] = coord_mins[1] + i*dlon;
    }
  }
  const int save = 1;
  status = nccf_def_lat_coord(2, dims, dimnames, lats, save, &coordids[0]);
  status = nccf_def_lon_coord(2, dims, dimnames, lons, save, &coordids[1]);
  status = nccf_def_grid(coordids, "latlon", gridid);
  free(lats);
  free(lons);
  return status;
}

/**
 * Create data object and fill in values
 * \param gridid grid ID
 * \param fct function of (lat, lon), can be NULL in which case data will 
 *            be set to missing value
 * \param dataid return ID
 */
int create_data(int gridid, double (*fct)(double, double), int *dataid) {
  int status;
  const char *standard_name = "";
  const char *units = "";
  // static data
  status = nccf_def_data(gridid, "data", standard_name, units, NULL, 
                         dataid);
  if (fct) {
    int dims[2];
    status = nccf_inq_data_dims(*dataid, dims);
    int coordids[2];
    status = nccf_inq_grid_coordids(gridid, coordids);
    double *lats, *lons;
    status = nccf_get_coord_data_pointer(coordids[0], &lats);
    status = nccf_get_coord_data_pointer(coordids[1], &lons);
    double *data = (double *) malloc(dims[0]*dims[1]*sizeof(double));
    int i, j, index;
    for (j = 0; j < dims[0]; ++j) {
      for (i = 0; i < dims[1]; ++i) {
        index = i + dims[1]*j;
        data[index] = fct(lats[index], lons[index]);
      }
    }
    const int save = 1;
    status = nccf_set_data_double(*dataid, data, save, NC_FILL_DOUBLE);
    free(data);
  }
  return status;
}

//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

  int status;
  int nprocs = 1; 
  int procid = 0;
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &procid);
#endif
  if (procid == 0) {
    printf("\nRegridding from a rectilinear grid to another in parallel.\n\n");
    printf("MPI number of procs: %d\n", nprocs);
  }

  // source grid
  const int nlat = 180 + 1;
  const int nlon = 360 + 1;
  const int dims[] = {nlat, nlon};
  const double coord_mins[] = {-90.0, -180.0};
  const double coord_maxs[] = {+90.0, +180.0};
  int coordids[2];
  int gridid;
  create_latlon(dims, coord_mins, coord_maxs, 
                coordids, &gridid);
  int dataid;
  create_data(gridid, data_fct, &dataid);

  // target grids, different for each proc
  const int nlat_tgt = 2*50 + 1;
  const int nlon_tgt = 2*60 + 1;
  const int dims_tgt[] = {nlat_tgt, nlon_tgt};
  int nx, ny;
  create_domain_decomp(nprocs, &ny, &nx);
  if (procid == 0) {
    printf("source grid is                              %d x %d\n", nlat, nlon);
    printf("domain decomposition of destination grid is %d x %d\n", ny, nx);
    printf("destination grid is                         %d x %d\n", nlat_tgt, nlon_tgt);
    printf("total number of interpolation points is     %d\n", nprocs*nlat_tgt*nlon_tgt);
  }
  double dy_dom = (coord_maxs[0] - coord_mins[0])/( (double)(ny) );
  double dx_dom = (coord_maxs[1] - coord_mins[1])/( (double)(nx) );
  int j = procid / nx;
  int i = procid % nx;
  double coord_mins_tgt[] = {coord_mins[0] + j*dy_dom, 
                             coord_mins[1] + i*dx_dom};
  double coord_maxs_tgt[] = {coord_mins[0] + (j+1)*dy_dom, 
                             coord_mins[1] + (i+1)*dx_dom};
  printf("[%d] min/max target corners are: %f %f / %f %f\n", procid,
         coord_mins_tgt[0], coord_mins_tgt[1], 
         coord_maxs_tgt[0], coord_maxs_tgt[1]);
  int coordids_tgt[2];
  int gridid_tgt;
  create_latlon(dims_tgt, coord_mins_tgt, coord_maxs_tgt, 
                coordids_tgt, &gridid_tgt);
  int dataid_tgt;
  create_data(gridid_tgt, data_fct_tgt, &dataid_tgt);
  
  // interpolation
  int regridid;
  if ((status = nccf_def_regrid(gridid, gridid_tgt, &regridid))) ERR;
  const int nitermax = 5;
  const double tolpos = 1.e-6;
  if ((status = nccf_compute_regrid_weights(regridid, nitermax, tolpos))) ERR;
  if ((status = nccf_apply_regrid(regridid, dataid, dataid_tgt))) ERR;
  
  // check
  const double tolfct = 0.01;
  int nfails = check_data(dataid_tgt, data_fct, tolfct);
  if (nfails > 0) {
    printf("[%d] number of interpolation failures: %d\n", procid, nfails);
  }
  int ntotfails = nfails;
#ifdef HAVE_MPI_H
  MPI_Reduce(&nfails, &ntotfails, 1, MPI_INT, MPI_SUM, 0, comm);
#endif
  if (procid == 0) {
    printf("total number of interpolation failures: %d\n", ntotfails);
  }

  // clean up
  if ((status = nccf_free_regrid(regridid))) ERR;
  if ((status = nccf_free_data(dataid))) ERR;
  if ((status = nccf_free_data(dataid_tgt))) ERR;
  if ((status = nccf_free_grid(gridid))) ERR;
  if ((status = nccf_free_grid(gridid_tgt))) ERR;
  if ((status = nccf_free_coord(coordids_tgt[0]))) ERR;
  if ((status = nccf_free_coord(coordids_tgt[1]))) ERR;
  if ((status = nccf_free_coord(coordids[0]))) ERR;
  if ((status = nccf_free_coord(coordids[1]))) ERR;

#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
  return 0;
}
