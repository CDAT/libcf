/**
 * Create a host file aggregation comprising static and time dependent 
 * data stored in different files on three tiles of a cubed-sphere. The time 
 * files contain multiple time steps. Time dependent and static data exist on 
 * multiple tiles.
 * 
 * "$Id: tst_cube_sphere_multi_dims_in_file.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_host.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <netcdf.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"
#include "nccf_mosaic.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

#define NAME_MOSAIC "tst_cube_sphere_multi_var_dims_in_file_mosaic"
#define NAME_GRID_0 "tst_cube_sphere_multi_var_dims_in_file_grid0"
#define NAME_GRID_1 "tst_cube_sphere_multi_var_dims_in_file_grid1"
#define NAME_GRID_2 "tst_cube_sphere_multi_var_dims_in_file_grid2"
#define NAME_GRID_0_FILE "tst_cube_sphere_multi_var_dims_in_file_grid0.nc"
#define NAME_GRID_1_FILE "tst_cube_sphere_multi_var_dims_in_file_grid1.nc"
#define NAME_GRID_2_FILE "tst_cube_sphere_multi_var_dims_in_file_grid2.nc"
#define NAME_MOSAIC_FILE "tst_cube_sphere_multi_var_dims_in_file_mosaic.nc"
#define NAME_HOST_FILE   "tst_cube_sphere_multi_var_dims_in_file_host.nc"
#define NAME_STATIC_DATA "tst_cube_sphere_multi_var_dims_in_file_stat_data"
#define NAME_TIME_DATA   "tst_cube_sphere_multi_var_dims_in_file_time_data"
#define NAME_STATIC_DATA_VARN "distance"
#define NAME_TIME_DATA_VARN_v   "V"
#define NAME_TIME_DATA_VARN_t   "T"
#define NAME_TIME_VARN        "time"

/* A structure for the time variable */
struct time_struct{
  int nVars;
  int nCells;
  int nTimes;
  int nTimesPerFile;
  int nGrids;

  float value;
  char *varname;
  char *st_name;
  char *units;
  char *time_units;
  char *time_stanname;
  char *time_longname;
} TIME;


void create_grids(const char *coordinates_id,
      int nCells, int nGrids,
      const int *faceVectors[], char **gridFiles,
      int gridids[], char **gridNames, int coordids[]) {

  int status = NC_NOERR;
  int i, iGrid;
  int dims[] = {nCells, nCells};
  double lon[dims[0] * dims[1]];
  double lat[dims[0] * dims[1]];
  const int save = 1;
  const char *dimNames[] = {"nlon", "nlat"};

  int globalids[nGrids];

  /* tile 0 */
  iGrid = 0;
  status = nccf_get_cubedsphere_grid(dims, faceVectors[0], lon, lat);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lon_coord(2, dims, dimNames, lon, save, &coordids[0 + 0*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lat_coord(2, dims, dimNames, lat, save, &coordids[1 + 0*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_grid(&coordids[0], NAME_GRID_0, &gridids[0]);
  if (status != NC_NOERR) ERR;
  if(( status = nccf_def_global( &globalids[iGrid] ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_GRIDNAME, NAME_GRID_0, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* tile 1 */
  iGrid = 1;
  status = nccf_get_cubedsphere_grid(dims, faceVectors[0], lon, lat);
  status = nccf_get_cubedsphere_grid(dims, faceVectors[1], lon, lat);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lon_coord(2, dims, dimNames, lon, save, &coordids[0 + 1*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lat_coord(2, dims, dimNames, lat, save, &coordids[1 + 1*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_grid(&coordids[1*2], NAME_GRID_1, &gridids[1]);
  if (status != NC_NOERR) ERR;
  if(( status = nccf_def_global( &globalids[iGrid] ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_GRIDNAME, NAME_GRID_1, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* tile 2 */
  iGrid = 2;
  status = nccf_get_cubedsphere_grid(dims, faceVectors[2], lon, lat);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lon_coord(2, dims, dimNames, lon, save, &coordids[0 + 2*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_lat_coord(2, dims, dimNames, lat, save, &coordids[1 + 2*2]);
  if (status != NC_NOERR) ERR;
  status = nccf_def_grid(&coordids[2*2], NAME_GRID_2, &gridids[2]);
  if (status != NC_NOERR) ERR;
  if(( status = nccf_def_global( &globalids[iGrid] ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_GRIDNAME, NAME_GRID_2, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalids[iGrid], CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;

  /* Create some file names */
  strncpy(gridFiles[0], NAME_GRID_0_FILE, STRING_SIZE);
  strncpy(gridFiles[1], NAME_GRID_1_FILE, STRING_SIZE);
  strncpy(gridFiles[2], NAME_GRID_2_FILE, STRING_SIZE);
  strncpy(gridNames[0], NAME_GRID_0, STRING_SIZE);
  strncpy(gridNames[1], NAME_GRID_1, STRING_SIZE);
  strncpy(gridNames[2], NAME_GRID_2_FILE, STRING_SIZE);

  /* write files */
  for (i = 0; i < nGrids; ++i) {
    int ncid;
    status = nc_create(gridFiles[i], NC_CLOBBER, &ncid);
    if (status != NC_NOERR) ERR;
    status = nccf_put_grid(gridids[i], ncid);
    status = nccf_put_global(globalids[i], ncid);
    
    /* Add global attributes */
    nc_redef( ncid );
    nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE, 
            strlen("Dummy data for GRIDSPEC")+1, 
            "Dummy data for GRIDSPEC" );
    nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation", 
                 "Multiple grids, variables, time files and times per file", 
                 "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
    nc_enddef( ncid );

    status = nc_close(ncid);
    if (status != NC_NOERR) ERR;
    status = nccf_free_global(globalids[i]);
  }
}

int create_stat_data( const char *data_id, const char *coordinates_id,  
                      int nCells, int nStatTiles, int nGrids, int gridids[],
                      char **staticFiles, int staticids[] ){

  int iStat, i, j, index, status, ncid, globalid;
  int save = 0;
  int toterr = 0;
  float  data_f[nCells * nCells];
  const char *df_sn = "static_data";
  const char *df_ut = "meters";
  char *name;
  name = ( char* )calloc( STRING_SIZE, sizeof(char));

  for( iStat = 0; iStat < nStatTiles; iStat++){

    /* Create the data */
    for( j = 0; j < nCells; j++ ){
      for( i = 0; i < nCells; i++ ){
        index = i + j * nCells;
        data_f[index] = iStat + i * (float)( j+1 );
      }
    }

    /* Define the data */
    int gid = iStat % nGrids;
    status += nccf_def_data( gridids[gid], NAME_STATIC_DATA_VARN,
                                       df_sn, df_ut,
                                       NULL, &staticids[iStat] );
    status = nccf_inq_grid_name(gridids[gid], name);
    toterr += status;
    status += nccf_set_data_float( staticids[iStat], data_f, 
					  save, NC_FILL_FLOAT );
    toterr += status;

    if(( status = nccf_def_global( &globalid ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_FILETYPE, 
                                       CF_GLATT_FILETYPE_STATIC_DATA, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_GRIDNAME, 
                                       name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_DATA_ID, 
                                       data_id, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID, 
                                       coordinates_id, 0 ))) ERR;

    /* Open the file */
    sprintf( staticFiles[iStat], "%s%d.nc", NAME_STATIC_DATA, iStat );
    status += nc_create( staticFiles[iStat], NC_CLOBBER, &ncid );

    /* Add global attributes */
    nc_redef( ncid );
    nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE, 
            strlen("Dummy data for GRIDSPEC")+1, 
            "Dummy data for GRIDSPEC" );
    nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation", 
                 "Multiple grids, variables, time files and times per file", 
                 "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
    nc_enddef( ncid );

    /* Write the data to disk */
    status += nccf_put_data(staticids[iStat], ncid);
    status += nccf_put_global(globalid, ncid);
    toterr += status;


    status += nc_close( ncid );
    status += nccf_free_global( globalid );
  }

  free( name );
  return toterr;

}

int create_time_data( const char *data_id, const char *coordinates_id, 
        int iVar, struct time_struct *time, int gridids[], char **timeFiles,
        int timeids[] ){

  int iTile, iTPF, i, j, ind, index, ii, status, ncid, iTime, globalid;
  int time_dimidp, time_varidp;
  int save = 0;
  int toterr = 0;
  long *times;
  float  data_f[time->nCells * time->nCells];
  char *name;
  name = ( char* )calloc( STRING_SIZE, sizeof(char));

  /* Initialize the time data */
  for( j = 0; j < time->nCells; j++ ){
    for( i = 0; i < time->nCells; i++ ){
      ind = i + j * time->nCells;
      data_f[ind] = 0.0;
    }
  }

  float div = time->nGrids + time->nTimes + time->nTimesPerFile + 
              time->nCells + time->nCells;
  div = div - 5;
  float div2  =  div / 2.0;

  for( iTime = 0; iTime < time->nTimes; iTime++ ){
    for( iTile = 0; iTile < time->nGrids; iTile++){

      /* Create the data */
      index = iTile + ( time->nGrids * ( iTime + ( time->nTimes * iVar )));

      /* Open the file*/
      sprintf( timeFiles[index], "%s_%d-%d-%d.nc", NAME_TIME_DATA, 
               iVar, iTime, iTile );
      status = nc_create( timeFiles[index], NC_CLOBBER, &ncid );
      if( status ) return status;

      /* Add the time coordinate */
      nc_redef( ncid );
      status = nccf_def_time(ncid, 
                    "time", 
                    time->nTimesPerFile, 
                    NC_LONG, 
                    time->time_units,
                    time->time_stanname, 
                    &time_dimidp, 
                    &time_varidp);
      nc_enddef( ncid );

      times = malloc( time->nTimesPerFile * sizeof(long));
      for( ii = 0; ii < time->nTimesPerFile; ii++ )
        times[ii] = ii + time->nTimesPerFile * iTime;
      if( status ) ERR;
      status = nc_put_var_long( ncid, time_varidp, times );
      if( status ) ERR;

      /* Define the time data for each tile */
      status = nccf_def_data( gridids[iTile], time->varname, 
                                         time->st_name, time->units,
                                         NAME_TIME_VARN, &timeids[index] );
      toterr += status;

      status = nccf_set_data_float( timeids[index], data_f, 
					   save, NC_FILL_FLOAT );
      toterr += status;

      /* Add global attributes */
      nc_redef( ncid );
      nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE, 
              strlen("Dummy data for GRIDSPEC")+1, 
              "Dummy data for GRIDSPEC" );
      nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation", 
                   "Multiple grids, variables, time files and times per file", 
                   "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
      nc_enddef( ncid );


      /* populate the data */
      for( iTPF = 0; iTPF < time->nTimesPerFile; iTPF++ ){
        for( j = 0; j < time->nCells; j++ ){
          for( i = 0; i < time->nCells; i++ ){
            ind = i + j * time->nCells;
            data_f[ind] = time->value + (( iTile + iTime + iTPF + j + i ) - div2);
          }
        }

        /* Write to disk*/
        status = nccf_put_data(timeids[index], ncid);
        toterr += status;
      }

      /* Global Attributes */
      status = nccf_inq_grid_name(gridids[iTile], name);
      status = nccf_def_global( &globalid );
      status = nccf_add_global_att(globalid, CF_GRIDNAME, name, 0);
      status = nccf_add_global_att(globalid, CF_FILETYPE, 
                                   CF_GLATT_FILETYPE_TIME_DATA, 0);
      status = nccf_add_global_att(globalid, CF_DATA_ID, data_id, 0);
      status = nccf_put_global(globalid, ncid);

      toterr += status;
      status = nc_close( ncid );
      status = nccf_free_global( globalid );
      if( status ) return status;
      free( times );
    }
  }

  free( name );
  return toterr;

}

int main(  ){

  int status = NC_NOERR;
  int i, j, ndims;

  /* create the coordinate id */
  char coordinates_id[36+1];
  char data_id[36+1];
  const int seed = 12345;
  if ((status += nccf_generate_id(seed, coordinates_id))) ERR;
  if ((status += nccf_generate_id(seed, data_id))) ERR;

  /* Dimensions */
  const int nGrids = 3;
  int nvars = 2;
  int nTimes = 2;
  int nTimesPerFile = 3;
  int nTimeFiles = nGrids * nTimes * nvars;
  int nStatFiles = nGrids;

  const int nCells = 10;

  struct time_struct local;

  local.nVars   = nvars;
  local.nCells  = nCells;
  local.nTimes  = nTimes;
  local.nTimesPerFile = nTimesPerFile;
  local.nGrids  = nGrids;
  local.value   = 273.15;
  local.varname = NAME_TIME_DATA_VARN_t;
  local.st_name = "temperature";
  local.units   = "K";
  local.time_stanname = "Time";
  local.time_longname = "Time";
  local.time_units    = "Days";

  /* create the grids */
  const int faceVect0[] = {1, 0, 0};
  const int faceVect1[] = {0, 1, 0};
  const int faceVect2[] = {0, 0, 1};
  const int *faceVectors[] = {faceVect0, faceVect1, faceVect2};
  char **gridFiles   = calloc(nGrids, sizeof(char *));
  char **gridNames   = calloc(nGrids, sizeof(char *));
  char **staticFiles = calloc(nStatFiles, sizeof(char *));
  char **timeFiles   = calloc(nTimeFiles, sizeof(char *));
  for (i = 0; i < nGrids; ++i) {
    gridFiles[i]   = calloc(STRING_SIZE, sizeof(char));
    gridNames[i]   = calloc(STRING_SIZE, sizeof(char));
  }
  for( i = 0; i < nStatFiles; i++ ) staticFiles[i] = calloc(STRING_SIZE, sizeof(char));
  for( i = 0; i < nTimeFiles; i++ ) timeFiles[i]   = calloc(STRING_SIZE, sizeof(char));
  int gridids[nGrids];
  int staticids[nStatFiles];
  int timeids[nTimeFiles];
  int coordids[nGrids * 2];
  create_grids(coordinates_id, nCells, nGrids, faceVectors,
         gridFiles, gridids, gridNames, coordids);

  /* create the mosaic */
  int mosaicid, globalid;
  const double periods[] = {0., 0., 0.};
  status += nccf_def_mosaic(nGrids, gridids, NAME_MOSAIC, &mosaicid);
  if (status != NC_NOERR) ERR;
  status += nccf_compute_mosaic_contacts( mosaicid, periods );
  status += nccf_def_global( &globalid );
  status += nccf_add_global_att(globalid, CF_COORDINATES_ID, coordinates_id, 0);
  status += nccf_add_global_att(globalid, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0);
  if (status != NC_NOERR) ERR;
  int ncid;
  status += nc_create(NAME_MOSAIC_FILE, NC_CLOBBER, &ncid);
  if (status != NC_NOERR) ERR;
  status += nccf_put_mosaic(mosaicid, ncid);
  status += nccf_put_global(globalid, ncid);
  status += nc_close(ncid);
  if (status != NC_NOERR) ERR;
  status += nccf_free_global(globalid);

  /* create static data */
  status += create_stat_data( data_id, coordinates_id, nCells, nStatFiles, nGrids,
                                            gridids, staticFiles, staticids );

 /* create time dependent data */

  status = create_time_data( data_id, coordinates_id, 0, &local, gridids, timeFiles, timeids);
  if( status ) ERR;
  local.value   = 10.0;
  local.varname = NAME_TIME_DATA_VARN_v; 
  local.st_name = "velocity";
  local.units = "m/s";
  status = create_time_data( data_id, coordinates_id, 1, &local, gridids, timeFiles, timeids);
  if( status ) ERR;

  /* create the host file */
  int hostid;
  if (( status = nccf_def_host(coordinates_id, data_id, nTimes, &hostid))) ERR;
  if (( status += nc_create( NAME_HOST_FILE, NC_CLOBBER, &ncid ))) ERR;
  for( i = 0; i < nGrids; i++ ){
    if (( status += nccf_add_host_file( hostid, gridFiles[i], 0 ))) ERR;
  }
  for( i = 0; i < nStatFiles; i++ ){
    if (( status += nccf_add_host_file( hostid, staticFiles[i], 0 ))) ERR;
  }

  /* Reverse the order to test nccf_li_insert */
  for( i = nTimeFiles-1; i >=  0; i-- ){
    if (( status += nccf_add_host_file( hostid, timeFiles[i], 0 ))) ERR;
  }
  if (( status += nccf_add_host_file( hostid, NAME_MOSAIC_FILE, 0 ))) ERR;

  /* Define the global attributes */
  if(( status = nccf_def_global( &globalid ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_FILETYPE, 
                                     CF_GLATT_FILETYPE_HOST, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_DATA_ID, 
                                     data_id , 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID, 
                                     coordinates_id, 0 ))) ERR;

  if (( status += nccf_put_host(hostid, ncid))) ERR;
  if (( status += nccf_put_global(globalid, ncid))) ERR;
  if (( status += nc_close( ncid ))) ERR;
  if (( status += nccf_free_global( globalid ))) ERR;

  /* clean up */
  if ((status += nccf_free_host(hostid))) ERR;
  if ((status += nccf_free_mosaic(mosaicid))) ERR;

  nccf_inq_grid_ndims( gridids[0], &ndims );
  int cid[ndims];
  for (i = 0; i < nGrids; ++i) {
    status += nccf_inq_grid_coordids( gridids[i], cid );

    for( j = 0; j < ndims; j++ )
      if ((status += nccf_free_coord(cid[j]))) ERR;
    if ((status += nccf_free_grid(gridids[i]))) ERR;
  }

  for( i = 0; i < nStatFiles; i++ ){
    if ((status += nccf_free_data(staticids[i]))) ERR;
    free( staticFiles[i] );
  }
  for( i = 0; i < nTimeFiles; i++ ){
    if ((status += nccf_free_data(timeids[i]))) ERR;
    free( timeFiles[i] );
  }
  for (i = 0; i < 2*nGrids; ++i) {
  }

  for (i = 0; i < nGrids; ++i) {
    free(gridFiles[i]);
    free(gridNames[i]);
  }

  free(gridFiles);
  free(gridNames);
  free(staticFiles);
  free(timeFiles);

  return status;

}
