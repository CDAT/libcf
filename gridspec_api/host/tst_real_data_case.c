/**
 * Use the data from modave/data/c48_data to create a mosaic and host file.
 * Then reread the host file.
 *
 * $Id: tst_real_data_case.c 928 2012-05-22 00:03:52Z dkindig $
 */

#include "nccf_host.h"
#include <netcdf.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "nccf_constants.h"
#include "nccf_utility_functions.h"
#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"

#define GRIDNAME "horizontal_grid"
#define DATANAME "19810101.atmos_month"
#define TILE     ".tile"
#define FILETYPE "gs_tile_file"
#define DATATYPE "gs_time_data_file"
#define MOSAICFILE "tst_real_data_case_mosaic.nc"
#define HOSTFILE "tst_real_data_case_host.nc"

int update_attrfiles(const char *name,   
		     const char *file, 
		     const char *coordinates_id, 
		     const char *data_id,
		     const char *filetype,
		     int data){

  int ncid, status;
  /* Update the grid files global attributes */
  status = nc_open(file, NC_WRITE, &ncid);
  switch(status){
    case 0: break;
    case 2: printf("%s\n",file);
    default: ERR;
  }

  if ((status = nc_redef(ncid))) ERR;
  if ((status = nc_put_att_text(ncid, NC_GLOBAL, CF_FILETYPE,
                   strlen(filetype), filetype))) ERR;
  if ((status = nc_put_att_text(ncid, NC_GLOBAL, CF_GRIDNAME,
                   strlen(name), name))) ERR;
  if ((status = nc_put_att_text(ncid, NC_GLOBAL, CF_COORDINATES_ID,
                   strlen(coordinates_id), coordinates_id))) ERR;
  if( data ){
    if (data_id) {
      if ((status = nc_put_att_text(ncid, NC_GLOBAL, CF_DATA_ID,
				   strlen(data_id), data_id))) ERR;
    }
    if (( status = nc_put_att_text( ncid, NC_GLOBAL, "run", 
                                  strlen("run_01"), "run_01" ))) ERR;
    if (( status = nc_put_att_text( ncid, NC_GLOBAL, "model", 
                                  strlen("model_01"), "model_01" ))) ERR;
    if (( status = nc_put_att_text( ncid, NC_GLOBAL, "Institution", 
                                  strlen("GFDL"), "GFDL" ))) ERR;
  }
  if(( status = nc_close( ncid ))) ERR;
  return NC_NOERR;
}

int main(int argc, char **argv){


  char *dirname;

  if (argc < 2) {
    printf("Error\nUsage:\n");
    printf("%s <directory where 19810101.atmos_month data reside>\n", argv[0]);
    exit(1);
  } else {
    dirname = calloc(strlen(argv[1]) + 1, sizeof(char));
    strcpy(dirname, argv[1]);
  }

  int status = NC_NOERR, i, ncid, hostid, globalid;
  char datafile[STRING_SIZE];
  char gridfile[STRING_SIZE];
  char gridname[STRING_SIZE];
  char coordinates_id[STRING_SIZE];
  char data_id[STRING_SIZE];
  char *coordNames[] = {"x", "y"};
  int nGrids = 6, nDims = 2;
  int gridids[nGrids], mosaicid;
  double period[] = {0.0, 0.0};
  int force = 0;  // Don't force a file into the host file
  int seed  = 54321;

  if ((status = nccf_generate_id(seed, coordinates_id))) ERR;
  if ((status = nccf_generate_id(seed, data_id))) ERR;


  for( i = 1; i <=  nGrids; i++ ){
    sprintf( gridfile, "%s%s%s%d%s", dirname, GRIDNAME, TILE, i, ".nc");
    sprintf( gridname, "%s%s%d", GRIDNAME, TILE, i );
    if(( status = update_attrfiles(gridname, gridfile, 
				   coordinates_id, NULL, 
				   CF_GLATT_FILETYPE_GRID, 1))) ERR;
    if((status = nccf_def_grid_from_file( gridfile, nDims, 
                (const char**)coordNames, gridname, &gridids[i-1] ))) ERR;
  }
  for( i = 1; i <=  nGrids; i++ ){
    sprintf( datafile, "%s%s%s%d%s", dirname, DATANAME, TILE, i, ".nc");
    sprintf( gridname, "%s%s%d", GRIDNAME, TILE, i );
    if((status = update_attrfiles(gridname, datafile, 
				  coordinates_id, data_id, 
				  CF_GLATT_FILETYPE_TIME_DATA, 1))) ERR;
  }

  /* Create the mosaic file and add it to the host file */
  nccf_def_mosaic(nGrids, gridids, gridname, &mosaicid);
  status += nccf_compute_mosaic_contacts( mosaicid, period );
  nccf_def_global( &globalid );
  nccf_add_global_att(globalid, CF_COORDINATES_ID, coordinates_id, 0);
  nccf_add_global_att(globalid, CF_FILETYPE, CF_GLATT_FILETYPE_MOSAIC, 0);

  if (( status = nc_create( MOSAICFILE, NC_CLOBBER, &ncid ))) ERR;
  if (( status = nccf_put_mosaic(mosaicid, ncid))) ERR;
  if (( status = nccf_put_global(globalid, ncid))) ERR;
  if (( status = nc_close( ncid ))) ERR;
  nccf_free_global( globalid );

  nccf_def_host( coordinates_id, data_id, 1, &hostid );

  /* Add the mosaic file */
  nccf_add_host_file( hostid, MOSAICFILE, 1);

  /* Add the grid names and files */
  for( i = 1; i <= nGrids; i++  ){
    sprintf( gridfile, "%s%s%s%d%s", dirname, GRIDNAME, TILE, i, ".nc");
    nccf_add_host_file( hostid, gridfile, force );
  }

  /* Add the data files */
  for( i = 1; i <= nGrids; i++  ){
    sprintf( datafile, "%s%s%s%d%s", dirname, DATANAME, TILE, i, ".nc");
    nccf_add_host_file( hostid, datafile, force );
  }

  if (( status = nc_create( HOSTFILE, NC_CLOBBER, &ncid ))) ERR;
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

  nccf_free_host( hostid );
  nccf_free_mosaic( mosaicid );
  int cid[nDims], j;
  for( i = 0; i < nGrids; i++ ){
    nccf_inq_grid_coordids( gridids[i], cid );
    for( j = 0; j < nDims; j++ ) nccf_free_coord( cid[j] );
    nccf_free_grid( gridids[i] );
  }

  /* Read the host file */
  int ngrids, nstatdatafiles, ntimedata, ntimes, hostid2;
  if (( status = nccf_def_host_from_file( HOSTFILE, &hostid2 ))) ERR;

  if (( status = nccf_inq_host_ngrids( hostid2, &ngrids ))) ERR;
  if (( status = nccf_inq_host_nstatdatafiles( hostid2, &nstatdatafiles ))) ERR;
  if (( status = nccf_inq_host_ntimedatafiles( hostid2, &ntimedata ))) ERR;
  if (( status = nccf_inq_host_ntimeslices( hostid2, &ntimes ))) ERR;

  /* Get the mosaic file name */
  char mosaicfilename[STRING_SIZE];
  int mosaicid2;
  nccf_inq_host_mosaicfilename( hostid2, mosaicfilename );
  nccf_def_mosaic_from_file( mosaicfilename, "", &mosaicid2 );
  nccf_inq_mosaic_ndims( mosaicid2, &nDims );
  char **cdNm;
  cdNm = (char**)calloc( STRING_SIZE, 2 * sizeof(char*));
  cdNm[0] = (char*)calloc( STRING_SIZE, sizeof(char));
  cdNm[1] = (char*)calloc( STRING_SIZE, sizeof(char));
  nccf_inq_mosaic_coordnames( mosaicid2, cdNm );

  /* Get the grid file names */
  int gridid2[ngrids], ii;
  if( ngrids != 0 ){
    for( ii = 0; ii < ngrids; ++ii ){
      nccf_inq_host_gridfilename( hostid2, ii, gridfile );
      gridid2[ii] = ii;
    }
  }

  /* Assert the host_from_file is behaving */
  assert( ngrids == 6 );
  assert( nstatdatafiles == 0 );
  assert( ntimedata == 1 );
  assert( ntimes == 1 );

  int tfindx, vfindx, gfindx;
  char fname[STRING_SIZE];
  /* time var file names */
  for (tfindx = 0; tfindx < ntimes; ++tfindx) {
    for (vfindx = 0; vfindx < ntimedata; ++vfindx) {
      for (gfindx = 0; gfindx < ngrids; ++gfindx) {
	      status = nccf_inq_host_timefilename(hostid2, tfindx, vfindx, gfindx,
	      				    fname);
	      if (status) ERR;
	      printf("time data file name [%d][%d][%d]: %s\n", tfindx, vfindx, 
	       gfindx, fname);
      }
    }
  }
  /* static var file names */
  for (vfindx = 0; vfindx < nstatdatafiles; ++vfindx) {
    for (gfindx = 0; gfindx < ngrids; ++gfindx) {
      status = nccf_inq_host_statfilename(hostid2, vfindx, gfindx,
					  fname);
      if (status) ERR;
      printf("stat data file name [%d][%d]: %s\n", vfindx, gfindx, fname);
    }
  }

  /* Clean up */
  nccf_free_mosaic( mosaicid2 );
  nccf_free_host( hostid2 );
  free(dirname);

  return status;

}
