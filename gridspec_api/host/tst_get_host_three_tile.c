/*
 * Test the opening of a host file and retriveing the data from the files
 * contained therein.
 *
 * "$Id: tst_get_host_three_tile.c 824 2011-09-13 18:38:07Z dkindig $"
 */

#include "nccf_host.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <netcdf.h>

#define GRIDPREF "tst_three_tile_cubed_sphere_grid"

int main(  ){

  int status = NC_NOERR;
  const char* host_filename = "tst_three_tile_cubed_sphere_host.nc";
  int hostid, ngrids, nstatdatafiles, ntimedata, nTimes;
  int nStatFiles = 0, nTimeFiles = 0;
  int mosaicid;
  int ndims, i;
  int gridix, statvix, timevix, index, itime;
  char staticfilename[STRING_SIZE];
  char timefilename[STRING_SIZE];
  char *name = NULL;
  char gridfilename[STRING_SIZE];

  status = nccf_def_host_from_file( host_filename, &hostid );
  if (status != NC_NOERR) {
    printf("Non-existant or corrupt file %s! Bailing out\n", host_filename);
    return 0;
  }
  if ((status = nccf_inq_host_ngrids( hostid, &ngrids ))) ERR;
  if ((status = nccf_inq_host_nstatdatafiles( hostid, &nstatdatafiles ))) ERR;
  if ((status = nccf_inq_host_ntimedatafiles( hostid, &ntimedata ))) ERR;
  if ((status = nccf_inq_host_ntimeslices( hostid, &nTimes ))) ERR;

  nStatFiles = nstatdatafiles * ngrids;
  nTimeFiles = ntimedata * ngrids * nTimes;

  int staticdataids[nStatFiles];
  int timedataids[nTimeFiles];

  char mosaicfilename[STRING_SIZE];
  if ((status = nccf_inq_host_mosaicfilename( hostid, mosaicfilename ))) ERR;
  nccf_def_mosaic_from_file( mosaicfilename, "", &mosaicid );
  int gridids[ngrids];
  nccf_inq_mosaic_ndims( mosaicid, &ndims );
  char **coordnames = calloc( ndims, STRING_SIZE * sizeof(char*));
  for( i = 0; i < ndims; i++ ) 
      coordnames[i] = calloc( STRING_SIZE, sizeof(char));
  nccf_inq_mosaic_coordnames( mosaicid, coordnames );

  for( gridix = 0; gridix < ngrids; ++gridix ){
    name = (char*)malloc( STRING_SIZE * sizeof(char));
    sprintf( name, "%s%d", GRIDPREF, gridix );

    if((status = nccf_inq_host_gridfilename( hostid, gridix, gridfilename ))) ERR;
    if((status = nccf_def_grid_from_file( gridfilename, ndims, 
        (const char**)coordnames, name, &gridids[gridix] ))) ERR;
    free( name );

    for( statvix = 0; statvix < nstatdatafiles; ++statvix ){
      
      if((status = nccf_inq_host_statfilename( hostid, statvix, gridix, 
                                                staticfilename ))) ERR;
      index = gridix + ngrids * statvix;
     
      if((status = nccf_def_data_from_file( staticfilename, gridix, 
        "distance", 1, &staticdataids[index]))) ERR;
      
    }
 
    for( itime = 0; itime < nTimes; ++itime ){
      for( timevix = 0; timevix < ntimedata; ++timevix ){
        
        if((status = nccf_inq_host_timefilename( hostid, itime, timevix, gridix,
                                                 timefilename ))) ERR;
        index = gridix + ngrids * timevix;
        if((status = nccf_def_data_from_file( timefilename, gridix, 
          "v", 1, &timedataids[index]))) ERR;
      }
    }
  }

  /* Assertions */
//  assert( ngrids == 3 );
//  assert( nstatdatafiles == 1 );
//  assert( ntimedata == 1 );

  /* Clean up */
  for( i = 0; i < nTimeFiles; i++ ) 
    nccf_free_data( timedataids[i] );
  for( i = 0; i < nStatFiles; i++ ) 
    nccf_free_data( staticdataids[i] );
  int cid[ndims], j;
  for( i = 0; i < ngrids; i++ ){
    nccf_inq_grid_coordids( gridids[i], cid );
    for( j = 0; j < ndims; j++) nccf_free_coord( cid[j] );
    nccf_free_grid( gridids[i] );
  }

  nccf_free_mosaic( mosaicid );
  nccf_free_host( hostid );

  for( i = 0; i < ndims; i++ ) free( coordnames[i] );
  free( coordnames );

  return 0;
}
