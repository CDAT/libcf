/*
 * Test for setting contacts
 *
 * $Id: tst_nccf_set_mosaic_contact.c 783 2011-07-11 15:31:22Z dkindig $
 */

#include <netcdf.h>
#include <stdio.h>
#include "nccf_mosaic.h"


int main(  ){

  int status = NC_NOERR, ngrids = 2, ndims = 2, ncid;
  int gridids[ngrids], mosaicid;
  char *coordnames[] = {"lon", "lat"};
  char mosaicname[] = "set_test";
  char *gridfile0 = "tst_two_tiles_grid0.nc";
  char *gridfile1 = "tst_two_tiles_grid1.nc";
  char *mosaicfile = "test_set_contact_mosaic.nc";
  int coordids[ngrids], i, j;

  /* Create/Define a pair of grids */
  status = nccf_def_grid_from_file( gridfile0, ndims, 
         (const char **)coordnames, "", &gridids[0] );
  if( status ) {
    printf( "file doesn't exist %s\n", gridfile0 );
    return status;
  }
  status = nccf_def_grid_from_file( gridfile1, ndims,
         (const char **)coordnames, "", &gridids[1] );
  if( status ) {
    printf( "file doesn't exist %s\n", gridfile1 );
    return status;
  }
  /* Define a mosaic */
  nccf_def_mosaic( ngrids, gridids, mosaicname, &mosaicid );

  /* Define the contacts manually */
  int grid0_beg_ind[] = {0, 9}, grid0_end_ind[] = {9, 9};
  int grid1_beg_ind[] = {0, 0}, grid1_end_ind[] = {9, 0};

  /* Populate the mosaic */
  nccf_set_mosaic_contact( mosaicid, ndims, 
       gridids[0], gridids[1],
       grid0_beg_ind, grid0_end_ind, 
       grid1_beg_ind, grid1_end_ind);

  /* Write to a file */
  nc_create( mosaicfile, NC_CLOBBER, &ncid );
  nccf_put_mosaic(mosaicid, ncid);
  nc_close( ncid );

  /* Clean up */
  nccf_free_mosaic( mosaicid );
  for( i = 0; i < ngrids; i++ ){
    if(( status = nccf_inq_grid_coordids( gridids[i], coordids ))) ERR;
    if(( status = nccf_free_grid( gridids[i] ))) ERR;
    for( j = 0; j < ndims; j++ )
      if(( status = nccf_free_coord( coordids[j] ))) ERR;
  }

  if( status ) return status;

  /* Verify the contents of the mosaic file */

  return status;

}
