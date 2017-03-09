/**
 * 
 * A test for the add_host method.
 *
 * "$Id: tst_nccf_add_host.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_host.h"
#include <assert.h>
#include <netcdf.h>
#include "nccf_global.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

int main(  ){

  const char coordinates_id[] = "696c568a-31aa-4cdc-9c1a-c32821cc25ad";
  const char data_id[] = "29fb9258-2458-11e0-8ea4-5c260a1834c1";
  int hostid, globalid;
  int status;
  char *mosaic_filename = "tst_three_tile_cubed_sphere_mosaic.nc";
  char *host_filename = "tst_nccf_add_host.nc";

  if ((status = nccf_def_host(coordinates_id, data_id, 1, &hostid))) ERR;
  status = nccf_add_host_file( hostid, mosaic_filename, 1 );
  int ncid;

  /* Define the global attributes */
  if(( status = nccf_def_global( &globalid ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_FILETYPE, 
                                     CF_GLATT_FILETYPE_HOST, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_DATA_ID, 
                                     data_id , 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID, 
                                     coordinates_id, 0 ))) ERR;

  if ((status = nc_create(host_filename, NC_CLOBBER, &ncid))) ERR;
  if (( status += nccf_put_host(hostid, ncid))) ERR;
  if (( status += nccf_put_global(globalid, ncid))) ERR;
  if (( status += nc_close( ncid ))) ERR;
  if (( status += nccf_free_global( globalid ))) ERR;
  if (( status += nccf_free_host( hostid ))) ERR;

  return status;

}
