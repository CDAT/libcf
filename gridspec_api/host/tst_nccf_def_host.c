/**
 * Test the API for defining a host file.
 *
 * "$Id: tst_nccf_def_host.c 767 2011-06-06 23:20:19Z pletzer $"
 */

#include "nccf_host.h"
#include <assert.h>
#include <netcdf.h>
#include "nccf_global.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

int main(  ){

  char coordinates_id[36+1];
  char data_id[36+1];
  const int seed = 12345;
  int hostid, globalid;
  int status;
  if ((status = nccf_generate_id(seed, coordinates_id))) ERR;
  if ((status = nccf_generate_id(seed, data_id))) ERR;
  if ((status = nccf_def_host(coordinates_id, data_id, 0, &hostid))) ERR;
  int ncid;
  /* Define the global attributes */
  if(( status = nccf_def_global( &globalid ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_FILETYPE, 
                                     CF_GLATT_FILETYPE_HOST, 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_DATA_ID, 
                                     data_id , 0 ))) ERR;
  if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID, 
                                     coordinates_id, 0 ))) ERR;

  if (( status = nc_create( "tst_nccf_def_host.nc", NC_CLOBBER, &ncid )));
  if (( status = nccf_put_host(hostid, ncid))) ERR;
  if (( status = nccf_put_global(globalid, ncid))) ERR;
  if (( status = nc_close( ncid ))) ERR;
  if (( status = nccf_free_global( globalid ))) ERR;

  if (( status = nccf_free_host(hostid ))) ERR;

  return NC_NOERR;
}
