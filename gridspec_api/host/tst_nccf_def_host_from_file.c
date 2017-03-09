/*
 * Test the method nccf_def_host_from_file.c  
 * NOTE: This test should be run after tst_three_tile_cubed_sphere_host
 *
 * "$Id: tst_nccf_def_host_from_file.c 719 2011-04-26 17:39:51Z srinath22 $"
 *
 */

#include "nccf_host.h"
#include <string.h>
#include <stdio.h>

#include "nccf_handle_error.h"

int main(  ){

  int status = NC_NOERR;
  int hostid;
  char *filename = "tst_three_tile_cubed_sphere_host.nc";

  // If the file is missing just exit as a pass.
  status = nccf_def_host_from_file( filename, &hostid );
  if ((status == NC_ENOTNC )) status = NC_NOERR;
  else ERR;
  if ((status = nccf_free_host( hostid ))) ERR;

  return status;
}
