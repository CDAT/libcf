/**
 * Test the global attribute methoc.
 *
 * $Id: tst_global.c 767 2011-06-06 23:20:19Z pletzer $
 */

#include "nccf_global.h"
#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "libcf_src.h"
#include "nccf_varObj.h"
#include "nccf_errors.h"

int main( ){

  int status = NC_NOERR, globalId, globalId2, ncid;
  char value[STRING_SIZE];
  char *filename = "tst_global.nc";
  char attname[STRING_SIZE], attval[STRING_SIZE];
  int natts, i;

  status = nccf_def_global( &globalId );
  /* SET a global */
  status = nccf_add_global_att( globalId, "name1", "value1", 0 );
  assert( status == NC_NOERR );
  status = nccf_add_global_att( globalId, "name2", "value2", 0 );
  assert( status == NC_NOERR );
  /* Value other than 0, 1 or 2 */
  status = nccf_add_global_att( globalId, "name3", "value3", -1 );
  assert( status  == NC_NOERR );

  /* APPEND an attribute value to an existing attribute */
  status = nccf_add_global_att( globalId, "name3", "value4", 1 );
  assert( status == NC_NOERR );

  /* REPLACE an attribute value with a new value */
  status = nccf_add_global_att( globalId, "name2", "value5", 2 );
  assert( status == NC_NOERR );

  if(( status = nc_create( filename, NC_CLOBBER, &ncid ))) ERR;
  status = nccf_put_global(globalId, ncid);
  assert( status  == NC_NOERR );
  if(( status = nc_close( ncid ))) ERR;

  status = nccf_free_global( globalId );
  assert( status == NC_NOERR );

  status = nccf_def_global_from_file( filename, &globalId2 );
  assert( status == NC_NOERR );
  
  status = nccf_inq_global_att( globalId2, "name1", value );
  assert( strcmp( value, "value3" ) != 0 );
  assert( strcmp( value, "value1" ) == 0 );

  nccf_inq_global_natts( globalId2, &natts );
  for( i = 0; i < natts; i++ ){
    nccf_inq_global_attval( globalId2, i, attname, attval );
    switch( i ){
      case 0: assert( strcmp( "name1", attname ) == 0 ); 
              assert( strcmp( "value1", attval ) == 0 ); break;
      case 1: assert( strcmp( "name2", attname ) == 0 ); break;
              assert( strcmp( "value5", attval ) == 0 ); break;
      case 2: assert( strcmp( "name3", attname ) == 0 ); break;
              assert( strcmp( "value3", attval ) == 0 ); break;
    }

  }
  status = nccf_free_global( globalId2 );

  return NC_NOERR;
}
