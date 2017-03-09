/*
 * Test the creation of an 8 digit alpha-numeric id.
 *
 * "$Id: tst_nccf_generate_id.c 828 2011-09-14 20:05:08Z pletzer $"
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <netcdf.h>

#include "cf_config.h"
#include "nccf_constants.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

int main(  ){

  int status, compare_value;
  char *id1, *id2, *id3, *id4, *id5;

  id1 = ( char * )calloc( STRING_SIZE, sizeof( char ));
  id2 = ( char * )calloc( STRING_SIZE, sizeof( char ));
  id3 = ( char * )calloc( STRING_SIZE, sizeof( char ));
  id4 = ( char * )calloc( STRING_SIZE, sizeof( char ));
  id5 = ( char * )calloc( STRING_SIZE, sizeof( char ));

  if(( status = nccf_generate_id( 88, id1))) ERR;
//  printf( "The id is:\n%s\n", id );

  if(( status = nccf_generate_id( 7, id2 ))) ERR;
//  printf( "The id is:\n%s\n", id );

  if(( status = nccf_generate_id( 7, id3 ))) ERR;

  if(( status = nccf_generate_id( -6, id4 ))) ERR;
//  printf( "The id is:\n%s\n", id4);

  if(( status = nccf_generate_id( 6, id5 ))) ERR;
//  printf( "The id is:\n%s\n", id5);

/* Compare each id. They should all be different */
  compare_value = strcmp( id1, id2 );
  assert( compare_value != 0 );
  compare_value = strcmp( id1, id3 );
  assert( compare_value != 0 );
  compare_value = strcmp( id1, id4 );
  assert( compare_value != 0 );
  compare_value = strcmp( id1, id5 );
  assert( compare_value != 0 );
  compare_value = strcmp( id2, id3 );
#ifdef HAVE_UUID_H
  assert( compare_value != 0 );
#else
  assert( compare_value == 0 );
#endif
  compare_value = strcmp( id2, id4 );
  assert( compare_value != 0 );
  compare_value = strcmp( id3, id4 );
  assert( compare_value != 0 );
  compare_value = strcmp( id3, id5 );
  assert( compare_value != 0 );
  compare_value = strcmp( id4, id5 );

  free( id1 );
  free( id2 );
  free( id3 );
  free( id4 );
  free( id5 );

  return NC_NOERR;

}
