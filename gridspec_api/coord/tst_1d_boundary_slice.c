/**
 * Test the boundary indexing in 1D
 *
 * $Id: tst_1d_boundary_slice.c 719 2011-04-26 17:39:51Z srinath22 $
 * 
 */

#include "nccf_coord.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <nccf_handle_error.h>
#include <nccf_varObj.h>

int main() {
  
  double l0[3], l1[3];
  int *dim;
  int ngrids = 2, ndim = 1;
  int coordid[ngrids];
  int *normVec0, *normVec1;
  int i;
  int save = 1, status;
  char slice0[STRING_SIZE], slice1[STRING_SIZE];
  char *gs_slice_format = "F";
  const char *dimnames[] = {"ni"};

  dim = ( int * )malloc( ndim * sizeof( int ));
  normVec0 = ( int * )malloc( ndim * sizeof( int ));
  normVec1 = ( int * )malloc( ndim * sizeof( int ));
  dim[0] = 3;

  for( i = 0; i < 3; i++ ){
    l0[i] = i;
    l1[i] = i+2;
  }

/* Create a coordinate for grid0 */
  nccf_def_lon_coord( ndim, dim, dimnames, l0, save, &coordid[0] );

/* Create a coordinate for grid1 */
  nccf_def_lon_coord( ndim, dim, dimnames, l1, save, &coordid[1] );

  normVec0[0] =  1;
  normVec1[0] = -1;
  
  status = nccf_inq_coord_bound_slice(coordid[0],
                    normVec0, 0,
                    gs_slice_format, slice0);
  status = nccf_inq_coord_bound_slice(coordid[1],
                    normVec1, 0,
                    gs_slice_format, slice1);

  nccf_free_coord( coordid[0] );
  nccf_free_coord( coordid[1] );
  free( dim );
  free( normVec0 );
  free( normVec1 );
  return 0;
}
