/**
 * Test the boundary indexing in 1D and 2D
 *
 * $Id: tst_bound1.c 719 2011-04-26 17:39:51Z srinath22 $
 * 
 */

#include "nccf_coord.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <nccf_handle_error.h>

int main(int argc, char *argv[]) {

  const int save = 0;
  const int ndimsMax = 3;
  const int dims[] = {10, 11, 12};
  double *data = NULL; // no data
  int coordids[ndimsMax];
  int norm_vec[ndimsMax];
  int status, ndims;
  int is[ndimsMax], ie[ndimsMax];
  char slice[STRING_SIZE];

  const char *dimname1[] = {"ni"};
  const char *dimname2[] = {"ni", "nj", "nk"};
  const int notflipped = 0;

  printf("Running %s\n", argv[0]);

  /* 1d */
  ndims = 1;
  if ((status = nccf_def_lon_coord(ndims, dims, dimname1, data, save, 
					&coordids[ndims - 1]))) ERR;

  norm_vec[0] = -1; // low end
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], norm_vec,
					     is, ie))) ERR;
  assert(is[0] == 0);
  assert(ie[0] == 0);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,  
						   "C", slice))) ERR;
  assert(!strcmp(slice, "0:0"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped, 
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "1:1"));


  norm_vec[0] = +1; // high end
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], norm_vec,
					     is, ie))) ERR;
  assert(is[0] == dims[0] - 1);
  assert(ie[0] == dims[0] - 1);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,  
						   "C", slice))) ERR;
  assert(!strcmp(slice, "9:9"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped, 
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "10:10"));
  

  /* 2d */
  ndims = 2;
  if ((status = nccf_def_lon_coord(ndims, dims, dimname2, data, save, 
					&coordids[ndims - 1]))) ERR;

  norm_vec[0] = -1; norm_vec[1] = 0;
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], 
					     norm_vec,
					     is, ie))) ERR;
  assert(is[0] == 0);
  assert(is[1] == 0);
  assert(ie[0] == 0);
  assert(ie[1] == dims[1] - 1);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,  
						   "C", slice))) ERR;
  assert(!strcmp(slice, "0:0 0:10"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped, 
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "1:11 1:1"));

  norm_vec[0] = +1; norm_vec[1] = 0;
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], 
					     norm_vec,
					     is, ie))) ERR;
  assert(is[0] == dims[0] - 1);
  assert(is[1] == 0);
  assert(ie[0] == dims[0] - 1);
  assert(ie[1] == dims[1] - 1);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped, 
						   "C", slice))) ERR;
  assert(!strcmp(slice, "9:9 0:10"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped, 
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "1:11 10:10"));
  
  norm_vec[0] = 0; norm_vec[1] = -1;
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], 
					     norm_vec,
					     is, ie))) ERR;
  assert(is[0] == 0);
  assert(is[1] == 0);
  assert(ie[0] == dims[0] - 1);
  assert(ie[1] == 0);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,
						   "C", slice))) ERR;
  assert(!strcmp(slice, "0:9 0:0"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "1:1 1:10"));

  norm_vec[0] = 0; norm_vec[1] = +1;
  if ((status = nccf_inq_coord_bound(coordids[ndims - 1], 
					     norm_vec,
					     is, ie))) ERR;
  assert(is[0] == 0);
  assert(is[1] == dims[1] - 1);
  assert(ie[0] == dims[0] - 1);
  assert(ie[1] == dims[1] - 1);
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,
						   "C", slice))) ERR;
  assert(!strcmp(slice, "0:9 10:10"));
  if ((status = nccf_inq_coord_bound_slice(coordids[ndims - 1], 
						   norm_vec, notflipped,
						   "Fortran", slice))) ERR;
  assert(!strcmp(slice, "11:11 1:10"));


  if ((status = nccf_free_coord(coordids[1]))) ERR;
  if ((status = nccf_free_coord(coordids[0]))) ERR;
  
  return 0;
}
