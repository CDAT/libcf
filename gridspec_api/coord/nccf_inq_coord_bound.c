/**
 * $Id: nccf_inq_coord_bound.c 737 2011-05-06 17:48:14Z edhartnett $
 * 
 */

#include "nccf_coord.h"
#include <stdlib.h>
#include <nccf_utility_functions.h>

/**
 * \ingroup gs_coord_grp
 * Get boundary contact ranges.
 *
 * \param coordid coordinate ID
 * \param norm_vect an array of -1, 0, 1 elements uniquely denoting the boundary (e.g. 1, 0 for north)
 * \param start_indices (output) start set of indices in the coordinate data
 * \param end_indices (output) inclusive end set of indices in the coordinate data
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
int nccf_inq_coord_bound(int coordid, const int norm_vect[], 
				 int *start_indices, int *end_indices ) {

  int ndims, status;
  int *dims;

  status = nccf_inq_coord_ndims(coordid, &ndims);
  dims = (int *) malloc(ndims * sizeof(int));
  status = nccf_inq_coord_dims(coordid, dims);

  /* 
     Compute the starting/ending positions in index space
     The endingIndices are inclusive (is,...ie)
   */
  nccf_get_start_end_bound_indices(ndims, dims, norm_vect, 0,
				      start_indices, end_indices);
  
  free(dims);

  return NC_NOERR;
}
