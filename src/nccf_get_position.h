/**
 * $Id: nccf_get_position.h 905 2011-12-29 04:56:48Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

int
nccf_get_position_TYPE(int ndims, const int dims[],
                       const _TYPE **coordData,
                       const _TYPE coord_periodicity[],
                       const _TYPE pos_ref[],
                       const _TYPE dindices[],
                       _TYPE position[]) {

  int status, i;
  int totError = 0;
  _TYPE *weights = (_TYPE *) malloc( pow(2, ndims) * sizeof(_TYPE) );

  status = nccf_get_linear_weights_TYPE(ndims, dims, 
					dindices,
					NULL,
					weights);
  totError += abs(status);

  for (i = 0; i < ndims; ++i) {
    if (pos_ref && coord_periodicity[i] < _HUGE_TYPE) {
      status = nccf_linear_interp_mod_TYPE(ndims, dims, 
                                           coordData[i],
                                           pos_ref[i],
                                           coord_periodicity[i], 
                                           dindices, weights, 
                                           _FILL_VALUE, 
                                           &position[i]);
    }
    else {
      status = nccf_linear_interp_TYPE(ndims, dims, 
                                       coordData[i], 
                                       dindices, weights, 
                                       _FILL_VALUE, 
                                       &position[i]);      
    }
    totError += abs(status);
  }

  free(weights);

  return totError;
}
