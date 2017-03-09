/**
 * $Id: nccf_find_next_indices.h 905 2011-12-29 04:56:48Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

int
nccf_find_next_indices_TYPE(int ndims, const int dims[], 
			    const _TYPE **coordData,
                            const _TYPE coord_periodicity[],
                            const _TYPE targetPos[],
                            const _TYPE dIndices_in[],
                            _TYPE dIndices_out[], 
                            _TYPE position_out[]) {

  int status, i, j, k;
  int totError = 0;

  _TYPE jac[ndims * ndims];
  int loCornerIndx[ndims];
  int bigIndex0;
  int indx[ndims];
  _TYPE idxDispl[ndims];
  _TYPE posDispl[ndims];

  status = nccf_get_position_TYPE(ndims, dims, coordData, 
                                  coord_periodicity, targetPos,
                                  dIndices_in, position_out);
  totError += abs(status);

  /* compute the position */
  for (i = 0; i < ndims; ++i) {

    posDispl[i] = targetPos[i] - position_out[i];
    if (coord_periodicity[i] < _HUGE_TYPE) {
      if (fabs((double)(posDispl[i] - coord_periodicity[i])) < 
          fabs((double)posDispl[i])) {
        posDispl[i] -= coord_periodicity[i];
      }
      if (fabs((double)(posDispl[i] + coord_periodicity[i])) < 
          fabs((double)posDispl[i])) {
        posDispl[i] += coord_periodicity[i];
      }
    }

    loCornerIndx[i] = (int) floor(dIndices_in[i]);
    loCornerIndx[i] = (loCornerIndx[i] >= dims[i] - 1? dims[i] - 2: loCornerIndx[i]);
    loCornerIndx[i] = (loCornerIndx[i] < 0? 0: loCornerIndx[i]);

    indx[i] = loCornerIndx[i];
  }
  bigIndex0 = nccf_get_flat_index(ndims, dims, loCornerIndx);

  /* Jacobian */
  int bigIndexPlus;
  for (i = 0; i < ndims; ++i) {
    for (j = 0; j < ndims; ++j) {
      indx[j] += 1;
      bigIndexPlus = nccf_get_flat_index(ndims, dims, indx);
      k = j + ndims*i;
      jac[k] = (coordData[i][bigIndexPlus] - coordData[i][bigIndex0]);
      if (coord_periodicity[i] < _HUGE_TYPE) {
        if (fabs((double)(jac[k] - coord_periodicity[i])) < 
            fabs((double)jac[k])) {
          jac[k] -= coord_periodicity[i];
        }
        if (fabs((double)(jac[k] + coord_periodicity[i])) < 
            fabs((double)jac[k])) {
          jac[k] += coord_periodicity[i];
        }
      }
      indx[j] -= 1;
    }
  }

  /* compute the increment */
  status = nccf_solve_TYPE(ndims, jac, posDispl, idxDispl);
  totError += abs(status);

  /* compute the new index position */
  for (i = 0; i < ndims; ++i) {
    dIndices_out[i] = dIndices_in[i] + idxDispl[i];
  }
  
  return totError;
}
