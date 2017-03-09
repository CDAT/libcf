/**
 * $Id: nccf_get_linear_weights.h 891 2011-12-21 17:22:02Z pletzer $
 */

int nccf_get_linear_weights_TYPE(int ndims, const int dims[], 
                                 const _TYPE dindices[],
                                 const int imask[],
                                 _TYPE weights[]) {
  int i, j, k;
  int num_invalid;
  int all_weights_zero = 1;
  int bi, neigh_bi;
  int ier;
  int num_invalid_corners;

  // find the lower corner set of indices,
  // the number of nodes in the cell, and the
  // total number of nodes
  int loCornerIndx[ndims];
  int nNodes = 1;
  int ntot = 1;
  for (i = 0; i < ndims; ++i) {
    loCornerIndx[i] = (int) floor(dindices[i]);
    nNodes *= 2;
    ntot *= dims[i];
  }
  
  // iterate over all nodes with the cell
  int displ[ndims]; // intra cell displacements from loCornerIndx
  int indx[ndims];  // index set
  int indxBar; // 1 - indx
  int neigh_jindx[ndims];  // the node index 0...2^ndims of the neighbor node
  int neigh_indx[ndims][ndims]; // the neighbor index sets
  int valid_nodes[nNodes]; // list of valid nodes in cell
  int offset[ndims]; // offsets from reference node, aka (1, 0, 0)

  _TYPE xis[ndims]; // triangular parametrization
  _TYPE rhs[ndims]; // target position - reference node position
  _TYPE mat[ndims*ndims]; // matrix of (neigh position - ref node position)
  _TYPE sum_xi;

  int prodDims[ndims];
  prodDims[ndims - 1] = 1; 
  for (i = ndims - 2; i >= 0; --i) {
    prodDims[i] = prodDims[i + 1] * 2;
  }

  num_invalid = 0;
  for (j = 0; j < nNodes; ++j) {

    valid_nodes[j] = 1;

    // compute the index set of the node, the index set opposite to 
    // the node, and the quad interp weights
    weights[j] = 1.0;
    for (i = 0; i < ndims; ++i) {
      displ[i] = j / prodDims[i] % 2; // vector of 0s and 1s
      indx[i] = loCornerIndx[i] + displ[i];
      indxBar = loCornerIndx[i] + 1 - displ[i];
      weights[j] *= fabs(dindices[i] - (_TYPE)indxBar);
    }

    if (imask != NULL) {
      bi = nccf_get_flat_index(ndims, dims, indx);
      bi = (bi > ntot - 1? ntot - 1: bi);
      bi = (bi < 0? 0: bi);
      if (imask[bi] == 0) {
	num_invalid++;
	valid_nodes[j] = 0;
      }
    }
  } // end of cell node iteration

  if (imask != NULL && num_invalid !=0) {
    
    for (j = 0; j < nNodes; ++j) {
      // default is no interpolation
      weights[j] = 0;
    } 

    if (nNodes - num_invalid >= ndims + 1) {
      // try tri interpolation 
      for (j = 0; j < nNodes; ++j) {
	if (valid_nodes[j] == 1) {
	  // there is hope that tri interpolation will work. However,
	  // will need to have all the closest neighbors to be also 
	  // valid
	  for (i = 0; i < ndims; ++i) {
	    displ[i] = j / prodDims[i] % 2;
	    indx[i] = loCornerIndx[i] + displ[i];
	  }
	  bi = nccf_get_flat_index(ndims, dims, indx);
	  bi = (bi > ntot - 1? ntot - 1: bi);
	  bi = (bi < 0? 0: bi);
	  num_invalid_corners = 0;
	  for (k = 0; k < ndims; ++k) {
	    // iterate over directions
	    for (i = 0; i < ndims; ++i) {
	      offset[i] = 0;
	      neigh_indx[k][i] = indx[i];
	    }
	    // if displ[k] == 0 then offset[k] = 1
	    // if displ[k] == 1 then offset[k] = -1
	    offset[k] = 1 - 2*displ[k];
	    neigh_jindx[k] = 0;
	    for (i = 0; i < ndims; ++i) {
	      neigh_indx[k][i] += offset[i];
	      neigh_jindx[k] += (displ[i] + offset[i])*prodDims[i];
	    }
	    neigh_bi = nccf_get_flat_index(ndims, dims, neigh_indx[k]);
	    neigh_bi = (neigh_bi > ntot - 1? ntot - 1: neigh_bi);
	    neigh_bi = (neigh_bi < 0? 0: neigh_bi);
	    if (imask[neigh_bi] == 0) {
	      num_invalid_corners++;
	      // no good, cannot interpolate, must have ndims valid neighbors
	      break; // k loop
	    }
	  } // end of k loop
	  if (num_invalid_corners == 0) {
	    // if we reached this point then we may be able to
	    // apply tri interpolation. solve for the xi cell parametrization.
	    // The xis must be within [0, 1] and sum xi <=1.
	    for (i = 0; i < ndims; ++i) {
	      rhs[i] = dindices[i] - indx[i];
	      for (k = 0; k < ndims; ++k) {
		mat[i*ndims + k] = neigh_indx[k][i] - indx[i];
	      }
	    }
	    ier = nccf_solve_TYPE(ndims, mat, rhs, xis);
	    int test_interp = 1;
	    // sum of all xis <= 1 for the target point to be in the triangle
	    sum_xi = 0;
	    for (i = 0; i < ndims; ++i) {
	      sum_xi += xis[i];
	      if (xis[i] < 0 - _EPS || xis[i] > 1 + _EPS) {
		test_interp = 0;
	      }
	    }
	    if (test_interp == 1 && sum_xi < 1 + _EPS) {
	      // now set the weights
	      weights[j] =  1 - sum_xi;
	      for (k = 0; k < ndims; ++k) {
		weights[neigh_jindx[k]] = xis[k];
	      }
	      all_weights_zero = 0; // success
	      break; // finish j loop
	    }
	  }   // num_invalid_corners test
	} // valid_nodes test
      } // j loop
    } // number of invalid nodes test
  } // mask test

  if (num_invalid != 0 && all_weights_zero != 0) {
    // no interpolation was possible (too many invalid data in
    //                                in this cell)
    return -num_invalid; 
  }
  else {
    // 0=all data are valid, >0 triangular interpolation was successful
    return num_invalid;
  }
}
