/**
 * Perform line search minimization in index space
 * $Id: nccf_linesearch_indices.h 891 2011-12-21 17:22:02Z pletzer $
 *
 * \author Alexander Pletzer
 */

int 
nccf_linesearch_indices_TYPE(int ndims, const int dims[], 
			const _TYPE **coordData,
                        const _TYPE coord_periodicity[],
			const _TYPE targetPos[],
			_TYPE indexTol, 
			const _TYPE dindices[], 
			_TYPE dindices_new[]) {

  int status, i;
  int totErr = 0;

  _TYPE leftIndices[ndims];
  _TYPE righIndices[ndims];
  _TYPE middIndices[ndims];
  _TYPE leftErr, middErr, righErr, indexErr;

  for (i = 0; i < ndims; ++i) {
    leftIndices[i] = dindices[i];
    righIndices[i] = dindices_new[i];
    middIndices[i] = 0.5 * (leftIndices[i] + righIndices[i]);
  }

  status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                           coord_periodicity,
                                           targetPos, leftIndices, 
                                           &leftErr);
  totErr += abs(status);
  status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                           coord_periodicity,
                                           targetPos, righIndices, 
                                           &righErr);
  totErr += abs(status);
  status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                           coord_periodicity,
                                           targetPos, middIndices, 
                                           &middErr);
  totErr += abs(status);
  
  indexErr = 0.0;
  for (i = 0; i < ndims; ++i) {
    indexErr += (righIndices[i] - leftIndices[i]) * 
      (righIndices[i] - leftIndices[i]);
  }
  indexErr = sqrt(indexErr);

  while (indexErr > indexTol) {

    if (leftErr < righErr) {
      for (i = 0; i < ndims; ++i) {
	righIndices[i] = middIndices[i];
      }
    }
    else {
      for (i = 0; i < ndims; ++i) {
	leftIndices[i] = middIndices[i];
      }
    }

    for (i = 0; i < ndims; ++i) {
      middIndices[i] = 0.5 * (leftIndices[i] + righIndices[i]);
    }
    
    status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                             coord_periodicity,
                                             targetPos, leftIndices, 
                                             &leftErr);
    totErr += abs(status);
    status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                             coord_periodicity,
                                             targetPos, righIndices, 
                                             &righErr);
    totErr += abs(status);
    status = nccf_get_distance_in_coord_TYPE(ndims, dims, coordData,
                                             coord_periodicity,
                                             targetPos, middIndices, 
                                             &middErr);
    totErr += abs(status);
  
    indexErr = 0.0;
    for (i = 0; i < ndims; ++i) {
      indexErr += (righIndices[i] - leftIndices[i]) * 
        (righIndices[i] - leftIndices[i]);
    }
    indexErr = sqrt(indexErr);
  }

  for (i = 0; i < ndims; ++i) {
    dindices_new[i] = middIndices[i];
  }

  return totErr;
}
