/**
 * $Id: nccf_find_indices.h 898 2011-12-22 23:42:03Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

int
nccf_find_indices_TYPE(int ndims, const int dims[], 
		       const _TYPE **coordData, 
                       const _TYPE coord_periodicity[],
		       const _TYPE targetPos[],
		       int *niter, _TYPE *tolpos, 
		       void (*adjustFunc)(int, const int[], _TYPE []),
		       _TYPE dindices[], int hit_bounds[]) {

  int status, i, iter;
  int totError = 0;
  int iterError = 0;
  int boundError = 0;
  int numUpperBoundErrors = 0;
  int numLowerBoundErrors = 0;

  // fix input, if need be
  int nIterMax = abs(*niter);
  _TYPE tolPos = fabs(*tolpos);

  _TYPE dindices_new[ndims];
  _TYPE xPos[ndims];
  _TYPE xNormErr, xNormErrNew;
  _TYPE dIndicesBest[ndims];
  _TYPE xNormErrBest = HUGE_VAL;
  _TYPE alpha = 1.0;

  /* no tolerance should be smaller than this */
  const _TYPE eps = 0.1 * tolPos;
  const _TYPE indexTol = 0.1; // for linesearch

  
  for (i = 0; i < ndims; ++i) {
    dIndicesBest[i] = dindices[i];
  }

#ifdef LOGGING
  printf("-----------\ntargetPos = ");
  for (i = 0; i < ndims; ++i) printf("%f, ", targetPos[i]);
  printf("\n");
  printf("iteration alpha   xNormErr   position/indices/displ \n");
#endif

  iter = -1;
  const int nIterMax2 = (nIterMax/2 > 2? nIterMax/2: 2);
  while (iter < nIterMax && 
         numUpperBoundErrors < nIterMax2 &&
         numLowerBoundErrors < nIterMax2) {

    iter++;

    status = nccf_get_distance_in_coord_TYPE(ndims, dims, 
                                             coordData,
                                             coord_periodicity,
                                             targetPos, 
                                             dindices, 
                                             &xNormErr);
    totError += abs(status);

    // check if done
    if (xNormErr < tolPos) break;

    // compute the increment
    status = nccf_find_next_indices_TYPE(ndims, dims, coordData, 
                                         coord_periodicity,
                                         targetPos,
                                         dindices, dindices_new, xPos);
    totError += abs(status);

#ifdef LOGGING
    printf("%4d      %5.3f  %10.3g ", iter, alpha, xNormErr);
    for (i = 0; i < ndims; ++i) {
      printf("%10.6f, ", xPos[i]);
    }
    for (i = 0; i < ndims; ++i) {
      printf("%10.4f, ", dindices[i]);
    }
    for (i = 0; i < ndims; ++i) {
      printf("%10.4f, ", dindices_new[i] - dindices[i]);
    }
    printf("\n");
#endif

    // adjust indices if these fall outside the domain
    if (adjustFunc) {
      adjustFunc(ndims, dims, dindices);
    }

    // error handling
    boundError = 0;
    for (i = 0; i < ndims; ++i) {

      if (dindices_new[i] < -eps) {
        /* a negative value means we went beyond the left boundary */
        hit_bounds[i] = -1;
        numUpperBoundErrors++;
        boundError++;
      }
      else if (dindices_new[i] > dims[i] - 1 + eps) {
        /* a positive value means we went beyond the right boundary */
        hit_bounds[i] = 1;
        numLowerBoundErrors++;
        boundError++;
      }
      else {
        hit_bounds[i] = 0;
      }
    }

    // correct if overshooting
    for (i = 0; i < ndims; ++i) {
      dindices_new[i] = (dindices_new[i] < 0? 0: dindices_new[i]);
      dindices_new[i] = (dindices_new[i] > dims[i] - 1? dims[i] - 1: dindices_new[i]);
    }

    status = nccf_get_distance_in_coord_TYPE(ndims, dims, 
                                             coordData,
                                             coord_periodicity,
                                             targetPos, 
                                             dindices_new,
                                             &xNormErrNew);
    totError += abs(status);

    // do a linesearch
    if (xNormErrNew > 0.99 * xNormErr) {
      status = nccf_linesearch_indices_TYPE(ndims, dims, coordData,
                                            coord_periodicity,
                                            targetPos, indexTol,
                                            dindices, dindices_new);
      totError += abs(status);
    }
    
    /* if (bad) { */
    /*   // generally means we left the domain */
    /*   break; */
    /* } */

    // update
    for (i = 0; i < ndims; ++i) {
      dindices[i] += alpha*(dindices_new[i] - dindices[i]);
    }

    // are we improving?
    status = nccf_get_distance_in_coord_TYPE(ndims, dims, 
                                             coordData,
                                             coord_periodicity,
                                             targetPos, 
                                             dindices,
                                             &xNormErrNew);
    totError += abs(status);

    // store the best result so far
    if (xNormErrNew < xNormErrBest) {
      for (i = 0; i < ndims; ++i) {
        dIndicesBest[i] = dindices[i];
      }
      xNormErrBest = xNormErrNew;
    }
    
    /*     if (xNormErrNew > 0.99 * xNormErr) { */
    /*       // add some relaxation */
    /*       alpha *= 0.9; */
    /* /\*       for (i = 0; i < ndims; ++i) { *\/ */
    /* /\* 	dindices[i] += eps * (_TYPE)random()/ (_TYPE)RAND_MAX; *\/ */
    /* /\*       } *\/ */
    /*     } */
    /*     else if (xNormErrNew < 0.7 * xNormErr) { */
    /*       alpha = 1.0; */
    /*     } */
  }
    
  // set the error code
  if (iter >= nIterMax || boundError) {
    iterError = 1;
  }

  // make sure we return our best guess so far...
  for (i = 0; i < ndims; ++i) {
    dindices[i] = dIndicesBest[i];
  }

  *niter = iter;
  *tolpos = xNormErrBest;
  
#ifdef LOGGING
  if (iterError || totError || boundError) {
    printf("*** Failure in nccf_find_indices_TYPE: totError=%d iterError=%d ", 
           totError, iterError);
    for (i = 0; i < ndims; ++i) {
      if (hit_bounds[i] > 0) {
        printf(" hit high bound %d", i);
      }
      else if (hit_bounds[i] < 0) {
        printf(" hit low bound %d", i);
      }
    }
    printf("\n");
  }
#endif
  return iterError;
}

