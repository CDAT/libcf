/**
 * Test linear masked interpolation
 * $Id: tst_nccf_linear_interp_mask.c 891 2011-12-21 17:22:02Z pletzer $
 */

#include <assert.h>
#include <math.h>
#include <cf_config.h>
#include <nccf_constants.h>
#include <nccf_utility_functions.h>

const double tol = 1.e-6;

int main() {
#ifdef HAVE_LAPACK_LIB
  {
    /* 1d tests */
    int status;
    const int ndims = 1;
    const int nnodes = 2;
    const int dims[] = {2};
    double weights[nnodes];

    // all data valid
    const double dindices[] = {0.2};
    int imask[] = {1, 1};
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 0);
    assert(fabs(weights[0] + weights[1] - 1) < tol);
    assert(fabs(weights[0] - 0.8) < tol);
    assert(fabs(weights[1] - 0.2) < tol);
    
    // one invalid datum, not possible to interpolate
    imask[0] = 0; imask[1] = 1;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -1);
    assert(fabs(weights[0] + weights[1]) < tol);

    // no valid data
    imask[0] = 0; imask[1] = 0;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -2);
    assert(fabs(weights[0] + weights[1]) < tol); 
  }
  
  /* 2d tests */
  {
    int status;
    const int ndims = 2;
    const int nnodes = 4;
    const int dims[] = {2, 2};
    double weights[nnodes];

    // all data valid
    double dindices[] = {0.2, 0.6};
    int imask[] = {1, 1, 1, 1};
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 0);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3]- 1) < tol);
    assert(fabs(weights[0] - (1-dindices[0])*(1-dindices[1])) < tol);
    assert(fabs(weights[1] - (1-dindices[0])*(0+dindices[1])) < tol);
    assert(fabs(weights[2] - (0+dindices[0])*(1-dindices[1])) < tol);

    // 1 missing datum, ok to interpolate
    imask[0] = 1; imask[1] = 1; imask[2] = 0; imask[3] = 1;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 1);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3] - 1) < tol);
    assert(fabs(weights[0] - (1-0.6)) < tol);
    assert(fabs(weights[1] - (1 - (1-0.6) - 0.2)) < tol);
    assert(fabs(weights[2] -   0) < tol);
    assert(fabs(weights[3] - 0.2) < tol);
    
    // 1 missing datum, cannot interpolate
    imask[0] = 0; imask[1] = 1; imask[2] = 1; imask[3] = 1;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -1);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3]) < tol);

    // 1 missing datum, ok to interpolate
    imask[0] = 1; imask[1] = 1; imask[2] = 1; imask[3] = 0;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 1);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3] - 1) < tol);
    assert(fabs(weights[0] - (1 - 0.6 - 0.2)) < tol);
    assert(fabs(weights[1] - 0.6) < tol);
    assert(fabs(weights[2] - 0.2) < tol);
    assert(fabs(weights[3] - 0) < tol);

    // 2 missing data, cannot interpolate
    imask[0] = 0; imask[1] = 1; imask[2] = 0; imask[3] = 1;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -2);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3]) < tol);
  }

  /* 3d tests */
  {
    int status;
    const int ndims = 3;
    const int nnodes = 8;
    const int dims[] = {2, 2, 2};
    double weights[nnodes];

    // all data valid
    double dindices[] = {0.2, 0.4, 0.6};
    int imask[] = {1, 1, 1, 1, 1, 1, 1, 1};
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 0);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3] +
                weights[4] + weights[5] + weights[6] + weights[7]
                - 1) < tol);
    assert(fabs(weights[0] - 
                (1-dindices[0])*(1-dindices[1])*(1-dindices[2]) < tol));
    assert(fabs(weights[1] - 
                (1-dindices[0])*(1-dindices[1])*(0+dindices[2]) < tol));
    assert(fabs(weights[2] - 
                (1-dindices[0])*(0+dindices[1])*(1-dindices[2]) < tol));
    assert(fabs(weights[3] - 
                (1-dindices[0])*(0+dindices[1])*(0+dindices[2]) < tol));
    assert(fabs(weights[4] - 
                (0+dindices[0])*(1-dindices[1])*(1-dindices[2]) < tol));
    assert(fabs(weights[5] - 
                (0+dindices[0])*(1-dindices[1])*(0+dindices[2]) < tol));
    assert(fabs(weights[6] - 
                (0+dindices[0])*(0+dindices[1])*(1-dindices[2]) < tol));
    assert(fabs(weights[7] - 
                (0+dindices[0])*(0+dindices[1])*(0+dindices[2]) < tol));

    // missing values, ok to interpolate
    imask[0] = 1; imask[1] = 1; imask[2] = 0; imask[3] = 1;
    imask[4] = 0; imask[5] = 1; imask[6] = 0; imask[7] = 0;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == 4);
    assert(fabs(weights[0] + weights[1] + weights[2] + weights[3] +
                weights[4] + weights[5] + weights[6] + weights[7] 
                - 1) < tol);
    assert(fabs(weights[5] - 0.2) < tol);
    assert(fabs(weights[3] - 0.4) < tol);
    assert(fabs(weights[0] - 0.4) < tol);
    
    // too many missing values
    imask[0] = 0; imask[1] = 1; imask[2] = 0; imask[3] = 1;
    imask[4] = 0; imask[5] = 1; imask[6] = 0; imask[7] = 0;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -5);

    imask[0] = 0; imask[1] = 0; imask[2] = 0; imask[3] = 1;
    imask[4] = 0; imask[5] = 1; imask[6] = 0; imask[7] = 0;
    status = nccf_get_linear_weights_double(ndims, dims, 
                                            dindices,
                                            imask, weights);
    assert(status == -6);
  }
#endif // HAVE_LAPACK_LIB

  return 0;
}
