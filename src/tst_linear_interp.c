/**
 * $Id: tst_linear_interp.c 892 2011-12-21 19:48:51Z pletzer $
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <netcdf.h>
#include <nccf_constants.h>
#include <nccf_utility_functions.h>
#include <nccf_handle_error.h>

#ifdef HAVE_CONFIG_H
#include <cf_config.h>
#endif

void test1d() {
  int status;
  const int ndims = 1;
  const int dims[] = {11};
  const double xmin[] = {0.0};
  const double xmax[] = {1.0};
  // set very large periocity lengths for non-periodic case
  const double coord_periodicity[] = {CF_HUGE_DOUBLE};
  // target
  double x0[ndims];
  int hit_bounds[ndims];
  x0[0] = xmin[0] + 0.45 * (xmax[0] - xmin[0]);
  const double tolpos = 1.e-3;
  const int nitermax = 20;
  double weights[2];

  double f;
  double dindices[ndims];
  double **xs;
  double *fs;
  int ntot = 1;
  int i, j;
  int inx[ndims];
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  xs = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    xs[i] = (double *) malloc(ntot * sizeof(double));
  }
  fs = (double *) malloc(ntot * sizeof(double));
  for (i = 0; i < ndims; ++i) {
    for (j = 0; j < ntot; ++j) {
      nccf_get_multi_index(ndims, dims, j, inx);
      xs[i][j] = xmin[i] + inx[i] * (xmax[i] - xmin[i])/(dims[i] - 1.0);      
    }
  }
  for (j = 0; j < ntot; ++j) {
    fs[j] = xs[0][j];
  }
  dindices[0] = 0;
  int niter = nitermax;
  double tol = tolpos;
  status = nccf_find_indices_double(ndims, dims, (const double **) xs, 
                                    coord_periodicity,
                                    x0, 
				    &niter, &tol, 
                                    NULL, dindices, hit_bounds);
  if (status) ERR;

  status = nccf_get_linear_weights_double(ndims, dims,
					  dindices, NULL, weights);
  if (status) ERR;

  status = nccf_linear_interp_double(ndims, dims,
                                     (const double *) fs, 
				     dindices, weights, NC_FILL_DOUBLE,
				     &f);
  if (status) ERR;

  free(fs);
  for (i = 0; i < ndims; ++i) {
    free(xs[i]);
  }
  free(xs);

  // check
  assert( fabs(f - x0[0]) < 1.e-6 );
}

//////////////////////////////////////////////////////////////////////

void test2d() {
  int status;
  const int ndims = 2;
  const int dims[] = {11, 21};
  const double xmin[] = {0.0, 0.0};
  const double xmax[] = {1.0, 2.0};
  // set very large periocity lengths for non-periodic case
  const double coord_periodicity[] = {CF_HUGE_DOUBLE, CF_HUGE_DOUBLE};
  // target
  double x0[ndims];
  int hit_bounds[ndims];
  x0[0] = xmin[0] + 0.43*(xmax[0] - xmin[0]);
  x0[1] = xmin[1] + 0.34*(xmax[1] - xmin[1]);
  const double tolpos = 1.e-3;
  const int nitermax = 20;
  double weights[2 * 2];

  double f;
  double dindices[ndims];
  double **xs;
  double *fs;
  int ntot = 1;
  int i, j;
  int inx[ndims];
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  xs = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    xs[i] = (double *) malloc(ntot * sizeof(double));
  }
  fs = (double *) malloc(ntot * sizeof(double));
  for (i = 0; i < ndims; ++i) {
    for (j = 0; j < ntot; ++j) {
      nccf_get_multi_index(ndims, dims, j, inx);
      xs[i][j] = xmin[i] + inx[i] * (xmax[i] - xmin[i])/(dims[i] - 1.0);      
    }
  }
  for (j = 0; j < ntot; ++j) {
    fs[j] = xs[0][j]*xs[1][j];
  }

  dindices[0] = 0;
  dindices[1] = 0;
  int niter = nitermax;
  double tol = tolpos;
  status = nccf_find_indices_double(ndims, dims, (const double **) xs,
                                    coord_periodicity,
                                    x0, 
				    &niter, &tol, 
                                    NULL, dindices, hit_bounds);
  if (status) ERR;

  status = nccf_get_linear_weights_double(ndims, dims,
					  dindices, NULL, weights);
  if (status) ERR;

  status = nccf_linear_interp_double(ndims, dims,
                                     (const double *) fs, 
				     dindices, weights, NC_FILL_DOUBLE,
				     &f);
  if (status) ERR;

  free(fs);
  for (i = 0; i < ndims; ++i) {
    free(xs[i]);
  }
  free(xs);

  // check
  assert( fabs(f - x0[0]*x0[1]) < 1.e-6 );
}

//////////////////////////////////////////////////////////////////////

void adjustIndicesFuncForPolarGrid(int ndims, const int dims[], double dindices[]) {

  if (dindices[0] < 0.0) {
    // when rho goes negative...

    // rho adjustment
    dindices[0] = fabs( dindices[0] );
    // theta adjustment
    dindices[1] += (dims[1] - 1.0) / 2.0; 
    dindices[1] = fmod(dindices[1], dims[1] - 1.0);
    if (dindices[1]< 0.0) dindices[1] += dims[1] - 1.0;
  }
  if (dindices[0] >= dims[0] - 1) {
    // when rho is beyond the last point
    dindices[0] = dims[0] - 1; // want to be sure that we are inside the domain
  }
}

void test2dPolar() {
  int status;
  const int ndims = 2;
  const int dims[] = {41, 129};
  // rho, theta, must avoid the singular point @ rho = 0
  const double xmin[] = {0.001, 0.0};
  const double xmax[] = {1.0, 2.0 * M_PI};
  // set very large periocity lengths for non-periodic case
  const double coord_periodicity[] = {CF_HUGE_DOUBLE, M_2_PI};
  // target
  double xTarget[ndims];
  int hit_bounds[ndims];
  double rhoTarget, theTarget;

  // tolerance
  const double tolpos = 1.e-3;
  // max number of iterations
  const int nitermax = 20;
  double weights[2 * 2]; // 2^ndims

  double f;
  double dindices[ndims];
  double **xs;
  double *fs;
  int ntot = 1;
  int i, j;
  int inx[ndims];
  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  xs = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    xs[i] = (double *) malloc(ntot * sizeof(double));
  }
  fs = (double *) malloc(ntot * sizeof(double));

  double dr = (xmax[0] - xmin[0]) / (double)(dims[0] - 1.0);
  double dt = (xmax[1] - xmin[1]) / (double)(dims[1] - 1.0);

  for (j = 0; j < ntot; ++j) {
    nccf_get_multi_index(ndims, dims, j, inx);
    double rho = xmin[0] + inx[0] * dr;
    double the = xmin[1] + inx[1] * dt;
    xs[0][j] = rho * cos(the); // x      
    xs[1][j] = rho * sin(the); // y
    fs[j] = 1.0 + rho*cos(2.0*the);
  }

  FILE *file = fopen("tst_linear_interp_polar2d.py", "w");
  fprintf(file, "from matplotlib import pylab\n");
  fprintf(file, "t = [\n");
  
  int iCase;
  const int nCases = 1000;
  for (iCase = 0; iCase < nCases; ++iCase) {
  
    // set the target point
    rhoTarget = xmin[0] + (dims[0] - 1) * dr * (double)random()/ (double)RAND_MAX;
    theTarget = xmin[1] + (dims[1] - 1) * dt * (double)random()/ (double)RAND_MAX;
    xTarget[0] = rhoTarget * cos(theTarget);
    xTarget[1] = rhoTarget * sin(theTarget);

    // initial guess
    dindices[0] = dims[0]/3.1235;
    dindices[1] = dims[1]/2.5678;

    // update dindices
    int niter = nitermax;
    double tol = tolpos;
    status = nccf_find_indices_double(ndims, dims, (const double **) xs,
                                      coord_periodicity,
                                      xTarget, 
				      &niter, &tol, 
				      adjustIndicesFuncForPolarGrid, 
				      dindices, hit_bounds);
    fprintf(file, "[%lf, %lf,  %d],\n", xTarget[0], xTarget[1], status);
    if (status) {
      printf("warning: nccf_find_indices_double did not converge for case %d (status = %d)\n",
	     iCase, status);
      printf("achieved tol = %lf\n", tol);
    }

    if (!status) {
      status = nccf_get_linear_weights_double(ndims, dims,
					      dindices, NULL, weights);
      if (status) ERR;

      status = nccf_linear_interp_double(ndims, dims,
                                         (const double *) fs, 
					 dindices, weights, NC_FILL_DOUBLE,
					 &f);
      // check
      double fx = 1.0 + rhoTarget*cos(2.0*theTarget);
#ifdef LOGGING
      printf("case: %d rho = %f the = %f f = %g fx = %g\n", iCase, 
	     rhoTarget, theTarget, f, fx);
#endif
      assert( fabs(f - fx) < 1.e-2 );
    }
  }
  fprintf(file, "]\n");
  fprintf(file, "for i in range(len(t)):\n");
  fprintf(file, "    if t[i][-1] == 0: pylab.plot([t[i][0]], [t[i][1]], 'go')\n"); // success
  fprintf(file, "    if t[i][-1] != 0: pylab.plot([t[i][0]], [t[i][1]], 'rx')\n"); // failure
  fprintf(file, "pylab.show()\n");
  fclose(file);

  free(fs);
  for (i = 0; i < ndims; ++i) {
    free(xs[i]);
  }
  free(xs);

}

//////////////////////////////////////////////////////////////////////

void test2dTripolar() {
  int status;
  const int ndims = 2;
  const int dims[] = {41, 41};
  // set very large periocity lengths for non-periodic case
  const double coord_periodicity[] = {CF_HUGE_DOUBLE, CF_HUGE_DOUBLE};
  // target
  double xTarget[ndims];
  const double latPerim = 60.0;

  // tolerance
  const double tolpos = 1.e-3;
  // max number of iterations
  const int nitermax = 20;
  double weights[2 * 2]; // 2^ndims

  double f;
  double dindices[ndims];
  int hit_bounds[ndims];
  double **xs;
  double *fs;
  int ntot = 1;
  int i, j;
  int failureCounter = 0;

  for (i = 0; i < ndims; ++i) {
    ntot *= dims[i];
  }

  xs = (double **) malloc(ndims * sizeof(double *));
  for (i = 0; i < ndims; ++i) {
    xs[i] = (double *) malloc(ntot * sizeof(double));
  }
  fs = (double *) malloc(ntot * sizeof(double));

  if ((status = nccf_get_bipolar_cap(dims, latPerim, -90.0, xs[1], xs[0]))) ERR;

  for (j = 0; j < ntot; ++j) {
    fs[j] = sin(2.0 * xs[0][j]/180.) * cos(xs[1][j]/180.0);
  }

  FILE *file = fopen("tst_linear_interp_tripolar.py", "w");
  fprintf(file, "from matplotlib import pylab\n");
  fprintf(file, "t = [\n");
  
  int iCase;
  const int nCases = 1000;
  for (iCase = 0; iCase < nCases; ++iCase) {
    // set the target point
    xTarget[0] = latPerim + (90.0 - latPerim) * (double)random()/ (double)RAND_MAX;
    xTarget[1] = -180 + 360.0 * (double)random()/ (double)RAND_MAX;

    // initial guess
    dindices[0] = dims[0]/3.1235;
    dindices[1] = dims[1]/2.5678;

    // update dindices
    int niter = nitermax;
    double tol = tolpos;
    status = nccf_find_indices_double(ndims, dims, (const double **) xs,
                                      coord_periodicity,
                                      xTarget, 
				      &niter, &tol, 
				      NULL, 
				      dindices, hit_bounds);
    fprintf(file, "[%lf, %lf,  %d],\n", xTarget[0], xTarget[1], status);
    if (status) {
      printf("warning: nccf_find_indices_double did not converge for case %d (status = %d)\n",
	     iCase, status);
      printf("achieved tol = %lf\n", tol);
      failureCounter++;
    }

    if (!status) {
      status = nccf_get_linear_weights_double(ndims, dims,
					      dindices, NULL, weights);
      if (status) ERR;

      status = nccf_linear_interp_double(ndims, dims,
                                         (const double *) fs, 
					 dindices, weights, NC_FILL_DOUBLE,
					 &f);
      // check
      double fx = sin(2.0 * xTarget[0]/180.) * cos(xTarget[1]/180.0);;
#ifdef LOGGING
      printf("case: %d lat = %f lon = %f f = %g fx = %g\n", iCase, 
	     xTarget[0], xTarget[1], f, fx);
#endif
      //assert( fabs(f - fx) < 1.e-2 );
    }
  }
  fprintf(file, "]\n");
  fprintf(file, "for i in range(len(t)):\n");
  fprintf(file, "    if t[i][-1] == 0: pylab.plot([t[i][1]], [t[i][0]], 'go')\n"); // success
  fprintf(file, "    if t[i][-1] != 0: pylab.plot([t[i][1]], [t[i][0]], 'rx')\n"); // failure
  fprintf(file, "pylab.show()\n");
  fclose(file);
  printf("%d failures in test2dTripolar\n", failureCounter);

  free(fs);
  for (i = 0; i < ndims; ++i) {
    free(xs[i]);
  }
  free(xs);

}

//////////////////////////////////////////////////////////////////////

int main() {

#ifdef HAVE_LAPACK_LIB

  /* 1d interpolation */
  //test1d();

  /* 2d interpolation */
  //test2d();
  //test2dPolar();
  test2dTripolar();

#endif

  return 0;
}
