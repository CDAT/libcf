/**
 * Utility functions
 * $Id: nccf_utility_functions.h 915 2012-01-09 16:58:08Z pletzer $
 */

#ifndef _NCCF_GLOBAL_FUNCTIONS
#define _NCCF_GLOBAL_FUNCTIONS

#ifdef __cplusplus
extern "C" {
#endif

/*! \defgroup gs_tools_grp Tools
    \ingroup gridspec_grp

Convenience functions for generating coordinate data for
cubed-sphere and tripolar grids, for obtaining index ranges 
for given tile boundaries, mapping sets of 
indices to a single ``big'' index, generating a unique 
coordinate system Ids, and for interpolating.

@{*/

/** 
 * Convert an index (0... 3^ndims - 1) to a corner index vector in ndims
 * whose elements are -1, 0, or 1. -1 refers to the lower side and +1 to
 * the upper side while 0 is the interior. 
 *
 * \param index (input) from loop over 3^ndims 
 * \param ndims number of dimensions
 * \param vector (output) array of ndims elements containing 0, 1, or -1, all elements are zero except one
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_index_to_corner_vector(int index, int ndims, int vector[]);

/** 
 * Convert an index (0... 2^ndims - 1) to an index vector in ndims
 *
 * \param index input from loop over 2^ndims 
 * \param ndims number of dimensions
 * \param vector (output) array of ndims elements containing 0, 1, or -1, 
 *                all elements are zero except one
 * \return NC_NOERR on success
 *
 * \author David Kindig, Tech-X Corp.
 */
int nccf_index2vector( int index, int ndims, int vector[] );

/**
 * Compute the flat array index of a multi-dim index set.
 *
 * \param ndims number of elements of index set
 * \param dims dimensions of each axis
 * \param inx index array
 * \return flat array index
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
int nccf_get_flat_index(int ndims, const int dims[], const int inx[]);

/**
 * Compute the multi-index set from a flat array index.
 *
 * \param ndims number of elements of index set
 * \param dims dimensions of each axis
 * \param index flat index
 * \param inx (output) index array
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
void nccf_get_multi_index(int ndims, const int dims[], int index, int inx[]);

/**
 * Compute the start and end indices of a grid boundary.
 *
 * \param ndims number of elements of index set
 * \param dims dimensions of each axis
 * \param normalVect vector of -1, 0, and 1 identifying the boundary
 * \param exclusive set to 1 if endIndices are exclusive (1 past last), 0 otherwise
 * \param startIndices (output) index vector of start position
 * \param endIndices (output) index vector of end position
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
void nccf_get_start_end_bound_indices(int ndims, const int dims[], 
				      const int normalVect[],
				      int exclusive,
				      int startIndices[],
				      int endIndices[]);
/**
 * Compute the longitudes and latitudes of a cubed sphere grid.
 *
 * \param dims the 2 dimensions of each axis
 * \param faceVect 3-element vector perpendicular to the face
 * \param lon (output) longitudes in deg (memory managed by user)
 * \param lat (output) latitudes  in deg (memory managed by user)
 * \return NC_NOERR on success
 *
 * \author Matthew Wrobel, Alexander Pletzer, and David Kindig, Tech-X Corp.
 */
int nccf_get_cubedsphere_grid(const int dims[], const int faceVect[], 
			      double lon[], double lat[]);

/**
 * Compute the longitudes and latitudes of the cap portion of a tripolar grid.
 *
 * \param dims the # points on each axis
 * \param latPerim the latitude at which the north cap region begins
 * \param lonSing the longitude of one of the two poles, the other +-180 deg away
 * \param lon (output) longitudes in deg (memory managed by user)
 * \param lat (output) latitudes  in deg (memory managed by user)
 * \return NC_NOERR on success
 *
 * \author Matthew Wrobel and Alexander Pletzer, Tech-X Corp.
 */
int nccf_get_bipolar_cap(const int dims[], 
			 double latPerim, double lonSing, 
			 double lon[], double lat[]);

/**
 * Compute the longitudes and latitudes of entire tripolar grid.
 *
 * \param dims the # points on each axis
 * \param capIndx the index at which the north cap region begins
 * \param lon (output) longitudes in deg (memory managed by user)
 * \param lat (output) latitudes  in deg (memory managed by user)
 * \return NC_NOERR on success
 *
 * \author Matthew Wrobel, Tech-X Corp.
 */
 int nccf_get_tripolar_grid(const int dims[], int capIndx,
    double lon[], double lat[]);

/**
 * Compute the longitudes and latitudes of the east and west halves of entire tripolar grid.
 *
 * \param wchHalf which half, 0 is [-180,0] and 1 is [0,180]
 * \param dims the # points on each axis
 * \param capIndx the index at which the north cap region begins
 * \param lon (output) longitudes in deg (memory managed by user)
 * \param lat (output) latitudes  in deg (memory managed by user)
 * \return NC_NOERR on success
 *
 * \author Matthew Wrobel, Tech-X Corp.
 */
int nccf_get_tripolar_halfgrid(const int dims[], int wchHalf, int capIndx,
    double lon[], double lat[]);

/**
 * Write multi-dimensional array of strings by blocks of size dims[ndims - 1].
 *
 * \param ncid netcdf file id
 * \param varId variable id
 * \param ndims number of dimensions
 * \param dims sizes for each dimension
 * \param string multi-dim string array to write
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
int nccf_put_var_slice_text(int ncid, int varId, int ndims, const int *dims, 
			const char *string);

/**
 * Create an pseudo-random alpha-numeric id.
 *
 * \param seed_adj a user given value to adjust the seed somewhat. This works only if include/uuid/uuid.h does not exist.
 * \param id (output) a pseudo-random 32 digit string in uuid format (need to allocate for 36 characters at minimum)
 * \return NC_NOERR on success
 *
 * \author David Kindig, Tech-X Corp.
 */
int nccf_generate_id(int seed_adj, char *id);

/**
 * Add the id from nccf_generate_id to a list of files.
 *
 * \param id id generated by nccf_generate_id
 * \param nfiles number of files
 * \param tile_names char array of tile_names 
 * \param file_types char array of file_types
 * \param list_o_files Char array of filenames
 * \return NC_NOERR on success
 *
 * \author David Kindig, Tech-X Corp.
 */
int nccf_add_id_to_files(const char *id, int nfiles, const char **tile_names, 
                         const char **file_types, const char **list_o_files);

/**
 * Compute the next adjacent index to a given index.
 * \param ndims number of elements of index set
 * \param dims dimensions of each axis
 * \param index flat index 
 * 
 * \return 1 if the adjacent index is inside the domain, 0 otherwise
 * \note the return flag can be used to terminate a while loop.
 */

int nccf_find_next_adjacent(int ndims, const int dims[], int index[]);

/**
 * Find next index position iterate. No checks are performed to 
 * ensure that the returned index position is within the domain. 
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param coordData list of grid vertex arrays, dimensioned [ndims][dims[0]*dims[1]...]
 * \param coord_periodicity periodicity length in coordinate space, for each axis
 * \param targetPos target position, dimension [ndims]
 * \param dIndices_in initial guess for position in index space on input. Array of size [ndims].
 * \param dIndices_out output index position. Array of size [ndims].
 * \param position_out output position. Array of size [ndims].
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
  int nccf_find_next_indices_double(int ndims, const int dims[], 
				    const double **coordData,
                                    const double coord_periodicity[],
				    const double targetPos[], 
				    const double dIndices_in[], 
				    double dIndices_out[],
				    double position_out[]);
  int nccf_find_next_indices_float(int ndims, const int dims[], 
				   const float **coordData, 
                                   const float coord_periodicity[],
				   const float targetPos[], 
				   const float dIndices_in[], 
				   float dIndices_out[],
				   float position_out[]);

/**
 * Compute distance between position given by indices and target position
 * in coordinate space.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param coordData list of grid vertex arrays, dimensioned [ndims][dims[0]*dims[1]...]
 * \param coord_periodicity periodicity length in coordinate space, for each axis
 * \param targetPos target position, dimension [ndims]
 * \param dindices previous guess for position in index space. Array of size [ndims].
 * \param distance distance from targetPos (output).
 * \return NC_NOERR on success
 */

int nccf_get_distance_in_coord_double(int ndims, const int dims[], 
                                      const double **coordData, 
                                      const double coord_periodicity[],
                                      const double targetPos[],
                                      const double dindices[],
                                      double *distance);
int nccf_get_distance_in_coord_float(int ndims, const int dims[], 
                                     const float **coordData, 
                                     const float coord_periodicity[],
                                     const float targetPos[],
                                     const float dindices[],
                                     float *distance);

/**
 * Perform a line search minimization in the direction of 
 * dindices_new - dindices.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param coordData list of grid vertex arrays, dimensioned [ndims][dims[0]*dims[1]...]
 * \param coord_periodicity periodicity length in coordinate space, for each axis
 * \param targetPos target position, dimension [ndims]
 * \param indexTol search stops when distance in index space < indexTol
 * \param dindices previous guess for position in index space. Array of size [ndims].
 * \param dindices_new (input and output) last guess for position in index space, will be corrected on output. Array of size [ndims].
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_linesearch_indices_double(int ndims, const int dims[], 
				   const double **coordData, 
                                   const double coord_periodicity[],
				   const double targetPos[],
				   double indexTol, 
				   const double dindices[], 
				   double dindices_new[]);
int nccf_linesearch_indices_float(int ndims, const int dims[], 
				  const float **coordData, 
                                  const float coord_periodicity[],
				  const float targetPos[],
				  float indexTol, 
				  const float dindices[], 
				  float dindices_new[]);

/**
 * Find index position.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param coordData list of grid vertex arrays, dimensioned 
 *                  [ndims][dims[0]*dims[1]...]
 * \param coord_periodicity periodicity length in coordinate space, for each axis
 * \param targetPos target position, dimension [ndims]
 * \param niter max number of iterations on input, used number of iterations
 *                 on output.
 * \param tolpos max tolerance in position space on input, achieved
 *               tolerance on output.
 * \param adjustFunc function supplied by the user, which is called to 
 *                   adjust the indices when index leaves the domain. 
 *                   Supply NULL if you want to exit when indices leave 
 *                   the domain. 
 * \param dindices initial guess for position in index space on input, 
 *                 index position on output. Array of size [ndims].
 * \param hit_bounds array of -1 (target is below original grid boundary, 
                     0 (target inside boundary), and 1 (target is above 
                     original grid boundary). 
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_find_indices_double(int ndims, const int dims[], 
                             const double **coordData, 
                             const double coord_periodicity[],
                             const double targetPos[], 
                             int *niter, double *tolpos, 
                             void (*adjustFunc)(int, const int[], double []), 
                             double dindices[], int hit_bounds[]);
int nccf_find_indices_float(int ndims, const int dims[], 
                            const float **coordData, 
                            const float coord_periodicity[],
                            const float targetPos[], 
                            int *niter, float *tolpos, 
                            void (*adjustFunc)(int, const int[], float []), 
                            float dindices[], int hit_bounds[]);
 

/**
 * Get the position in coordinate space given a set of indices.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param coordData list of grid vertex arrays, dimensioned 
 *                  [ndims][dims[0]*dims[1]...]
 * \param coord_periodicity periodicity length in coordinate space, for 
 *                          each axis
 * \param pos_ref modulo operations will be applied for periodic positions so 
 *                as to minimize the distance to pos_ref (passing NULL will 
 *                not apply modulo)
 * \param dindices index position returned by nccf_find_indices_double
 * \param position (output) position in coordinate space
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_get_position_double(int ndims, const int dims[],
			     const double **coordData,
                             const double coord_periodicity[],
                             const double pos_ref[],
			     const double dindices[],
			     double position[]);
int nccf_get_position_float(int ndims, const int dims[],
                            const float **coordData,
                            const float coord_periodicity[],
                            const float pos_ref[],
                            const float dindices[],
                            float position[]);

/**
 * Get linear interpolation weights.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param dindices index position returned by nccf_find_indices_double
 * \param imask optional mask array (1=valid data, 0=invalid data), can be NULL (meaning all data are valid)
 * \param weights (output) interpolation weights (array of dimension 2^ndims)
 *
 * \return NC_NOERR on success, a positive number indicates that invalid 
 * data were detected, a negative number indicates that there were too many
 * invalid data to allow interpolation.
 * 
 * \see nccf_find_indices_double
 * \note the sum of the weights typically amount to 1. In the case of invalid 
 * data, zero weights may be returned. 
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_get_linear_weights_double(int ndims, const int dims[], 
				   const double dindices[], 
				   const int imask[],
				   double weights[]);
int nccf_get_linear_weights_float(int ndims, const int dims[], 
				  const float dindices[], 
				  const int imask[],
				  float weights[]);

/**
 * Perform a linear interpolation for given index position.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param f_nodes grid values of the function to interpolate
 * \param dindices interpolation point in index space (returned by 
 *                 nccf_find_indices_double/nccf_find_indices_float)
 * \param weights interpolation weights (returned by 
 *                nccf_get_linear_interp_weights)
 * \param fill_value the value f_interp will be set to if too many nodal 
 *                   function values are invalid.
 * \param f_interp (output) interpolated value
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_linear_interp_double(int ndims, const int dims[],
			      const double *f_nodes,
			      const double dindices[],
			      const double weights[],
			      double fill_value,
			      double *f_interp);
int nccf_linear_interp_float(int ndims, const int dims[], 
                             const float *f_nodes,
                             const float dindices[],
                             const float weights[],
                             float fill_value,
                             float *f_interp);

/**
 * Perform a linear interpolation for given index position, if necessary 
 * correcting the indices to account for periodicity.
 *
 * \param ndims number of space dimensions
 * \param dims dimensions along each axis
 * \param f_nodes grid values of the function to interpolate
 * \param dindices interpolation point in index space (returned by 
 *                 nccf_find_indices_double/nccf_find_indices_float)
 * \param weights interpolation weights (returned by 
 *                nccf_get_linear_interp_weights)
 * \param fill_value the value f_interp will be set to if too many nodal 
 *                   function values are invalid.
 * \param f_interp (output) interpolated value
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_linear_interp_mod_double(int ndims, const int dims[],
                                  const double *f_nodes,
                                  double f_ref, 
                                  double f_periodicity, 
                                  const double dindices[],
                                  const double weights[],
                                  double fill_value,
                                  double *f_interp);
int nccf_linear_interp_mod_float(int ndims, const int dims[], 
                                 const float *f_nodes,
                                 float f_ref, 
                                 float f_periodicity, 
                                 const float dindices[],
                                 const float weights[],
                                 float fill_value,
                                 float *f_interp);

/**
 * Solve a dense float linear system.
 *
 * \param ndims number of space dimensions
 * \param mat matrix containing ndims*ndims elements
 * \param rhs right-hand side vector of dimension ndims
 * \param sol solution vector of size ndims
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */
int nccf_solve_float(int ndims, float mat[], const float rhs[], 
		      float sol[]);
int nccf_solve_double(int ndims, double mat[], const double rhs[], 
		      double sol[]);

/*!@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NCCF_GLOBAL_FUNCTIONS */


