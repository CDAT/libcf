/* $Id: nccf_get_distance_in_coord.h 905 2011-12-29 04:56:48Z pletzer $ 
 *  \author Alexander Pletzer, Tech-X Corp.
 */

int
nccf_get_distance_in_coord_TYPE(int ndims, const int dims[], 
				const _TYPE **coordData,
                                const _TYPE coord_periodicity[],
				const _TYPE targetPos[],
				const _TYPE dindices[],
                                _TYPE *distance) {
  int i;
  _TYPE xPos[ndims];
  nccf_get_position_TYPE(ndims, dims, coordData, 
                         coord_periodicity, targetPos, 
                         dindices, xPos);
  *distance = 0;
  _TYPE diff;
  for (i = 0; i < ndims; ++i) {
    diff = targetPos[i] - xPos[i];
    *distance += diff * diff;
  }
  *distance = sqrt(*distance);
  return 0;
}
