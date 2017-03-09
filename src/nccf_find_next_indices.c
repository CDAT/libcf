/**
 * $Id: nccf_find_next_indices.c 903 2011-12-28 21:46:08Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <nccf_utility_functions.h>
#include <nccf_constants.h>

#include <cf_config.h>

#define nccf_find_next_indices_TYPE nccf_find_next_indices_double
#define _TYPE double
#define nccf_get_position_TYPE nccf_get_position_double
#define nccf_solve_TYPE nccf_solve_double
#define _HUGE_TYPE CF_HUGE_DOUBLE
#include <nccf_find_next_indices.h>
#undef _HUGE_TYPE
#undef nccf_solve_TYPE
#undef nccf_get_position_TYPE
#undef _TYPE
#undef nccf_find_next_indices_TYPE

#define nccf_find_next_indices_TYPE nccf_find_next_indices_float
#define _TYPE float
#define nccf_get_position_TYPE nccf_get_position_float
#define nccf_solve_TYPE nccf_solve_float
#define _HUGE_TYPE CF_HUGE_FLOAT
#include <nccf_find_next_indices.h>
#undef _HUGE_TYPE
#undef nccf_solve_TYPE
#undef nccf_get_position_TYPE
#undef _TYPE
#undef nccf_find_next_indices_TYPE



