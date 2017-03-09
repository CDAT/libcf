/**
 * API to define a structured grid for the gridspec convention to libcf.
 *
 * $Id: nccf_coord.h 895 2011-12-22 13:12:38Z pletzer $
 */

#ifndef _NCCF_COORDINATE_H
#define _NCCF_COORDINATE_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

extern struct CFLISTITEM *CFLIST_COORDINATE;

#define CF_COORD_LON_VARNAME "lon"
#define CF_COORD_LON_STNAME "longitude"
#define CF_COORD_LON_UNITS "degrees_east"
#define CF_COORD_LAT_VARNAME "lat"
#define CF_COORD_LAT_STNAME "latitude"
#define CF_COORD_LAT_UNITS "degrees_north"

struct nccf_coord_type {

  /* name of this coordinate, eg "lon" */
  char *coord_name;

  /* netcdf-like object */
  struct nccf_var_obj *coordVar;

  /* pointer to the data, only double supported */
  double *data;

  /* whether or not data are copy-saved (!= 0) by the object. If copy
     saved then the object owns the data, the data will be deallocated
     when destroying the object. Otherwise, the caller owns the data
     and is responsible for their destruction (after freeing this object.)
  */
  int save;
};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_coord(int ndims, const int *dims, const char **dimnames,
		   const double *data, int save, const char *name, 
		   const char *standard_name, const char *units,
		   int *coordid);

int nccf_def_lon_coord(int ndims, const int *dims, const char **dimnames,
		       const double *data, int save, int *coordid );

int nccf_def_lat_coord(int ndims, const int *dims, const char **dimnames,
		       const double *data, int save, int *coordid);

int nccf_def_coord_from_file(const char *filename, 
			     const char *coord_name, int *coordid);

int nccf_def_coord_from_axes(int ndims, const int *axisids, 
                             int index_pos, const char *name, 
                             const char *standard_name, 
                             const char *units, int *coordid);

int nccf_free_coord(int coordid);

int nccf_put_coord(int ncid, int coordid);

int nccf_add_coord_att(int coordid, const char *name, const char *value);

int nccf_inq_coord_bound(int coordid, const int norm_vect[],
         int *start_indices, int *end_indices);

int nccf_inq_coord_bound_slice(int coordid, const int norm_vect[],
			       int flip, const char *format, char *slice);

int nccf_get_coord_data_pointer(int coordid, double **data);

int nccf_inq_coord_ndims(int coordid, int *ndims);

int nccf_inq_coord_dims(int coordid, int *dims);

int nccf_inq_coord_dimnames(int coordid, char **dimnames);

int nccf_inq_coord_name(int coordid, char *coord_name);

int nccf_inq_coord_units(int coordid, char *unit);

int nccf_is_coord_lon(int coordid, int *yesno);

int nccf_is_coord_lat(int coordid, int *yesno);

int nccf_is_coord_vert(int coordid, int *yesno);

/*!@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _NCCF_COORDINATE_H */

