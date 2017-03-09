/*
 * API for axes.
 *
 * $Id: nccf_axis.h 834 2011-09-14 21:21:07Z pletzer $
 */

#ifndef _NCCF_AXIS_H
#define _NCCF_AXIS_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

extern const int nccf_num_axis_types;
extern const char* nccf_axis_type_name[];

extern struct CFLISTITEM *CFLIST_AXIS;

struct nccf_axis_type {

  /* dimension length */
  int len;

  /* name of this coordinate, eg "lon" */
  char *axis_name;

  /* netcdf-like object */
  struct nccf_var_obj *axisVar;

  /* pointer to the data */
  void *data;

  /* type of data */
  nc_type xtype;
};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_axis(const char *name, int len, nc_type xtype, 
                  const void *data, const char *standard_name,
                  const char *units, int cdm_axis_type, 
                  const char *axis, int positive_up,  
                  const char *formula_terms, int *axisid);

int nccf_def_axis_from_file(const char *filename, const char *name,
                            int *axisid);

int nccf_free_axis(int axisid);

int nccf_put_axis(int axisid, int ncid);

int nccf_inq_axis_datatype(int axisid, nc_type *xtype);

int nccf_inq_axis_name(int axisid, char *name);

int nccf_inq_axis_len(int axisid, int *len);

int nccf_add_axis_att(int axisid, const char *name, const char *value);

int nccf_get_axis_datapointer(int axisid, void **datapp);
                    
/*!@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _NCCF_AXIS_H */

