/**
 * API for data on curvilinear coordinates
 *
 * $Id: nccf_data.h 923 2012-03-22 18:44:09Z dkindig $
 */

#ifndef _NCCF_STRUCTURED_DATA_H
#define _NCCF_STRUCTURED_DATA_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

extern struct CFLISTITEM *CFLIST_STRUCTURED_DATA;

struct nccf_data_type {

  /* grid object to which the data are attach to */
  int gridid;

  /* type of data (compatible with netcdf types) */
  nc_type dataType;

  /* name of the data */
  char *name;

  /* number of dimensions */
  int ndims;

  /* netcdf-like variable object */
  struct nccf_var_obj *dataVar;

  /* whether or not data are copy-saved by the object */
  int save;

  /* number of written records */
  int numRecords;

  /* pointer to the data */
  void *data;
};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_data(int gridid, const char *name, 
		  const char *standard_name,
		  const char *units, const char *time_dimname,
		  int *dataid);

int nccf_def_data_from_file(const char *filename, 
			    int gridid, const char *varname, 
			    int read_data, int *dataid);

int nccf_free_data(int dataid);

int nccf_inq_data_gridid(int dataid, int *gridid);

int nccf_set_data_double(int dataid, const double *data, 
			 int save, double fill_value);

int nccf_set_data_float(int dataid, const float *data, 
			int save, float fill_value);

int nccf_set_data_int(int dataid, const int *data, 
		      int save, int fill_value);

int nccf_set_data_short(int dataid, const short *data, 
			int save, short fill_value);

int nccf_put_data(int ncid, int dataid);

int nccf_add_data_att(int dataid, 
		      const char *name, const void *value);

int nccf_inq_data_type(int dataid, nc_type *dataType);

int nccf_get_data_pointer(int dataid, nc_type *xtypep, 
			  void **dataPtr, const void **fill_valuep);

int nccf_inq_data_ndims(int dataid, int *nDataDims);
  
int nccf_inq_data_dims(int dataid, int *dataDims);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // _NCCF_STRUCTURED_DATA_H

