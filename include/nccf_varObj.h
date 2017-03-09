/**
 * $Id: nccf_varObj.h 895 2011-12-22 13:12:38Z pletzer $
 *
 * Methods to populate a netcdf-like variable
 */

#ifndef _NCCF_VAROBJ_H
#define _NCCF_VAROBJ_H

#include <netcdf.h>
#include "nccf_keyvalue.h"

/**
 * Struct mimicking a netcdf variable but stored in memory
 * A variable has data (double. float , ..) and attributes,
 * currently of either string, double, float, or int type.
 * The attributes are accessed by string name. The variable 
 * knows how to write itself to a netcdf file and how to load 
 * itself from a netcdf file.
 */
struct nccf_var_obj{

  /* name */
  char *name;

  /* list of dimension names */
  char **dimnames;

  /* number of dimensions */
  int ndims;

  /* list of dimensions */
  int *dims;

  /* pointer to the data */
  void *data;

  /* netcdf data type (eg NC_DOUBLE) */
  nc_type data_type;

  /* map string to string */
  struct CFLISTITEM *attr;

  /* whether or not this object owns the data */
  int save;

  /* number of written records for variables with a time-like axis */
  int numWrittenRecords;

  /* time dimension (-1 if static) */
  int time_dimension;
};

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Constructor
 * handle - object handle
 * name   - name of variable, if empty will denote a global attribute
 *          name will be strdup-ed, so the callere owns the pointer.
 * \return 0 on success
 */
int nccf_varCreate(struct nccf_var_obj **handle, const char* name);

/**
 * Constructor from netcdf file
 * handle       - object handle
 * name         - name of variable, if empty will denote a global attribute
 * ncid         - NetCDF file id 
 * readData     - set to 1 if data should be read
 * castToDouble - set to 1 if data should be cast into a double array
 * \return 0 on success
 */
int nccf_varCreateFromFile(struct nccf_var_obj **handle, const char* name, 
			   int ncid, int readData, int castToDouble);

/**
 * Destructor
 * handle - object handle
 */
int nccf_varDestroy(struct nccf_var_obj **handle);

/**
 * Retrieve the variable name
 * handle  - object handle
 * varname - variable name
 *           the pointer is owned by the struct handle.
 * \return 0 on success
 */
int nccf_varGetVarNamePtr(struct nccf_var_obj **handle, const char **varname);

/**
 * name   - name of variable, if empty will denote a global attribute
 *          name will be strdup-ed, so the caller owns the pointer.
 * \return 0 on success
 */
int nccf_varSetVarName(struct nccf_var_obj **handle, const char *varname);

/**
 * Set the attribute iterator to the beginning
 * handle - object handle
 * \return 0 on success
 */
int nccf_varAttribIterBegin(struct nccf_var_obj **handle);

/**
 * Set the attribute iterator to the next value
 * handle - object handle
 * \return 1 if the next attribute exists, 0 otherwise
 */
int nccf_varAttribIterNext(struct nccf_var_obj **handle);

/**
 * Inquire the current attribute name
 * handle - object handle
 * name   - name of the attribute
 * \return 0 on success
 */
int nccf_varInqAttribNamePtr(struct nccf_var_obj **handle, const char **name);

/**
 * Attach an attribute to variable
 * handle - object handle
 * name   - attribute name
 * value  - value of the attribute
 * \return 0 on success
 */
int nccf_varSetAttribText(struct nccf_var_obj **handle, const char *name, 
			  const char *value);
int nccf_varSetAttribDouble(struct nccf_var_obj **handle, const char *name, 
			    double value);
int nccf_varSetAttribFloat(struct nccf_var_obj **handle, const char *name, 
			   float value);
int nccf_varSetAttribInt(struct nccf_var_obj **handle, const char *name, 
			 int value);
int nccf_varSetAttribShort(struct nccf_var_obj **handle, const char *name, 
			   short value);
/**
 * Attach an array attribute to variable
 * handle - object handle
 * name   - attribute name
 * n      - number of elements
 * values - attribute values
 * \return 0 on success
 */
int nccf_varSetAttribDoubleArray(struct nccf_var_obj **handle, 
                                 const char *name,  int n, 
                                 const double values[]);
int nccf_varSetAttribFloatArray(struct nccf_var_obj **handle, 
                                const char *name, int n, 
                                const float values[]);
int nccf_varSetAttribIntArray(struct nccf_var_obj **handle, 
                              const char *name,  int n, 
                              const int values[]);
int nccf_varSetAttribShortArray(struct nccf_var_obj **handle, 
                                const char *name,  int n, 
                                const short values[]);

/**
 * Inquire about an attribute value
 * handle - object handle
 * name   - attribute name
 * xtypep - a valid netcdf type
 * lenp   - number of elements
 * \return 0 on success
 */
int nccf_varInqAttrib(struct nccf_var_obj **handle, const char *name,
		      nc_type *xtypep, int *lenp);

/**
 * Get the attribute pointer value by name
 * handle - object handle
 * name   - name of the attribute
 * value  - The value of the attribute. The caller does not own the pointer
 *          returned in *value. This pointer will be freed when *v is freed.
 * \return 0 on success
 */
int nccf_varGetAttribPtr(struct nccf_var_obj **handle, 
			  const char *name, const void **value);

/**
 * Get number of data values per time slice
 * handle  - object handle
 * ntot    - (output) number of data values
 * \return 0 on success
 */
int nccf_varGetNumValsPerTime(struct nccf_var_obj **handle, int *ntot);

/**
 * Get number of dimensions for a variable
 * handle  - object handle
 * numDims - (output) number of dimension attached to the variable
 * \return 0 on success
 */
int nccf_varGetNumDims(struct nccf_var_obj **handle, int *numDims);

/**
 * Set the dimensions for a variable
 * handle   - object handle
 * numDims  - number of dimension attached to the variable
 * dims     - array of dimensions
 * dimnames - the names of each dimension (will be copied)
 * \return 0 on success
 */
int nccf_varSetDims(struct nccf_var_obj **handle, int numDims, 
 	       const int dims[], const char **dimname);

/**
 * Get the dimension names for a variable as a pointer
 * handle   - object handle
 * index    - index of the dimension
 * dimname  - the name of the dimension (object owns the pointer)
 * \return 0 on success
 */
int nccf_varGetDimNamePtr(struct nccf_var_obj **handle, int index, 
                         const char **dimname);

/**
 * Get a pointer to the dimensions of a variable object.
 * handle - object handle
 * dim    - the dimensions
 * \return 0 on success
 */
int nccf_varGetDimsPtr(struct nccf_var_obj **handle, int **dims);

/** 
 * Set single time slice data pointer (caller owns the pointer)
 * handle   - object handle
 * dataType - netcdf var type (e.g. NC_FLOAT)
 * val      - values
 * \return 0 on success
 * \note should be large enough to accommodate tyhe values 
 * for a single time slice
 */
int nccf_varSetDataPtr(struct nccf_var_obj **handle, 
                       nc_type dataType, void *val);

/**
 * Read a time slice
 * handle       - object handle
 * ncid         - netcdf file handle
 * time_index   - time index (0 if static)
 * castToDouble - set to 1 if data should be read as doubles
 */
  int nccf_varReadData(struct nccf_var_obj **handle, int ncid, 
                       int time_index, int castToDouble);

/** 
 * Set single time slice data (will make a copy)
 * handle   - object handle
 * val      - values
 * \return 0 on success
 */
int nccf_varSetDataDouble(struct nccf_var_obj **handle, const double val[]);

/** 
 * Set single time slice data (will make a copy)
 * handle   - object handle
 * val      - values
 * \return 0 on success
 */
int nccf_varSetDataFloat(struct nccf_var_obj **handle, const float val[]);

/** 
 * Set single time slice data (will make a copy)
 * handle   - object handle
 * val      - values
 * \return 0 on success
 */
int nccf_varSetDataInt(struct nccf_var_obj **handle, const int val[]);

/** 
 * Set single time slice data (will make a copy)
 * handle   - object handle
 * val      - values
 * \return 0 on success
 */
int nccf_varSetDataShort(struct nccf_var_obj **handle, const short val[]);

/** 
 * Set single time slice data (will make a copy)
 * handle   - object handle
 * val      - values
 * \return 0 on success
 */
int nccf_varSetDataChar(struct nccf_var_obj **handle, const char val[]);

/**
 * Get the data type
 * handle   - object handle
 * dataType - (output) e.g. NC_FLOAT
 * \return 0 on success
 */
int nccf_varGetDataType(struct nccf_var_obj **handle, nc_type *dataType);

/** 
 * Get data pointer to a single time slice (object owns the pointer)
 * handle   - object handle
 * val      - (output) value pointer
 * \return 0 on success
 */
int nccf_varGetDataPtr(struct nccf_var_obj **handle, void **val);

/** 
 * Write a list of variables to a NetCDF file
 * 
 * ncid     - NetCDF file handle
 * numVars  - number of var object handles
 * ...      - var object handles
 * \return 0 on success
 */
int nccf_writeListOfVars(int ncid, int numVars, ...);

/** 
 * Write data attached to give variables. This will write a new record for 
 * variables that have the NC_UNLIMITED dimension. Use nccf_writeListOfVars 
 * to write the header.
 *
 * ncid     - NetCDF file handle
 * numVars  - number of var object handles
 * ...      - var object handles
 * \return 0 on success
 */
int nccf_writeListOfVarData(int ncid, int numVars, ...);

#ifdef __cplusplus
}
#endif

#endif // _NCCF_VAROBJ_H
