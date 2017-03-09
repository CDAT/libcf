/*
 * $Id: nccf_def_axis.c 851 2011-11-08 14:37:20Z pletzer $
 */

#include "nccf_axis.h"
#include <stdlib.h>
#include <string.h>

#include <nccf_handle_error.h>
#include <nccf_varObj.h>

struct CFLISTITEM *CFLIST_AXIS;

const int nccf_num_axis_types = 13;
const char* nccf_axis_type_name[] = {"", "Lat", "Lon", "GeoX", "GeoY", "GeoZ", 
                                     "Height", "Height", "Pressure", "Time", 
                                     "RadialAzimuth", "RadialElevation", 
                                     "RadialDistance"};


/** \defgroup gs_axis_grp Coordinate axes
    \ingroup gridspec_grp

Coordinate axes are one-dimensional objects with monotonically 
increasing or decreasing values. The dimension of a coordinate 
axis has the same name as the corresponding variable. 

*/

/**
 * \ingroup gs_axis_grp
 * Define a coordinate axis (constructor).
 *
 * \param name name of the object
 * \param len number of elements
 * \param xtype netCDF data type
 * \param data pointer to the data
 * \param standard_name CF standard_name (or NULL if empty)
 * \param units CF units (or NULL if empty)
 * \param cdm_axis_type one of NCCF_NOAXISTYPE, NCCF_LONGITUDE, NCCF_LATITUDE, NCCF_GEOX, NCCF_GEOY, NCCF_GEOZ, NCCF_PRESSURE, NCCF_HEIGHT_UP, NCCF_HEIGHT_DOWN, NCCF_TIME, NCCF_RADAZ, NCCF_RADEL, NCCF_RADDIST
 * \param axis axis attribute (e.g. "X", "Y", "Z", or "T", or NULL)
 * \param positive_up 1 if pointing up, 0 otherwise
 * \param formula_terms optional attribute, can be NULL
 * \param axisid object ID (output)
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer (Tech-X Corp.) and Ed Hartnett (UCAR).
 */

int nccf_def_axis(const char *name, int len, nc_type xtype, 
                  const void *data, const char *standard_name,
                  const char *units, int cdm_axis_type, 
                  const char *axis, int positive_up,  
                  const char *formula_terms, int *axisid) {

  int i;  

  /* allocate objet */
  struct nccf_axis_type *self;
  self = (struct nccf_axis_type *) malloc(sizeof(struct nccf_axis_type));
  
  /* default initialization */
  self->axis_name = NULL;
  self->len = 0;
  self->axisVar = NULL;
  self->data = NULL;
  self->xtype = NC_NAT;

  /* fill in object */

  self->axis_name = strdup(name);
  nccf_varCreate(&self->axisVar, name);
  nccf_varSetDims(&self->axisVar, 1, &len, &name); // dim has same name as var
  /* fill in attributes */
  if (standard_name != NULL) {
    nccf_varSetAttribText(&self->axisVar, CF_ATTNAME_STANDARD_NAME, 
			  standard_name);
  }
  if (units != NULL) {
    nccf_varSetAttribText(&self->axisVar, CF_ATTNAME_UNITS, 
			  units);
  }
  if (cdm_axis_type)
    {
      /* Is the axis type valid? */
      if (cdm_axis_type < 0 || cdm_axis_type > NCCF_RADDIST)
        return NC_EINVAL;
      nccf_varSetAttribText(&self->axisVar, COORDINATE_AXIS_TYPE, 
                            nccf_axis_type_name[cdm_axis_type]);
    }
  if ( (axis && (strcmp(axis, "X") == 0)) || (strcmp(axis, "Y") == 0) || 
      (strcmp(axis, "Z") == 0) || (strcmp(axis, "T") == 0) ) 
    {
      nccf_varSetAttribText(&self->axisVar, CF_AXIS, axis);
    }
  if (positive_up)
    {
      nccf_varSetAttribText(&self->axisVar, CF_POSITIVE, CF_UP);
    }
  else
    {
      nccf_varSetAttribText(&self->axisVar, CF_POSITIVE, CF_DOWN);
    }
  if (formula_terms)
    {
      nccf_varSetAttribText(&self->axisVar, CF_FORMULA_TERMS, formula_terms);
    }

   /* copy the data */
  self->xtype = xtype;
  self->len = len;
  if (xtype == NC_DOUBLE) {
    double *dDouble = (double *) malloc(len*sizeof(double));
    double *d = (double *) data;
    for (i = 0; i < len; ++i) {
      dDouble[i] = d[i];
    }
    self->data = (void *) dDouble;
  }
  else if (xtype == NC_FLOAT) {
    float *dFloat = (float *) malloc(len*sizeof(float));
    float *d = (float *) data;
    for (i = 0; i < len; ++i) {
      dFloat[i] = d[i];
    }
    self->data = (void *) dFloat;
  }
  else if (xtype == NC_INT) {
    int *dInt = (int *) malloc(len*sizeof(int));
    int *d = (int *) data;
    for (i = 0; i < len; ++i) {
      dInt[i] = d[i];
    }
    self->data = (void *) dInt;
  }
  else if (xtype == NC_SHORT) {
    short *dShort = (short *) malloc(len*sizeof(short));
    short *d = (short *) data;
    for (i = 0; i < len; ++i) {
      dShort[i] = d[i];
    }
    self->data = (void *) dShort;
  }
  else {
    return CF_EBADTYPE;
  }
  /* attach the pointer */
  nccf_varSetDataPtr(&self->axisVar, xtype, self->data);

  /* add object to list */
  if (CFLIST_AXIS == NULL) nccf_li_new(&CFLIST_AXIS);
  *axisid = nccf_li_add(&CFLIST_AXIS, self);

  return NC_NOERR;
}
