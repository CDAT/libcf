/**
 * $Id: tst_nccf_varObj.c 922 2012-03-22 17:13:41Z pletzer $
 *
 * Test nccf_var_obj
 */

// std includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <netcdf.h>
// libcf
#include <cflistitem.h>

// tools
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

#include <assert.h>

int main() {

  int status;
  struct nccf_var_obj *ta, *u, *strVar, *global, *global_append;
  struct nccf_var_obj *u2, *ta2, *global2, *u3;
  double *taData;
  double *uData;
  const int ncontacts = 3;
  const int nd[] = {ncontacts, NC_MAX_NAME};
  const char *ndn[] = {"ncontacts", "string"};
  const int ndims = 2;
  int ndims2;
  int *dims2;
  nc_type dataType2;

  const int ni = 10;
  const int nj = 12;
  const int dims[] = {ni, nj};
  int i;
  int ncid;
  const char *dimname[] = {"lon", "lat"};
  char *names;
  names = (char *) calloc( ncontacts, STRING_SIZE * sizeof( char ));

  const short as[] = {1, 2, 3, 4};
  const int ai[] = {1, 2, 3, 4};
  const float af[] = {1.0f, 2.0f, 3.0f, 4.0f};
  const double ad[] = {1.0, 2.0, 3.0, 4.0};
  /* create variables */
  nccf_varCreate(&ta, "ta");
  nccf_varSetAttribText(&ta, "long_name", "air temperature");
  nccf_varSetAttribText(&ta, "units", "degrees K");
  nccf_varSetAttribShortArray(&ta, "as", 4, as);
  nccf_varSetAttribIntArray(&ta, "ai", 4, ai);
  nccf_varSetAttribFloatArray(&ta, "af", 4, af);
  nccf_varSetAttribDoubleArray(&ta, "ad", 4, ad);
  nccf_varSetDims(&ta, ndims, dims, dimname);

  nccf_varCreate(&u, "u");
  /* set attributes of a variety of types */
  nccf_varSetAttribText(&u, "long_name", "zonal wind");
  nccf_varSetAttribText(&u, "units", "m/s");
  nccf_varSetAttribDouble(&u, "u_attr_double", 1);
  nccf_varSetAttribFloat(&u, "u_attr_float", 2);
  nccf_varSetAttribInt(&u, "u_attr_int", 3);
  nccf_varSetAttribShort(&u, "u_attr_short", 4);
  nccf_varSetDims(&u, ndims, dims, dimname);

  int *u_attr_int_value_p;
  nccf_varGetAttribPtr(&u, "u_attr_int",
                       (const void**) &u_attr_int_value_p);
  assert(*u_attr_int_value_p == 3);

  /* check that the attributes are properly stored */
  nc_type xtype;
  int len;
  nccf_varInqAttrib(&u, "u_attr_double", &xtype, &len);
  assert(xtype == NC_DOUBLE);
  assert(len == 1);

  nccf_varCreate(&strVar, "nametest");
  nccf_varSetAttribText(&strVar, "cf_type_name", "name");
  nccf_varSetDims(&strVar, 2, nd, ndn);

  nccf_varCreate(&global, ""); // global attribute
  nccf_varSetAttribText(&global, "history", "a short history");
  
  nccf_varCreate(&global_append, ""); // global attribute
  nccf_varSetAttribText(&global_append, "history", "of human kind");
  
  /* define a variable with unlimited dimension */
  struct nccf_var_obj *timeDependentData;
  nccf_varCreate(&timeDependentData, "timeDependentData"); 
  const int timeDependentData_dims[] = {NC_UNLIMITED, ni};
  const char *timeDependentData_nms[] = {"ntime", "nx"};
  nccf_varSetDims(&timeDependentData, 2, timeDependentData_dims, timeDependentData_nms);
  float *timeDependentData_vals = (float *) malloc(ni * sizeof(float));
  for (i = 0; i < ni; ++i) {
    timeDependentData_vals[i] = i;
  }
  nccf_varSetDataPtr(&timeDependentData, NC_FLOAT, timeDependentData_vals);
  
  /* allocate data and fill in values */
  taData = (double *) malloc(sizeof(double) * ni * nj);
  uData = (double *) malloc(sizeof(double) * ni * nj);

  for (i = 0; i < ni * nj; ++i) {
    taData[i] = 192.0;
    uData[i] = i;
  }
  for (i = 0; i < ncontacts; ++i) {
    sprintf(&names[i * STRING_SIZE], "contact %d", i);
  }

  /* set data pointers */
  nccf_varSetDataPtr(&ta, NC_DOUBLE, taData);
  nccf_varSetDataPtr(&u, NC_DOUBLE, uData);
  nccf_varSetDataPtr(&strVar, NC_CHAR, names);

  /* write to disk */
  if ( (status = nc_create("tst_varObj.nc", NC_CLOBBER, &ncid)) ) ERR;
  nccf_writeListOfVars(ncid, 5, global, global_append, strVar, ta, u);

  /* write first record of time dependent data */
  nccf_writeListOfVars(ncid, 1, timeDependentData);
  for (i = 0; i < ni; ++i) {
    timeDependentData_vals[i] = 2*i;
  }
  /* use this to write subsequent records */
  nccf_writeListOfVarData(ncid, 1, timeDependentData);
  if ( (status = nc_close(ncid)) ) ERR;

  /* Read from a variable */ 
  const char *varname;
  nccf_varGetVarNamePtr( &ta, &varname );

  double *tat;
  nccf_varGetDataPtr( &ta, (void **)&tat );

  int *new_dims;
  nccf_varGetDimsPtr( &ta, &new_dims );

  // read from file 
  if ( (status = nc_open("tst_varObj.nc", NC_NOWRITE, &ncid)) ) ERR;

  /* With data, cast data as doubles */
  nccf_varCreateFromFile(&u2, "u", ncid, 1, 1);
  /* Without data */
  nccf_varCreateFromFile(&u3, "u", ncid, 0, 0);

  /* Create time dependent data, without reading the data */
  struct nccf_var_obj *timeDependentData2;
  nccf_varCreateFromFile(&timeDependentData2, "timeDependentData", ncid, 0, 0);
  /* Check the number of written time slices */
  int *tdims;
  nccf_varGetDimsPtr(&timeDependentData2, &tdims);
  /* Get the last data record */
  nccf_varReadData(&timeDependentData2, ncid, tdims[0]-1, 0);
  float *fd;
  nccf_varGetDataPtr(&timeDependentData2, (void **) &fd);
  for (i = 0; i < ni; ++i) {
    assert(fd[i] == (float)2*i);
  }
 
  nccf_varGetNumDims(&u2, &ndims2);
  assert(ndims2 == ndims);
  nccf_varGetNumDims(&u3, &ndims2);
  assert(ndims2 == ndims);

  nccf_varGetDimsPtr(&u2, &dims2);
  for (i = 0; i < ndims; ++i) {
    assert(dims2[i] == dims[i]);
  }
  nccf_varGetDimsPtr(&u3, &dims2);
  for (i = 0; i < ndims; ++i) {
    assert(dims2[i] == dims[i]);
  }

  nccf_varGetDataType(&u2, &dataType2);
  assert(dataType2 == NC_DOUBLE);

  double *uData2;
  nccf_varGetDataPtr(&u2, (void **) &uData2);
  for (i = 0; i < ni * nj; ++i) {
    assert(uData2[i] == (double) uData[i]);
  }
  nccf_varGetDataPtr(&u3, (void **) &uData2);
  assert(uData2 == NULL); 

  nccf_varAttribIterBegin(&u2);
  while ( nccf_varAttribIterNext(&u2) ) {
    const char *attrName;
    nccf_kv_get_key(&(u2->attr), &attrName);
    int nelem;
    nc_type type;
    const char *attrVal;
    nccf_kv_get_value(&(u2->attr), attrName, &type, &nelem, 
		      (const void **) &attrVal);
    if (strcmp(attrName, "name") == 0) {
      assert(strcmp(attrVal, "u") == 0);
    }
    if (strcmp(attrName, "long_name") == 0) {
      assert(strcmp(attrVal, "zonal wind") == 0);
    }
    if (strcmp(attrName, "units") == 0) {
      assert(strcmp(attrVal, "m/s") == 0);
    }
  }

  nccf_varCreateFromFile(&ta2, "ta", ncid, 1, 0);

  nccf_varInqAttrib(&ta2, "as", &xtype, &len);
  assert(xtype == NC_SHORT);
  assert(len == 4);
  short *as2;
  nccf_varGetAttribPtr(&ta2, "as", (const void **)&as2);
  for (i = 0; i < (int)len; ++i) {
    assert(as2[i] == as[i]);
  }

  nccf_varInqAttrib(&ta2, "ai", &xtype, &len);
  assert(xtype == NC_INT);
  assert(len == 4);
  int *ai2;
  nccf_varGetAttribPtr(&ta2, "ai", (const void **)&ai2);
  for (i = 0; i < (int)len; ++i) {
    assert(ai2[i] == ai[i]);
  }

  nccf_varInqAttrib(&ta2, "af", &xtype, &len);
  assert(xtype == NC_FLOAT);
  assert(len == 4);
  float *af2;
  nccf_varGetAttribPtr(&ta2, "af", (const void **)&af2);
  for (i = 0; i < (int)len; ++i) {
    assert(af2[i] == af[i]);
  }

  nccf_varInqAttrib(&ta2, "ad", &xtype, &len);
  assert(xtype == NC_DOUBLE);
  assert(len == 4);
  double *ad2;
  nccf_varGetAttribPtr(&ta2, "ad", (const void **)&ad2);
  for (i = 0; i < (int)len; ++i) {
    assert(ad2[i] == ad[i]);
  }

  nccf_varCreateFromFile(&global2, "", ncid, 1, 0);

  nccf_kv_begin(&(global2->attr));
  while ( nccf_kv_next(&(global2->attr)) ) {
    const char * attrName;
    nccf_kv_get_key(&(global2->attr), &attrName);
    int nelem;
    nc_type type;
    const char * attrVal;
    nccf_kv_get_value(&(global2->attr), attrName, &type, &nelem, 
		      (const void **) &attrVal);
    if (strcmp(attrName, "history") == 0) {
      assert(strcmp(attrVal, "a short history of human kind") == 0);
    }
  }

  /* close the file */
  if ( (status = nc_close(ncid)) ) ERR;

  /* reclaim memory */
  free(uData);
  free(taData);
  free(names);
  free(timeDependentData_vals);

  nccf_varDestroy(&global);
  nccf_varDestroy(&global_append);
  nccf_varDestroy(&u);
  nccf_varDestroy(&ta);
  nccf_varDestroy(&strVar);
  nccf_varDestroy(&global2);
  nccf_varDestroy(&u2);
  nccf_varDestroy(&u3);
  nccf_varDestroy(&ta2);
  nccf_varDestroy(&timeDependentData);
  nccf_varDestroy(&timeDependentData2);
  return 0;
}
