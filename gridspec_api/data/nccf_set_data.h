/**
 * $Id: nccf_set_data.h 770 2011-06-08 03:08:35Z pletzer $
 */

// generic data setter written with the idea that the preprocessor
// will fill in anything that starts with _ and ends with _

  int status = 0;
  int totError = 0;
  _type_ * dataBuff;

  struct nccf_data_type *self;
  self = nccf_li_find(&CFLIST_STRUCTURED_DATA, dataid);

  /* set the netcdf data type */
  self->dataType = _nc_type_;

  /* set the default fill value */
  _nccf_varSetAttribTYPE_(&self->dataVar, "_FillValue", fill_value);

  nc_type dataType;
  nccf_varGetDataType(&self->dataVar, &dataType);
  if (!save) {
    self->data = (void *) data;
    nccf_varSetDataPtr(&self->dataVar, _nc_type_, self->data);
  }
  else if(dataType == NC_NAT || dataType == _nc_type_) {
    /* make sure the type is not changing */
    /* number of space dimensions */
    int ndims;
    status = nccf_inq_grid_ndims(self->gridid, &ndims);
    totError += abs(status);
    /* dimensions */
    int coordIds[ndims];
    status = nccf_inq_grid_coordids(self->gridid, coordIds);
    totError += abs(status);
    int dims[ndims];
    status = nccf_inq_coord_dims(coordIds[0], dims);
    totError += abs(status);
    /* compute tot number of elements */
    int ntot = 1;
    int i;
    for (i = 0; i < ndims; ++i) {
      ntot *= dims[i];
    }
    /* allocate and copy the data */
    dataBuff = (_type_ *) malloc( ntot * sizeof(_type_) );
    for (i = 0; i < ntot; ++i) {
      dataBuff[i] = data[i];
    }
    self->data = (void *) dataBuff;
    nccf_varSetDataPtr(&self->dataVar, _nc_type_, self->data);
    self->save = 1;
  }
  else {
    /* no op, cannot reset data of a different type */
    totError += 1;
    return totError;
  }

  return totError;

