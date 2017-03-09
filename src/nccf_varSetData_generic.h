/* $Id: nccf_varSetData_generic.h 840 2011-09-23 03:25:52Z pletzer $ */

int nccf_varSetData_TYPE_(struct nccf_var_obj **v, 
			  const _TYPE_ val[]) {
  int status = NC_NOERR;
  int i, ntot;
  _TYPE_ *d;
  (**v).data_type = _NC_TYPE_;
  status = nccf_varGetNumValsPerTime(v, &ntot);
  d = (_TYPE_ *) malloc(sizeof(_TYPE_) * ntot);
  for (i = 0; i < ntot; ++i) {
    d[i] = val[i];
  }
  // set pointer
  (**v).data = (void *) d;
  (**v).save = 1;
  return status;
}
