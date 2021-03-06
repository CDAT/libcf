/**
 * $Id: nccf_get_data_pointer.c 738 2011-05-06 22:26:09Z edhartnett $
 */

#include <nccf_data.h>

/**
 * \ingroup gs_data_grp
 * Get pointer to the data.
 *
 * \param dataid data ID
 * \param xtypep (output) pointer to the data type value
 * \param dataPtr (output) pointer to the data
 * \param fill_valuep (output) pointer to the fill (missing) value
 * \return NC_NOERR on success
 *
 * \author Alexander Pletzer and David Kindig, Tech-X Corp.
 */
int nccf_get_data_pointer(int dataid, nc_type *xtypep, 
				 void **dataPtr, const void **fill_valuep){

  int status;
  int totError = 0;
  struct nccf_data_type *self;
  self = nccf_li_find(&CFLIST_STRUCTURED_DATA, dataid);

  int len;
  status = nccf_varInqAttrib(&self->dataVar, "_FillValue", xtypep, &len);
  totError += abs(status);
  status = nccf_varGetAttribPtr(&self->dataVar, "_FillValue", fill_valuep);
  totError += abs(status);
  status = nccf_varGetDataPtr(&self->dataVar, dataPtr);
  totError += abs(status);

  return totError;
}
