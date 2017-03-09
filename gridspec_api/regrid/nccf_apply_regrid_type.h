/* $Id: nccf_apply_regrid_type.h 924 2012-03-23 16:47:05Z pletzer $ */

_TYPE *ori_data = NULL;
_TYPE *tgt_data = NULL;
_TYPE datum;
_TYPE old_tgt_value;
double sum_weights;
nc_type xtype;
int k, k2;
const void *fill_value = NULL;
const void *ori_data_fill_value = NULL;

double *weights;
status = nccf_varGetDataPtr(&self->weights_stt, 
                            (void **)&weights);
totErr += abs(status);

int *lower_corner_indices;
status = nccf_varGetDataPtr(&self->lower_corner_indices_stt, 
                            (void **)&lower_corner_indices);
totErr += abs(status);

const int *displ;
status = nccf_varGetAttribPtr(&self->lower_corner_indices_stt,
                              CF_INDEX_OFFSETS, 
                              (const void **)&displ);

const int *ori_dims;
status = nccf_varGetAttribPtr(&self->lower_corner_indices_stt,
                              CF_ORI_DIMS, 
                              (const void **)&ori_dims);
int ori_ntot = 1;
for (i = 0; i < self->ndims; ++i) {
  ori_ntot *= ori_dims[i];
 }
                              
char *inside_domain;
status = nccf_varGetDataPtr(&self->inside_domain_stt, 
                            (void **)&inside_domain);
totErr += abs(status);

status = nccf_get_data_pointer(ori_data_id, &xtype, 
              (void **)&ori_data, &ori_data_fill_value);
totErr += abs(status);
status = nccf_get_data_pointer(tgt_data_id, &xtype, 
              (void **)&tgt_data, &fill_value);
totErr += abs(status);

/* Set the tgt fill_value to the original fill value */
status = nccf_add_data_att(tgt_data_id, "_FillValue", 
                           (const void*)(_TYPE*)ori_data_fill_value);
totErr += abs(status);
    
/* Make sure at least one target point is in the domain */
if (self->nvalid > 0) {
  for (i = 0; i < self->ntargets; ++i) {
    /* Compute sum of weights * field values when inside_domain == 1
       Don't change the value if inside_domain == 0 
       This will set tgt_data[i] to zero if valid and keep 
       tgt_data[i] otherwise
    */
    tgt_data[i] = (_TYPE)(1 - inside_domain[i])*tgt_data[i];

    /* Linear interpolation */
    
    sum_weights = 0;
    old_tgt_value = tgt_data[i];
    for (j = 0; j < nNodes; ++j) {
      k = i*nNodes + j;
      k2 = lower_corner_indices[i] + displ[j];
      // make sure index is valid, weights should be 
      // zero if k2 is outside.
      k2 = (k2 > ori_ntot - 1? lower_corner_indices[i]: k2);
      datum = ori_data[k2];
      tgt_data[i] += (_TYPE)(weights[k]*datum*inside_domain[i]);
      sum_weights += weights[k];
    }

    if (sum_weights == 0) {
      /* cell with invalid data */
      if (ori_data_fill_value) {
        tgt_data[i] = *(_TYPE *)ori_data_fill_value;
      }
      else {
  /* re-instate the old value, the user does not want to use any
     fill_value (== NULL). This amounts to a no-op. */
        tgt_data[i] = old_tgt_value;
      }
    }
  }
  if( totErr ) ERR;
}
