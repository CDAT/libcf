/*
 * API to perform regridding operations on structured grids in libcf.
 *
 * $Id: nccf_regrid.h 915 2012-01-09 16:58:08Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

#ifndef _NCCF_REGRID_H
#define _NCCF_REGRID_H

#include <libcf_src.h>
#include <cflistitem.h>

/* Names we have picked */
#define CF_REGRID_WEIGHTS            "regrid_weights"
#define CF_REGRID_INDICES            "regrid_indices"
#define CF_REGRID_INSIDE_DOMAIN      "regrid_inside_domain"
#define CF_INDEX_OFFSETS             "flat_index_offsets_from_corner_node"
#define CF_ORI_DIMS                  "original_grid_dimensions"

extern struct CFLISTITEM *CFLIST_REGRID;

struct nccf_regrid_type {

   /* list of boxes and their lo/hi sets of indices */
   struct CFLISTITEM *box_lohi;

   /* NetCDF like variable objects */
   struct nccf_var_obj *weights_stt;
   struct nccf_var_obj *lower_corner_indices_stt;
   struct nccf_var_obj *inside_domain_stt;

   /* original grid id */
   int ori_grid_id;

   /* target grid id */
   int tgt_grid_id;
  
   /* number of space dimensions */
   int ndims;

   /* number of target points (product of the target dims) */
   int ntargets;

   /* number of nodes for given cell (should be 2**ndims) */
   int nnodes;

   /* number of valid target values (number of times cell search succeded) */
   int nvalid;
};


#ifdef __cplusplus
extern "C" {
#endif

  int nccf_def_regrid(int ori_grid_id, int tgt_grid_id, int *regrid_id);

  int nccf_free_regrid(int regrid_id);				  

  int nccf_compute_regrid_weights(int regrid_id, int nitermax, double tolpos);

  int nccf_apply_regrid(int regrid_id, int ori_data_id, int tgt_data_id);

  int nccf_inq_regrid_weights(int regrid_id, const int tgt_indices[], 
                              int ori_nodes[], double weights[]);

  int nccf_get_regrid_weights_pointer(int regrid_id, double **datap);

  int nccf_inq_regrid_ntargets(int regrid_id, int *ntargets);

  int nccf_inq_regrid_nvalid(int regrid_id, int *nvalid);

  int nccf_inq_regrid_nnodes(int regrid_id, int *nnodes);

  int nccf_add_regrid_forbidden(int regrid_id, const int lo[], const int hi[]);

  int nccf_put_regrid(int ncid, int regrid_id);

  int nccf_def_regrid_from_file(const char *filename,
                                int *regrid_id);
  
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NCCF_REGRID_H */
