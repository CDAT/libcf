/**
 * API to define a global attribute var_obj for the gridspec convention to
 * libcf.
 *
 * $Id: nccf_global.h 834 2011-09-14 21:21:07Z pletzer $
 */

#ifndef _NCCF_GLOBAL_H
#define _NCCF_GLOBAL_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

extern struct CFLISTITEM *CFLIST_GLOBAL;

/**
 * Global attributes structure
 */
struct nccf_global_type {
  /* Remember the global attribute takes a "" for a variable name
   * Not NC_GLOBAL */

  /* Global attributes */
  struct nccf_var_obj *global;

};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_global(int *globalid);

int nccf_def_global_from_file(const char *filename, int *globalid);

int nccf_free_global(int globalid);

int nccf_put_global(int ncid, int globalid);

int nccf_add_global_att(int globalid, const char *attname, 
                        const char *attvalue, int actionflag);

int nccf_inq_global_att(int globalid, const char *name, char *value);

int nccf_inq_global_natts(int globalid, int *natts);

int nccf_inq_global_attval(int globalid, int attid,
                           char *attname, char *attval);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NCCF_GLOBAL_H */


