/*
 * API to define a structured grid for the gridspec convention to libcf.
 *
 * $Id: nccf_grid.h 891 2011-12-21 17:22:02Z pletzer $
 */

#ifndef _NCCF_STRUCTURED_GRID_H
#define _NCCF_STRUCTURED_GRID_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_constants.h>

extern struct CFLISTITEM *CFLIST_STRUCTURED_GRID;

struct nccf_struct_grid_type {

  /* Coordinate handles */
  int *coordids;

  /* Number of space dimensions */
  int ndims;

  /* Logical name of the grid */
  char *gridname;

  /* (inverse) mask array */
  int *imask;

  /* coordinate periodicity lengths */
  double *coord_periodicity;
};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_grid(const int coordids[], const char *gridname, 
			 int *gridid);

int nccf_def_grid_from_file(const char *filename, int ndims, 
               const char **coordnames, char *gridname, int *gridid);

int nccf_free_grid(int gridid);

int nccf_put_grid(int ncid, int gridid);

int nccf_set_grid_validmask(int gridid, const int *imask);

int nccf_get_grid_mask_pointer(int gridid, int **imask_ptr);

int nccf_save_grid_scrip(int gridid, const char *filename);

int nccf_inq_grid_coordids(int gridid, int coordids[]);

int nccf_inq_grid_ndims(int gridid, int *ndims);

int nccf_inq_grid_coordnames(int gridid, char *coordnames);

int nccf_inq_grid_name(int gridid, char *gridname);

int nccf_inq_grid_periodicity(int gridid, 
                              double coord_periodicities[]);

// private methods
int nccf_detect_grid_periodicity(struct nccf_struct_grid_type *self);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NCCF_STRUCTURED_GRID_H */


