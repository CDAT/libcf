/*
 * API for mosaic objects
 *
 * $Id: nccf_mosaic.h 835 2011-09-15 01:00:55Z pletzer $
 *
 */

#ifndef _NCCF_MOSAIC_H
#define _NCCF_MOSAIC_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>
#include <nccf_varObj.h>
#include <nccf_handle_error.h>
#include <nccf_grid.h>
#include <nccf_coord.h>
#include <nccf_constants.h>

#define CF_MOSAIC_COORDINATE_NAME    "gridspec_coordinate_names"
#define CF_MOSAIC_TILE_NAMES         "tile_names"
#define CF_MOSAIC_TILE_CONTACTS      "tile_contacts"
#define CF_MOSAIC_CONTACT_MAP        "contact_map"

extern struct CFLISTITEM *CFLIST_MOSAIC;

struct nccf_mosaic_type {

  /* Name of object */
  char *name;

  /* NetCDF like variable objects */
  struct nccf_var_obj *coordnames;
  struct nccf_var_obj *gridNames;
  struct nccf_var_obj *gridToGrid;
  struct nccf_var_obj *contactIndex;

  /* Permanent data containers
     for the above nccf_var_obj objects */
  char *coordnamesbuffer;
  char *gridnamesbuffer;
  char *gridtogridbuffer;
  char *contactindexbuffer;

  struct CFLISTITEM *coordnameslist;
  struct CFLISTITEM *gridnameslist;
  struct CFLISTITEM *gridtogridlist;
  struct CFLISTITEM *contactindexlist;

  /* Number of dimensions */
  int ndims;

  /* number of grids (or tiles) */
  int ngrids;

  /* number of mosaic */
  int ncontacts;

  /* Slice Format */
  char *gs_slice_format;

  /* unique grid identifiers (for checking consistency) */
  int *gridids;

};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_mosaic( int ngrids, const int gridids[], const char *name,
                     int *mosaicid);

int nccf_def_mosaic_from_file(const char *filename, const char *name,
      int *mosaicid);

int nccf_free_mosaic(int mosaicid);

int nccf_put_mosaic(int ncid, int mosaicid);

int nccf_set_mosaic_contact( int mosaicid, int ndims, 
				 int gridid0, int gridid1,
				 int grid0_beg_ind[], int grid0_end_ind[], 
				 int grid1_beg_ind[], int grid1_end_ind[]);

int nccf_compute_mosaic_contacts( int mosaicid, const double period[] );

int nccf_inq_mosaic_ngrids(int mosaicid, int *ngrids);

int nccf_inq_mosaic_ndims(int mosaicid, int *ndims);

int nccf_inq_mosaic_coordnames( int mosaicid, char **coordnames );

int nccf_inq_mosaic_gridranges( int mosaicid, int index, 
				 int *gridid0, int *gridid1,
				 int grid0_beg_ind[], int grid0_end_ind[], 
				 int grid1_beg_ind[], int grid1_end_ind[]);

int nccf_inq_mosaic_gridids(int mosaicid, int gridids[]);

int nccf_inq_mosaic_gridname(int mosaicid, int index, char *file);

int nccf_add_mosaic_att( int mosaicid, const char *name,
                          const char *value );

int nccf_inq_mosaic_ncontacts( int ncid, int *ncontacts );

int nccf_inq_mosaic_contactmap(int mosaicid, int index, char *contact_map);

int nccf_inq_mosaic_tilecontact(int mosaicid, int index, char *contact_map);

int nccf_inq_mosaic_tileseparator(char *tile_separator);

int nccf_print_mosaic_as_polytopes(int mosaicid, const char *file_name);

#ifdef __cplusplus
}
#endif

#endif /* _NCCF_MOSAIC_H */


