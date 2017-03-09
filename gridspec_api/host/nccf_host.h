/*
 * API for host objects
 *
 * "$Id: nccf_host.h 835 2011-09-15 01:00:55Z pletzer $"
 */

#ifndef _NCCF_HOST_H
#define _NCCF_HOST_H

#include <netcdf.h>
#include <libcf_src.h>
#include <cflistitem.h>

#include <nccf_varObj.h>
#include <nccf_mosaic.h>
#include <nccf_data.h>
#include <nccf_grid.h>
#include <nccf_coord.h>
#include <nccf_constants.h>
#include <nccf_handle_error.h>

/* variable names we choose */

#define CF_HOST_MOSAIC_FILENAME      "mosaic_filename"
#define CF_HOST_TILE_NAMES           CF_MOSAIC_TILE_NAMES
#define CF_HOST_TILE_FILENAMES       "tile_filenames"
#define CF_HOST_STATDATA_FILENAME    "static_data_filenames"
#define CF_HOST_TIMEDATA_FILENAME    "time_data_filenames"

extern struct CFLISTITEM *CFLIST_HOST;

struct var_information {
  
  /* Name of variable */
  char *varname;

  /* Grid, Static or Time category */
  char *grid_stat_time_cat;

  /* Variable resides in these files */
  /* Each item is a nccf_varObj for each file. */
  struct CFLISTITEM *file_var_obj;

  /* Number of dimensions */
  int ndims;

  /* Dimension Ids for variable */
  int *dimids;

  /* Number of attributes */
  int natts;

};

struct nccf_host_type {

  /* NetCDF like variable objects */
  struct nccf_var_obj *global;

  struct CFLISTITEM *gridFiles;
  struct CFLISTITEM *gridNames;
  struct CFLISTITEM *timeDataFiles;
  struct CFLISTITEM *statDataFiles;

  /* List of var_information structures */
  /* Each variables item will be of var_information */
  struct CFLISTITEM *variables;
  int nVars;

  /* Permanent data containers for the above nccf_var_obj objects */
  char *mosaicFileBuffer;

  /* if mosaic file was added */
  int hasMosaic;

  /* number of grids (or tiles) */
  int nGrids;

  /* number of static data files (e.g. topography),
     the host files will contain nStatDataFiles*nGrids files */
  int nStatDataFiles;

  /* number of time dependent data files 
     the host file will contain nTimeDataFiles*nGrids*nTimeSlices files */
  int nTimeDataFiles;

  /* Number time slices in an indvidual time dependant variable file */
  int nTimeSlices;

  /* unique identifier binding the mosaci file layout */
  char *uuid;

  /* coordinates_id - Needed to match coordinates_id to other files */
  char *coordinates_id;

  /* data_id - Needed to match data_id to other files */
  char *data_id;

};

#ifdef __cplusplus
extern "C" {
#endif

int nccf_def_host(const char *coordinates_id, const char *data_id, 
                  int nTimeSlices, int *hostid);

int nccf_add_host_file(int hostid, const char *filename, int force);

int nccf_def_host_from_file(const char *filename, int *hostid);

int nccf_free_host(int hostid);

int nccf_put_host(int ncid, int hostid);

int nccf_inq_host_ngrids(int hostid, int *ngrids);

int nccf_inq_host_nstatdatafiles(int hostid, int *nstatdatafiles);

int nccf_inq_host_ntimedatafiles(int hostid, int *ntimedatafiles);

int nccf_inq_host_ntimeslices( int hostid, int *ntimeslices );

int nccf_inq_host_statfilename(int hostid, int vfindx, int gfindx, 
				 char *fname);

int nccf_inq_host_timefilename(int hostid, int tfindx, 
				 int vfindx, int gfindx, 
				 char *fname);

int nccf_inq_host_gridfilename( int hostid, int gfindx, char *filename );

int nccf_inq_host_gridname( int hostid, int gfindx, char *gridname );

int nccf_inq_host_mosaicfilename( int hostid, char *mosaicfilename );

int nccf_def_host_data(int hostid, const char *varname, int gfindx, 
		       int read_data, int *dataid);

int nccf_inq_host_statfileindex(int hostid, const char *varname, 
				   int *vfindx);	     

int nccf_inq_host_timefileindex(int hostid, const char *varname, 
				int *vfindx);	     

/*!@}*/

#ifdef __cplusplus
}
#endif 

#endif /* _NCCF_HOST_H */

