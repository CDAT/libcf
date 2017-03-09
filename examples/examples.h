/**
 * Header for ex2.c. This applies to ex_create_cube_grids.c, 
 * ex_create_static_data.c, and ex_create_time_data.c
 *
 * "$Id: examples.h 784 2011-07-14 19:53:33Z pletzer $"
 */


#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <netcdf.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"
#include "nccf_regrid.h"
#include "nccf_mosaic.h"
#include "nccf_host.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

/* A structure for the time variable */
struct time_struct{
  int nVars;
  int nCells;
  int nTimes;
  int nTimesPerFile;
  int nGrids;

  float value;
  char *varname;
  char *st_name;
  char *units;
  char *time_units;
  char *time_stanname;
  char *time_longname;
};

/* Grid */
#define NAME_GRID_0 "gsex2_grid0"
#define NAME_GRID_1 "gsex2_grid1"
#define NAME_GRID_2 "gsex2_grid2"
#define NAME_GRID_0_FILE "gsex2_grid0.nc"
#define NAME_GRID_1_FILE "gsex2_grid1.nc"
#define NAME_GRID_2_FILE "gsex2_grid2.nc"

/* Static data */
#define NAME_STATIC_DATA_VARN "distance"
#define NAME_STATIC_DATA "gsex2_static_data"

/* Time data */
#define NAME_TIME_DATA   "gsex2_time_data"
#define NAME_TIME_DATA_VARN_v   "V"
#define NAME_TIME_DATA_VARN_t   "T"
#define NAME_TIME_VARN        "time"

/* Host and mosaic */
#define NAME_MOSAIC "gsex2_mosaic"
#define NAME_MOSAIC_FILE "gsex2_mosaic.nc"
#define NAME_HOST_FILE   "gsex2_host.nc"

/**
 * Create nGrids of a cube sphere with nCells resolution.
 *
 * \param coordinates_id The coordinates_id
 * \param nCells The resolution of the grid. 
 * \param nGrids The number of cube-sphere grids. e.g. 3 out of 6
 * \param gridFiles The filenames for each of the grids
 * \param gridids The Id for each grid ( returned )
 */
int ex_create_cube_grids( const char *coordinates_id, int nCells, int nGrids,
      char **gridFiles, int gridids[] );

/**
 * Create a uniform longitude-latitude grid.
 *
 * \param dims The resolutions of the grid
 * \param coord_mins a 2 value tuple for the minimum lon/lat coordinate values
 * \param coord_maxs a 2 value tuple for the maximum lon/lat coordinate values
 * \param coordids The Id for each coordinate (returned)
 * \param gridid The Id of the grid (returned)
 */
int ex_create_lonlat(const int dims[], const double coord_mins[],
                     const double coord_maxs[], int coordids[], 
                     int *gridid);
/**
 * Create nGrids of timeic data dummy data.
 *
 * \param data_id The data_id
 * \param coordinates_id The coordinates_id
 * \param nCells The resolution of the grid. 
 * \param ntimeFiles The number of timeic data variable files
 * \param nGrids The number of cube-sphere grids. e.g. 3 of 6
 * \param gridids The Id for each grid ( returned )
 * \param timeicFiles The file name for each timeic data variable for each grid
 * \param timeicids The Id for each timeic file ( returned )
 */
int ex_create_static_data( const char *data_id, const char *coordinates_id,
      int nCells, int ntimeFiles, int nGrids, int gridids[],
      char **timeicFiles, int timeicids[] );

/**
 * Create nGrids of time data dummy data.
 *
 * \param data_id The data_id
 * \param coordinates_id The coordinates_id
 * \param iVar Which variable is being created
 * \param time Structure of time dependent information. See examples.h for contents
 * \param gridids The Id for each grid ( returned )
 * \param timeFiles The file name for each time data variable for each grid
 * \param timeids The Id for each time file ( returned )
 */
int ex_create_time_data( const char *data_id, const char *coordinates_id,
      int iVar, struct time_struct *time, int gridids[], char **timeFiles,
      int timeids[] );

