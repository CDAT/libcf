/**
 * Create a cubed sphere mosaic of nGrids with nCells dimensions.
 *
 * "$Id: ex_create_cube_grids.c 784 2011-07-14 19:53:33Z pletzer $"
 */

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <netcdf.h>

#include "nccf_grid.h"
#include "nccf_global.h"
#include "nccf_data.h"
#include "nccf_mosaic.h"
#include "nccf_host.h"
#include "nccf_utility_functions.h"
#include "nccf_handle_error.h"

#include "examples.h"

int ex_create_cube_grids( const char *coordinates_id, int nCells, int nGrids,
      char **gridFiles, int gridids[] ){

  int status = NC_NOERR;
  int i, iGrid;
  const int ndims = 2;
  int dims[] = {nCells, nCells};
  double lon[dims[0] * dims[1]];
  double lat[dims[0] * dims[1]];
  const int save = 1;
  const char *dimNames[] = {"nlon", "nlat"};
  int coordids[nGrids * 2];
  const char *grid_names[] = {NAME_GRID_0, NAME_GRID_1, NAME_GRID_2};

  const int faceVect0[] = {1, 0, 0};
  const int faceVect1[] = {0, 1, 0};
  const int faceVect2[] = {0, 0, 1};
  const int *faceVectors[] = {faceVect0, faceVect1, faceVect2};

  int globalids[nGrids];

  for (iGrid = 0; iGrid < nGrids; ++iGrid) {
    status = nccf_get_cubedsphere_grid(dims, faceVectors[iGrid], lon, lat);
    if (status != NC_NOERR) ERR;
    status = nccf_def_lon_coord(ndims, dims, dimNames, lon, save, &coordids[0 + iGrid*2]);
    if (status != NC_NOERR) ERR;
    status = nccf_def_lat_coord(ndims, dims, dimNames, lat, save, &coordids[1 + iGrid*2]);
    if (status != NC_NOERR) ERR;
    status = nccf_def_grid(&coordids[iGrid*2], grid_names[iGrid], &gridids[iGrid]);
    if (status != NC_NOERR) ERR;
    if(( status = nccf_def_global( &globalids[iGrid] ))) ERR;
    if(( status = nccf_add_global_att( globalids[iGrid], CF_FILETYPE, CF_GLATT_FILETYPE_GRID, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalids[iGrid], CF_GRIDNAME, grid_names[iGrid], 0 ))) ERR;
    if(( status = nccf_add_global_att( globalids[iGrid], CF_COORDINATES_ID, coordinates_id, 0 ))) ERR;
  }


  /* Create some file names */
  strncpy(gridFiles[0], NAME_GRID_0_FILE, STRING_SIZE);
  strncpy(gridFiles[1], NAME_GRID_1_FILE, STRING_SIZE);
  strncpy(gridFiles[2], NAME_GRID_2_FILE, STRING_SIZE);

  /* write files */
  char scrip_filename[STRING_SIZE];
  for (i = 0; i < nGrids; ++i) {
    sprintf(scrip_filename, "gsex2_grid%d_scrip.nc", i);
    status = nccf_save_grid_scrip(gridids[i], scrip_filename);
    if (status != NC_NOERR) ERR;
    int ncid;
    status = nc_create(gridFiles[i], NC_CLOBBER, &ncid);
    if (status != NC_NOERR) ERR;
    status = nccf_put_grid(gridids[i], ncid);
    if (status != NC_NOERR) ERR;
    status = nccf_put_global(globalids[i], ncid);
    if (status != NC_NOERR) ERR;

    /* Add global attributes */
    nc_redef( ncid );
    nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE,
            strlen("Dummy data for GRIDSPEC")+1,
            "Dummy data for GRIDSPEC" );
    nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation",
                 "Multiple grids, variables, time files and times per file",
                 "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
    nc_enddef( ncid );

    status = nc_close(ncid);
    if (status != NC_NOERR) ERR;
    status = nccf_free_global(globalids[i]);
  }

  return status;
}
