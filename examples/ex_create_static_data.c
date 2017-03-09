/**
 * Create static data on a cube grid mosaic.
 *
 * "$Id: ex_create_static_data.c 767 2011-06-06 23:20:19Z pletzer $"
 *
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

int ex_create_static_data( const char *data_id, const char *coordinates_id,
                      int nCells, int nStatTiles, int nGrids, int gridids[],
                      char **staticFiles, int staticids[] ){

  int iStat, i, j, index, status, ncid, globalid;
  int save = 0;
  int toterr = 0;
  float  data_f[nCells * nCells];
  const char *df_sn = "static_data";
  const char *df_ut = "meters";
  char *name;
  name = ( char* )calloc( STRING_SIZE, sizeof(char));

  for( iStat = 0; iStat < nStatTiles; iStat++){

    /* Create the data */
    for( j = 0; j < nCells; j++ ){
      for( i = 0; i < nCells; i++ ){
        index = i + j * nCells;
        data_f[index] = iStat + i * (float)( j+1 );
      }
    }

    /* Define the data */
    int gid = iStat % nGrids;
    status += nccf_def_data( gridids[gid], NAME_STATIC_DATA_VARN,
                                       df_sn, df_ut,
                                       NULL, &staticids[iStat] );
    status = nccf_inq_grid_name(gridids[gid], name);
    toterr += status;
    status += nccf_set_data_float( staticids[iStat], data_f, 
					  save, NC_FILL_FLOAT);
    toterr += status;

    if(( status = nccf_def_global( &globalid ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_FILETYPE,
                                       CF_GLATT_FILETYPE_STATIC_DATA, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_GRIDNAME,
                                       name, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_DATA_ID,
                                       data_id, 0 ))) ERR;
    if(( status = nccf_add_global_att( globalid, CF_COORDINATES_ID,
                                       coordinates_id, 0 ))) ERR;

    /* Open the file */
    sprintf( staticFiles[iStat], "%s%d.nc", NAME_STATIC_DATA, iStat );
    status += nc_create( staticFiles[iStat], NC_CLOBBER, &ncid );

    /* Add global attributes */
    nc_redef( ncid );
    nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE,
            strlen("Dummy data for GRIDSPEC")+1,
            "Dummy data for GRIDSPEC" );
    nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation",
                 "Multiple grids, variables, time files and times per file",
                 "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
    nc_enddef( ncid );

    /* Write the data to disk */
    status += nccf_put_data(staticids[iStat], ncid);
    status += nccf_put_global(globalid, ncid);
    toterr += status;


    status += nc_close( ncid );
    status += nccf_free_global( globalid );
  }

  free( name );
  return toterr;

}
