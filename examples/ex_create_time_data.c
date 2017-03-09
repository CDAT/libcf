/**
 * Create time dependent data on a cube grid mosaic.
 *
 * "$Id: ex_create_time_data.c 776 2011-06-20 17:44:48Z dkindig $"
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

int ex_create_time_data( const char *data_id, const char *coordinates_id,
                         int iVar, struct time_struct *time, int gridids[],
                         char **timeFiles, int timeids[] ){

  int iTile, iTPF, i, j, ind, index, ii, status, ncid, iTime, globalid;
  int time_dimidp, time_varidp;
  int save = 0;
  int toterr = 0;
  long *times;
  float  data_f[time->nCells * time->nCells];
  char *name;
  name = ( char* )calloc( STRING_SIZE, sizeof(char));

  /* Initialize the time data */
  for( j = 0; j < time->nCells; j++ ){
    for( i = 0; i < time->nCells; i++ ){
      ind = i + j * time->nCells;
      data_f[ind] = 0.0;
    }
  }

  float div = time->nGrids + time->nTimes + time->nTimesPerFile +
              time->nCells + time->nCells;
  div = div - 5;
  float div2  =  div / 2.0;

  for( iTime = 0; iTime < time->nTimes; iTime++ ){
    for( iTile = 0; iTile < time->nGrids; iTile++){

      /* Create the data */
      index = iTile + ( time->nGrids * ( iTime + ( time->nTimes * iVar )));

      /* Open the file*/
      sprintf( timeFiles[index], "%s_%d-%d-%d.nc", NAME_TIME_DATA,
               iVar, iTime, iTile );
      status = nc_create( timeFiles[index], NC_CLOBBER, &ncid );
      if( status ) return status;

      /* Add the time coordinate */
      nc_redef( ncid );
      status = nccf_def_time(ncid,
                    "time",
                    time->nTimesPerFile,
                    NC_LONG,
                    time->time_units,
                    time->time_stanname,
                    &time_dimidp,
                    &time_varidp);
      nc_enddef( ncid );

      times = malloc( time->nTimesPerFile * sizeof(long));
      for( ii = 0; ii < time->nTimesPerFile; ii++ )
        times[ii] = ii + time->nTimesPerFile * iTime;
      if( status ) ERR;
      status = nc_put_var_long( ncid, time_varidp, times );
      if( status ) ERR;

      /* Define the time data for each tile */
      status = nccf_def_data( gridids[iTile], time->varname,
                                         time->st_name, time->units,
                                         NAME_TIME_VARN, &timeids[index] );
      toterr += status;

      status = nccf_set_data_float( timeids[index], data_f, 
					   save, NC_FILL_FLOAT );
      toterr += status;

      /* Add global attributes */
      nc_redef( ncid );
      nc_put_att_text( ncid, NC_GLOBAL, CF_TITLE,
              strlen("Dummy data for GRIDSPEC")+1,
              "Dummy data for GRIDSPEC" );
      nccf_def_notes( ncid, NC_GLOBAL, "Tech-X Corporation",
                   "Multiple grids, variables, time files and times per file",
                   "Comment", "https://ice.txcorp.com/trac/modave/wiki/CFProposalGridspec" );
      nc_enddef( ncid );


      /* populate the data */
      for( iTPF = 0; iTPF < time->nTimesPerFile; iTPF++ ){
        for( j = 0; j < time->nCells; j++ ){
          for( i = 0; i < time->nCells; i++ ){
            ind = i + j * time->nCells;
            data_f[ind] = time->value + (( iTile + iTime + iTPF + j + i ) - div2);
          }
        }

        /* Write to disk*/
        status = nccf_put_data(timeids[index], ncid);
        toterr += status;
      }

      /* Global Attributes */
      status = nccf_inq_grid_name(gridids[iTile], name);
      status = nccf_def_global( &globalid );
      status = nccf_add_global_att(globalid, CF_GRIDNAME, name, 0);
      status = nccf_add_global_att(globalid, CF_FILETYPE,
                                   CF_GLATT_FILETYPE_TIME_DATA, 0);
      status = nccf_add_global_att(globalid, CF_DATA_ID, data_id, 0);
      status = nccf_add_global_att(globalid, CF_COORDINATES_ID, coordinates_id, 0);
      status = nccf_put_global(globalid, ncid);

      toterr += status;
      status = nc_close( ncid );
      status = nccf_free_global( globalid );
      if( status ) return status;
      free( times );
    }
  }

  free( name );
  return toterr;

}

