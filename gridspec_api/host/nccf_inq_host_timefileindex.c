/* $Id: nccf_inq_host_timefileindex.c 822 2011-09-13 14:39:33Z pletzer $ */

#include <netcdf.h>
#include <nccf_host.h>
#include <string.h>
#include <stdio.h>
#include <netcdf.h>
#include <nccf_errors.h>
#include <nccf_axis.h>
#include <nccf_coord.h>
#include <nccf_data.h>
#include <nccf_grid.h>
#include <nccf_constants.h>

/**
 * \ingroup gs_host_grp
 * Inquire the file index of a time dependent variable
 * 
 * \param hostid the ID for the host object
 * \param varname the variable name
 * \param vfindx file index (output)
 * \return NC_NOERR on success
 * \author Alexander Pletzer and Dave Kindig, Tech-X Corp
 */
int nccf_inq_host_timefileindex(int hostid, const char *varname, 
				int *vfindx) {

  int toterr = 0;
  int status;
  struct nccf_host_type *self;
  self = nccf_li_find(&CFLIST_HOST, hostid);
  
  int ntimedata;
  int gfindx = 0;
  int tfindx = 0;
  int ncid, ifile, ivar, nvars;
  
  char *fname = calloc(STRING_SIZE, sizeof(char));
  char *vname = calloc(STRING_SIZE, sizeof(char));

  *vfindx = -1;
  status = nccf_inq_host_ntimedatafiles(hostid, &ntimedata);
  toterr += abs(status);
  for (ifile = 0; ifile < ntimedata; ++ifile) {
    status = nccf_inq_host_timefilename(hostid, tfindx, 
					ifile, gfindx, fname);
    toterr += abs(status);
    status = nc_open(fname, NC_NOWRITE, &ncid);
    toterr += abs(status);
    status = nc_inq_nvars(ncid, &nvars);
    toterr += abs(status);
    for (ivar = 0; ivar < nvars; ++ivar) {
      status = nc_inq_varname(ncid, ivar, vname);
      toterr += abs(status);
      if ( strcmp(vname, varname) == 0 ) {
	/* found! */
	*vfindx = ivar;
	status = nc_close(ncid);
	toterr += abs(status);
        free(fname);
        free(vname);
	return toterr;
      }
    }
    status = nc_close(ncid);
    toterr += abs(status);
  }

  free(fname);
  free(vname);

  return NCCF_EBADVAR;
}
