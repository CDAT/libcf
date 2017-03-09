/**
 * Gobals attributes for all files
 * $Id: nccf_constants.h 917 2012-01-27 03:53:01Z pletzer $
 */

#include <math.h> /* for HUGE_VAL */
#include <netcdf.h>
#ifdef HAVE_UUID_H
#include <uuid/uuid.h>
#endif


#ifndef _NCCF_GLOBAL_ATTRS
#define _NCCF_GLOBAL_ATTRS

/* Some numbers */
#define CF_HUGE_INT                  2147483647
#define CF_HUGE_FLOAT                3.402823466e+38f
#define CF_HUGE_DOUBLE               HUGE_VAL

/* HOST FILE DEFINES */
#define CF_GS_HOST_MOSAIC_FILENAME   "gridspec_mosaic_filename"
//#define CF_GS_HOST_TILE_NAMES        "gridspec_tile_names"
#define CF_GS_HOST_TILE_FILENAMES    "gridspec_tile_filenames"
#define CF_GS_HOST_STATDATA_FILENAME "gridspec_static_data_filenames"
#define CF_GS_HOST_TIMEDATA_FILENAME "gridspec_time_data_filenames"
#define CF_GLATT_FILETYPE_HOST       "gridspec_host_file"

/* MOSAIC FILE DEFINES */
#define CF_GS_MOSAIC_CONTACT_MAP     "gridspec_contact_map"
#define CF_GS_MOSAIC_COORDINATE_NAME "gridspec_coord_names"
#define CF_GS_MOSAIC_TILE_NAMES      "gridspec_tile_names"
#define CF_GS_MOSAIC_TILE_CONTACTS   "gridspec_tile_contacts"
#define CF_GLATT_FILETYPE_MOSAIC     "gridspec_mosaic_file"
#define CF_CONTACT_FORMAT            "gridspec_slice_format"
#define CF_INDEX_SEPARATOR           " "
#define CF_TILE_SEPARATOR            " | "
#define CF_RANGE_SEPARATOR           ":"

/* GRID FILE DEFINES */
#define CF_GLATT_FILETYPE_GRID       "gridspec_tile_file"

/* STATIC DATA FILE DEFINES */
#define CF_GLATT_FILETYPE_STATIC_DATA "gridspec_static_data_file"

/* TIME DATA FILE DEFINES */
#define CF_GLATT_FILETYPE_TIME_DATA "gridspec_time_data_file"

/* Our choice for Convenience NETCDF DIMENSION names */
#define CF_DIMNAME_STRING        "nstring"
#define CF_DIMNAME_MOSAIC        "ncontact"
#define CF_DIMNAME_NDIMS         "ndims"
#define CF_DIMNAME_NTIMEDATA     "ntimedata"
#define CF_DIMNAME_NSTATDATA     "nstatdata"
#define CF_DIMNAME_NGRIDS        "ngrid"
#define CF_DIMNAME_NCONTACTS     "ncontacts"
#define CF_DIMNAME_NTARGETS      "ntargets"
#define CF_DIMNAME_NNODES        "nnodes"
#define CF_DIMNAME_NTIMES        "ntimes"

/* NETCDF VARIABLE ATTRIBUTE NAMES */
#define CF_ATTNAME_CF_TYPE_NAME  "gridspec_type_name"
#define CF_ATTNAME_STANDARD_NAME "standard_name"
#define CF_ATTNAME_LONG_NAME     "long_name"
#define CF_ATTNAME_UNITS         "units"
#define CF_ATTNAME_BOUNDS        "bounds"
#define CF_ATTNAME_AXIS          "axis"
#define CF_ATTNAME_POSITIVE      "positive"

/* GLOBAL ATTRIBUTE NAMES */
#define CF_GLOBAL                ""  // Must be empty
#define CF_NA                    "N/A"
#define CF_COORDINATES_ID        "gridspec_coordinates_id"
#define CF_DATA_ID               "gridspec_data_id"
#define CF_FILETYPE              "gridspec_file_type"
#define CF_GRIDNAME              "gridspec_tile_name"
#define CF_ORIGINALFILENAME      "originalfilename"
#define CF_GLOBAL_TITLE          "Title"
#define CF_GLOBAL_IDENTIFIER     "Identifier"
#define CF_GLOBAL_INSTITUTION    "Institution"
#define CF_GLOBAL_MODEL          "Model"
#define CF_GLOBAL_RUN            "Run"

#define CF_FILENAME_GRID         "_grid"
#define CF_FILENAME_CONTACTS     "_contacts"
#define CF_FILENAME_MOSAIC       "_mosaic"

#define STRING_SIZE NC_MAX_NAME

#endif /* _NCCF_GLOBAL_ATTRS */
