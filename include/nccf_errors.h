/* $Id: nccf_errors.h 821 2011-09-13 02:07:59Z pletzer $ */
#include <netcdf.h>

#ifndef _NCCF_ERRORS_
#define _NCCF_ERRORS_

#define NCCF_ENOAXISID (-1000) /* Not a valid axis id */
#define NCCF_ENOCOORDID (-1001) /* Not a valid coord id */
#define NCCF_ENOGRIDID (-1002) /* Not a valid grid iad */
#define NCCF_ENODATAID (-1003) /* Not a valid data id */
#define NCCF_ENOREGRIDID (-1004) /* Not a valid regrid id */
#define NCCF_ENOMOSAICID (-1005) /* Not a valid mosaic id */
#define NCCF_ENOHOSTID (-1006) /* Not a valid host id */

#define NCCF_VAROBJCREATE (-1007) /* Error when creating a var obj */
#define NCCF_VAROBJCREATEFROMFILE (-1008) /* Error when creating a var obj from file */
#define NCCF_VAROBJDESTROY (-1009) /* Error when destroying var obj */
#define NCCF_VAROBJSETATTRIB (-1010) /* Error when setting an attribute in var obj */
#define NCCF_VAROBJGETATTRIBPTR (-1011) /* Error when getting an attribute from var obj */
#define NCCF_VAROBJSETDIMS (-1012) /* Error when setting dims in var obj */
#define NCCF_VAROBJGETDIMNAMEPTR (-1013) /* Error when getting dims name pointer in var obj */
#define NCCF_VAROBJSETDATAPTR (-1014) /* Error when setting data pointer in var obj */
#define NCCF_VAROBJSETDATA (-1015) /* Error when setting data in var obj */
#define NCCF_VAROBJGETVARNAMEPTR (-1016) /* Error when getting var name pointer in var obj */
#define NCCF_VAROBJGETDATATYPE (-1017) /* Error when getting data type in var obj */
#define NCCF_VAROBJGETDATAPTR (-1018) /* Error when getting data pointer in var obj */
#define NCCF_VAROBJGETNUMDIMS (-1019) /* Error when getting num dims in var obj */
#define NCCF_VAROBJGETATTRIBLIST (-1020) /* Error when getting attrib list in var obj */
#define NCCF_VAROBJSETVARNAME (-1021) /* Error when setting var name in var obj */
#define NCCF_VAROBJWRITELISTOFVARS (-1022) /* Error occurred when writing var objs */

#define NCCF_EINDEXOUTOFRANGE (-1030) /* Error index out of range */
#define NCCF_EPARSERANGES (-1031) /* Error parsing string ranges/slices */

#define NCCF_ENOTHOSTFILE (-1040) /* Error the file given is not a host file */
#define NCCF_EPUTCOORD (-1041) /* Error occurred when writing coordinate object to file */
#define NCCF_EPUTGRID (-1042) /* Error occurred when writing grid object to file */
#define NCCF_EPUTDATA (-1043) /* Error occurred when writing data object to file */
#define NCCF_EPUTREGRID (-1044) /* Error occurred when writing regrid object to file */
#define NCCF_EPUTMOSAIC (-1045) /* Error occurred when writing mosaic object to file */
#define NCCF_EPUTHOST (-1046) /* Error occurred when writing host object to file */

#define NCCF_ENODATA (-1050) /* Error no data were set in data object */
#define NCCF_EVERTMISMATCH (-1051) /* Vertices mismatch in polytope viewer */
#define NCCF_ENGRIDMISMATCH (-1052) /* ngrid from grid don't match ngrid from data */
#define NCCF_EATTEXISTS (-1053) /* Attribute exists */
#define NCCF_EBADGRIDINDEX (-1054) /* Bad grid/tile index */
#define NCCF_EBADVAR (-1055) /* Bad variable name */
#define NCCF_ENDIMS (-1056) /* Number of dimensions should be > 0 */

#define NCCF_LISTITEMEXISTS ( -1060 ) /* Warning Inserting existing item into a list */


#endif /* _NCCF_ERRORS_ */
