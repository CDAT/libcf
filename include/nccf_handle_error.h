/*
 * "$Id: nccf_handle_error.h 373 2011-01-14 19:53:49Z pletzer $"
 */

#ifndef _NCCF_HANDLE_ERROR_H
#define _NCCF_HANDLE_ERROR_H

#define ERR nccf_handle_error(__FILE__,__LINE__,status)

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void nccf_handle_error(const char *filename, int linenumber, int status);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NCCF_HANDLE_ERROR_H */
