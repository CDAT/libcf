/**
 * C object behaving like an associative array
 * 
 * $Id: nccf_keyvalue.h 833 2011-09-14 20:50:14Z pletzer $
 */

#ifndef NCCF_KEYVALUE
#define NCCF_KEYVALUE

#include <netcdf.h>
#include "cflistitem.h"

struct nccf_kv {
  char *key;      // key used to access elements
  nc_type type;   // type of the element
  int nelem;      // number of element stored
  void *val;      // value attached to each element
};

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Constructor
 * \param lst handle to the object
 */
void nccf_kv_new(struct CFLISTITEM **lst);

/**
 * Destructor
 * \param lst handle to the object
 */
void nccf_kv_del(struct CFLISTITEM **lst);

/**
 * Access element by key
 * \param lst handle to the object
 * \param key key 
 * \param nc_type type of the element (output)
 * \param nelem number of elements (output)
 * \param val pointer to the value(s) (output)
 * 
 * \note NULL will be returned for val if the key is not found.
 * The caller does NOT own the pointer. 
 */
void nccf_kv_get_value(struct CFLISTITEM **lst, const char *key, 
		       nc_type *type, int *nelem, const void **val);

/**
 * Add an element, replacing the one already stored if need be.
 * \param lst handle to the object
 * \param key key 
 * \param nc_type type of the element
 * \param nelem number of elements
 * \param val pointer to the value(s)
 * 
 * \note the val values will be copied so caller is responsible 
 * deallocating the memory associated with keys and values. For 
 * char * strings, nelem will not be used. Use nelem = 1 for 
 * scalars. 
 */
void nccf_kv_add(struct CFLISTITEM **lst, const char *key, 
		 nc_type type, int nelem, const void *val);
/**
 * Set iterator to begin.
 * \param lst handle to the object
 */
void nccf_kv_begin(struct  CFLISTITEM **lst);

/**
 * Move the iterator by one step
 * \param lst handle to the object
 * \return 1 if the iterator is valid, 0 otherwise
 */
int nccf_kv_next(struct  CFLISTITEM **lst);

/**
 * Get the current key
 * \param lst handle to the object
 * \param key pointer to the key 
 * 
 * \note caller does NOT own the key.
 */
void nccf_kv_get_key(struct  CFLISTITEM **lst, const char **key);

#ifdef __cplusplus
}
#endif

#endif // NCCF_KEYVALUE
