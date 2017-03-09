/**
 * API for adding, removing, and extracting elements from a linked list
 * $Id: cflistitem.h 835 2011-09-15 01:00:55Z pletzer $
 *
 * \author Alexander Pletzer, Tech-X Corp.
 */

#ifndef _CFLISTITEM_H
#define _CFLISTITEM_H

#define CF_LIST_ITEM_HEAD_ID -1


#define _DATATYPE void

struct CFLISTITEM {

  /* a pointer to the head to the linked list 
     there is at least one element in the list (the head) 
   */
  struct CFLISTITEM *first;

  /* a pointer to the next element in the linked list */
  struct CFLISTITEM *next;

  /* a unique id number for each element */
  int id;

  /* pointer to the data 
     USERS MUST #define _DATATYPE BEFORE #include'ing this file
   */
  _DATATYPE *data;
};

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create new list, with a head
 *
 * \param lst list object
 */
void 
nccf_li_new(struct CFLISTITEM **lst);

/**
 * Destroy head of list. Assumes all elements have been removed.
 *
 * \param lst list object
 *
 * \note Call nccf_li_remove to remove each element in turn.
 */
void 
nccf_li_del(struct CFLISTITEM **lst);

/**
 * Reset the list to the beginning
 *
 * \param lst list object
 */
void 
nccf_li_begin(struct CFLISTITEM **lst);

/** 
 * Go to the next element (if it is not the last element)
 *
 * \param lst list object
 * \return 0 if it is the last element, otherwise return 1
 */
int 
nccf_li_next(struct CFLISTITEM **lst);

/** 
 * Get the id of the current element pointer by the linked list
 *
 * \param lst list object
 * \return id
 */
int 
nccf_li_get_id(struct CFLISTITEM **lst);

/**
 * Add an element to the list
 * 
 * \param lst list object
 * \param data pointer to the data to add
 * \return new id
 */
int
nccf_li_add(struct CFLISTITEM **lst, const _DATATYPE *data);


/**
 * Find the number of elements in the list - also is the maximum id.
 * 
 * \param lst list object
 * \return nelem - number of elements in list ( maximum id )
 */
int
nccf_li_get_nelem( struct CFLISTITEM **lst );

/**
 * Insert an element into the list after element with given id.
 * 
 * \param lst list object
 * \param data pointer to the data to add
 * \param id Id to insert element after. Use id = CF_LIST_ITEM_HEAD_ID to 
 *           insert at the beginning.
 * \return new Id of the inserted element. 
 */
int
nccf_li_insert_after(struct CFLISTITEM **lst, const _DATATYPE *data, 
		     const int id);

/**
 * Insert an element into a list. The element will be inserted before 
 * encountering the first list element that is larger than the 
 * inserted element.
 * 
 * \param lst list object
 * \param data pointer to the data to add
 * \param comparison User provided, external comparison function, should 
 *        return:
 *                       -1 if data < data2  
 *                        0 if data == data2
 *                        1 if data > data2
 * \param ifresult0 Boolean = 0, don't add item,  = 1, Add item to list
 * \return new Id of the inserted element
 */
int
nccf_li_insert( struct CFLISTITEM **lst, 
                const _DATATYPE *data,
                int (*comparison)(const _DATATYPE *data, const _DATATYPE *data2), 
                int ifresult0 );

/**
 * Remove element id
 *
 * \param lst list object
 * \param id id of element to remove
 * \return pointer to the data of the removed element
 *
 * \note To destroy the list, remove each element. Then call 
 *       nccf_li_del to remove the head.
 */
_DATATYPE *
nccf_li_remove(struct CFLISTITEM **lst, int id);

/**
 * Get data from element with given id, or NULL if no id was found
 *
 * \param lst list object
 * \param id  id of searched element
 * \return pointer to the data of the element
 */
_DATATYPE *
nccf_li_find(struct CFLISTITEM **lst, int id);

#ifdef __cplusplus
}
#endif

#endif /* _CFLISTITEM_H */
