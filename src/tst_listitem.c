/**
 * Unit test for list item API
 * $Id: tst_listitem.c 924 2012-03-23 16:47:05Z pletzer $
 */

/* std includes */
#include <stdio.h>
#include <stdlib.h>
#include <nc_tests.h>


/* some struct */
struct X {
  float member;
};

/* Comparison function for inserting data.
 * The user knows the data type and can therefore test equality */
int item_comparision( const void *data1, const void *data2 ){
  int result;
  struct X *d1;
  struct X *d2;

  d1 = (struct X *)data1;
  d2 = (struct X *)data2;

  if( d1->member <  d2->member ) result = -1;
  if( d1->member >  d2->member ) result =  1;
  if( d1->member == d2->member ) result =  0;
  return result;
}

#include <cflistitem.h>

int 
main(int argc, char **argv) 
{

  int id, i;
  struct CFLISTITEM *lst;
  
  /* elements */
  struct X *x0;
  struct X *x1;
  struct X *x2;
  struct X *x3;
  struct X *x4;
  struct X *x5;
  struct X *x;

  /* create elements */
  x0 = (struct X *) malloc(sizeof(struct X)); 
  x0->member = 0;
  x1 = (struct X *) malloc(sizeof(struct X)); 
  x1->member = 1;
  x2 = (struct X *) malloc(sizeof(struct X)); 
  x2->member = 2;
  x3 = (struct X *) malloc(sizeof(struct X)); 
  x3->member = 3;
  x4 = (struct X *) malloc(sizeof(struct X)); 
  x4->member = 4;

  const int outId1[] = {0, 1, 2};
  const int outId2[] = {0, 1, 3, 2};
  const int outId3[] = {4, 0, 1, 3, 2};
  const int outId4[] = {5, 4, 0, 1, 3, 2};
  const float outMem1[] = {0, 1, 2};
  const float outMem2[] = {0, 1, 3, 2};
  const float outMem3[] = {4, 0, 1, 3, 2};
  const float outMem4[] = {1.5, 4, 0, 1, 3, 2};
  
  /* create linked list */
  nccf_li_new(&lst);
  assert( nccf_li_get_nelem( &lst ) == 0 );
  
  /* attach elements */
  id = nccf_li_add(&lst, x0);
  assert(id == 0);

  id = nccf_li_add(&lst, x1);
  assert(id == 1);

  id = nccf_li_add(&lst, x2);
  assert(id == 2);

  /* iterate over elements */
  int index = 0;
  nccf_li_begin(&lst);
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_find(&lst, id);
    assert(outId1[index] == id);
    assert(x->member == outMem1[index++]);
  }

  assert( nccf_li_get_nelem( &lst ) == 3 );

  /****** Insert a new element ******/
  nccf_li_begin(&lst);
  while( nccf_li_next( &lst ) ){
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_find(&lst, id);
    if( id == 1 ){
      id = nccf_li_insert_after( &lst, x3, id);
    }
    x = (struct X *) nccf_li_find(&lst, id);
    assert(x->member == id);
  }
  assert( nccf_li_get_nelem( &lst ) == 4 );
  /* Check list */
  nccf_li_begin(&lst);
  index = 0;
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_find(&lst, id);
    assert(x->member == id);
    assert(x->member == outMem2[index]);
    assert(outId2[index++] == id);
  }
  /* Test insert at beginning of list */
  id = -1;
  id = nccf_li_insert_after( &lst, x4, id);
  x = (struct X *) nccf_li_find(&lst, id);
  assert(x->member == id);
  /* Check list - again */
  index = 0;
  nccf_li_begin(&lst);
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_find(&lst, id);
    assert(x->member == outMem3[index]);
    assert(outId3[index] == id);
    assert(x->member == id);
    index++;
  }
  /* Insert using comparison function pointer */
  x5 = (struct X *) malloc(sizeof(struct X)); 
  x5->member = 1.5;
  id = nccf_li_insert( &lst, x5, &item_comparision, 0 );
  x = (struct X *) nccf_li_find(&lst, id);
  assert( x->member == 1.5 );
  assert( id == 5 );
  /* Check list - again */
  index = 0;
  nccf_li_begin(&lst);
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_find(&lst, id);
    assert(x->member == outMem4[index]);
    assert(outId4[index++] == id);
  }

  /* remove elements */
  x = nccf_li_remove(&lst, 0);
  free(x); x = NULL;
  x = nccf_li_remove(&lst, 2);
  free(x); x = NULL;
  x = nccf_li_remove(&lst, 1);
  free(x); x = NULL;
  x = nccf_li_remove(&lst, 3);
  free(x); x = NULL;
  x = nccf_li_remove(&lst, 4);
  free(x); x = NULL;
  x = nccf_li_remove(&lst, 5);
  free(x); x = NULL;
  nccf_li_del(&lst);

  /* start over again */
  nccf_li_new(&lst);
  
  /* check that iterating over an empty list is ok */
  nccf_li_begin(&lst);
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
  }  

  /* checking adding/removing works */
  x0 = (struct X *) malloc(sizeof(struct X));
  id = nccf_li_add(&lst, x0);
  x = nccf_li_remove(&lst, id);
  free(x); x = NULL;
  nccf_li_del(&lst);

  /* try to add an item to a list that has been reset */
  nccf_li_new(&lst);
  for (i = 0; i < 10; ++i) {
    x = (struct X *) malloc(sizeof(struct X));
    x->member = i;
    id = nccf_li_add(&lst, x);
  }
  /* now reset */
  nccf_li_begin(&lst);
  /* add another element */
  x = (struct X *) malloc(sizeof(struct X));
  x->member = 10;
  /* now add */
  id = nccf_li_add(&lst, x);
  assert(id == 10);

  /* this is the proper way to destroy the list:
     Iterate over all elements and remove each element.
     The head will be removed with the destructor call.
  */
  nccf_li_begin(&lst);
  while (nccf_li_next(&lst)) {
    id = nccf_li_get_id(&lst);
    x = (struct X *) nccf_li_remove(&lst, id);
    free(x); x = NULL;
  }
  nccf_li_del(&lst);

  /* check that the head will be removed even when there
     are no elements in the list
   */
  nccf_li_new(&lst);
  nccf_li_del(&lst);

  SUMMARIZE_ERR;
  FINAL_RESULTS;
}


