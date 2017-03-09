/**
 * Test nccf_keyvalue
 * $Id: tst_nccf_keyvalue.c 664 2011-04-01 17:36:20Z pletzer $
 */

#include <nccf_keyvalue.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

int main() {
  
  int i;
  struct CFLISTITEM *lst;

  //
  // create attribut list 
  //
  nccf_kv_new(&lst);

  //
  // push attributes in
  //
  // add a string attribute, no need to pass the string length
  nccf_kv_add(&lst, "att_string", NC_CHAR, 0, "att_string value");
  // add a float attribute, a scalar
  float vf = 1;
  nccf_kv_add(&lst, "att_float", NC_FLOAT, 1, &vf);
  // add a double array
  double ad[] = {2, 3};
  nccf_kv_add(&lst, "att_double_array", NC_DOUBLE, 2, ad);

  //
  // access known attributes
  //
  int nelem;
  nc_type type;
  const void *ds;
  nccf_kv_get_value(&lst, "att_string", &type, &nelem, &ds);
  assert(type == NC_CHAR);
  assert(strcmp(ds, "att_string value") == 0);
  
  const float *df;
  nccf_kv_get_value(&lst, "att_float", &type, &nelem, 
		    (const void **) &df);
  assert(nelem == 1);
  assert(*df == 1);

  const double *dd;
  nccf_kv_get_value(&lst, "att_double_array", &type, &nelem, 
		    (const void **) &dd);
  assert(nelem == 2);
  assert(dd[0] == 2 && dd[1] == 3);

  // try to access a non-existant key 
  const void *dn;
  nccf_kv_get_value(&lst, "non existant", &type, &nelem, &dn);
  assert(dn == NULL);

  //
  // iterator 
  //
  nccf_kv_begin(&lst);
  while ( nccf_kv_next(&lst) ) {
    const char *k;
    const void *data;
    nccf_kv_get_key(&lst, &k);
    nccf_kv_get_value(&lst, k, &type, &nelem, &data);
    if (type == NC_CHAR) {
      printf("%s -> %s (NC_CHAR)\n", k, (char *) data);
    }
    else if (type == NC_DOUBLE) {
      const double *dd = data;
      printf("%s -> ", k);
      for (i = 0; i < nelem; ++i) {
	printf("%lf ", dd[i]);
      }
      printf(" (NC_DOUBLE) \n");
    }
    else if (type == NC_FLOAT) {
      const float *df = data;
      printf("%s -> ", k);
      for (i = 0; i < nelem; ++i) {
	printf("%f ", df[i]);
      }
      printf("(NC_FLOAT) \n");
    }
    else if (type == NC_INT) {
      const int *di = data;
      printf("%s -> ", k);
      for (i = 0; i < nelem; ++i) {
	printf("%d ", di[i]);
      }
      printf("(NC_INT) \n");
    }
    else {
      // error
      assert(1 == 0);
    }
  }
  
  //
  // destroy attribute list
  //
  nccf_kv_del(&lst);
  return 0;
}
