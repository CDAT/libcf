/*
 * Test conversion of a flat index to a vector
 * "$Id: tst_nccf_index2vector.c 560 2011-03-14 15:53:54Z pletzer $"
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <nccf_utility_functions.h>

int main(  ){

    int *vector;
    int ndims, i, ndims2;
    int status;

    /* 2D*/
    ndims = 2;
    ndims2 = ndims * 2;
    vector = ( int * )malloc( ndims * sizeof( int ));

    for( i = 0; i < ndims2; i++ ){
        status = nccf_index2vector( i, ndims, vector );
    }

    free( vector );

    /* 3D */
    ndims = 3;
    ndims2 = ndims * 2;
    vector = ( int * )malloc( ndims * sizeof( int ));

    for( i = 0; i < ndims2; i++ ){
        status = nccf_index2vector( i, ndims, vector );
    }

    free( vector );
    
    return 0;
}
