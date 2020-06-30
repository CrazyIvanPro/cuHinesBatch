/**
 *
 * 	@file cuThomasVBatch.cu
 *
 * 	@brief cuThomasVBatch.
 *
 * 	@author Ivan Martinez-Perez ivan.martinez@bsc.es
 * 	@author Pedro Valero-Lara   pedro.valero@bsc.es 
 *	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

#include "cuThomasVBatch.h"
#include <cstdio>


__global__ void cuThomasVBatch(
            double *U, double *L, double *D, double *RHS,
			int BATCH_COUNT, int MAX_SYZE ) 
{

    int tid = threadIdx.x + blockDim.x*blockIdx.x;
			
    if(tid < BATCH_COUNT) {

        int first = tid;
        int last  = BATCH_COUNT*(MAX_SYZE-1)+tid;
			
		U[first] /= D[first];
        RHS[first] /= D[first];

        for (int i = first + BATCH_COUNT; i < last; i+=BATCH_COUNT) {
            U[i] /= D[i] - L[i] * U[i-BATCH_COUNT];
            RHS[i] = ( RHS[i] - L[i] * RHS[i-BATCH_COUNT] ) / 
						( D[i] - L[i] * U[i-BATCH_COUNT] );
        }

        RHS[last] = ( RHS[last] - L[last] * RHS[last-BATCH_COUNT] ) / 
						( D[last] - L[last] * U[last-BATCH_COUNT] );

        for (int i = last-BATCH_COUNT; i >= first; i-=BATCH_COUNT) {
            RHS[i] -= U[i] * RHS[i+BATCH_COUNT];
        }
	}
}
