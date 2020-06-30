/**
 *
 * 	@file cuThomasVBatch.h
 *
 * 	@brief cuThomasVBatch header definition.
 *
 * 	@author Ivan Martinez-Perez ivan.martinez@bsc.es
 * 	@author Pedro Valero-Lara   pedro.valero@bsc.es
 *	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

__global__ void cuThomasVBatch( double *U, double *L, double *D, 
						 double *RHS,int BATCH_COUNT, int MAX_SIZE ); 
