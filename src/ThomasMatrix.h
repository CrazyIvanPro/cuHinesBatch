/**
 *
 *  @file ThomasMatrix.cu
 *
 *  @brief cuThomasBatch kernel implementaion.
 *
 *  @author Ivan Martinez-Perez ivan.martinez@bsc.es
 *  @author Pedro Valero-Lara   pedro.valero@bsc.es
* 	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <sstream>
#include "memm.h"


template<typename T>
void copyContainer(std::vector<T> &start,T* &end ){

    end = (T*) malloc(start.size()*sizeof start[0]);
    std::copy(start.begin(),start.end(),end);

}


double fRand(double fMin, double fMax){

    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);

}

ThomasMatrix loadThomasMatrixSyn(int size){

    ThomasMatrix tm;

    std::vector<double> u;
    std::vector<double> l;
    std::vector<double> d;
    std::vector<double> rhs;

    for (int i = 0; i < size; ++i)
    {
      
        u.push_back(fRand((double)-2.,(double)2.));
        l.push_back(fRand((double)-2.,(double)2.));
        d.push_back(fRand((double)5.,(double)10.));
        rhs.push_back(fRand((double)-2.,(double)2.));
      
    }

    copyContainer(u,tm.a);
    copyContainer(l,tm.b);
    copyContainer(d,tm.d);
    copyContainer(rhs,tm.rhs);
    tm.batchCount=1;
    tm.M=size;

    return tm;

}

ThomasMatrix loadThomasFromVectors(double *D, double *U, double *L, double *RHS, int size){
	
	ThomasMatrix tm;
	
	tm.a = (double*) malloc(size*sizeof(double));
	tm.b = (double*) malloc(size*sizeof(double));
	tm.d = (double*) malloc(size*sizeof(double));
	tm.rhs = (double*) malloc(size*sizeof(double));
	
	for(int i=0;i < size-1; ++i){
		tm.a[i] = U[i];
		tm.b[i] = L[i];
		tm.d[i] = D[i];
		tm.rhs[i] = RHS[i];
	}
	tm.d[size-1] = D[size-1];
	tm.rhs[size-1] = RHS[size-1];
	
	tm.batchCount=1;
	tm.M = size;

	return tm;
}


// Square tridiagonal matrices 
void tridiagonal_CSR_to_vector( int *CSR_ROW_PTR_A, int ROW_SIZE, 
						  int *CSR_COL_IND_A, int COL_SIZE, 
						  double *CSR_VAL_A, double  *D, double  *U, 
						  double  *L,int SIZE )
{

	printf("CSR_TO_VECTOR\n");
	printf("row_size : %d\n",ROW_SIZE);
	int i,j;
	int u_index,l_index,d_index=0;
	
	for(i=0;i<ROW_SIZE-1;++i){
		printf("Row %d from %d to %d\n",i,CSR_ROW_PTR_A[i],CSR_ROW_PTR_A[i+1]);
		for(j=CSR_ROW_PTR_A[i];j<CSR_ROW_PTR_A[i+1];++j){
			// digonal => row == column
			printf("row_index: %d\n",i);
			printf("row_ptr: %d\n",j);
			printf("column_index: %d\n",CSR_COL_IND_A[j-1]);
			printf("Val: %f\n",CSR_VAL_A[j-1]);
			
			if(CSR_COL_IND_A[j-1] == i+1){
				D[d_index]=CSR_VAL_A[j-1];
				++d_index;
				//printf("filling D\n");
			}
			// lower => row > column
			else if(CSR_COL_IND_A[j-1] < i+1){
				L[l_index]=CSR_VAL_A[j-1];
				++l_index;
				//printf("filling L\n");
			}
			// upper => row < column
			else if(CSR_COL_IND_A[j-1] > i+1){
				U[u_index]=CSR_VAL_A[j-1];
				++u_index;
				//printf("filling U\n");
			}
			else printf("ERROR: NOT TRIDIAGONAL MATRIX\n");
		}
	}	
}

void tridiagonal_vector_to_CSR( int *CSR_ROW_PTR_A, int ROW_SIZE, 
						  int *CSR_COL_IND_A, int COL_SIZE, 
						  double *CSR_VAL_A, double  *D, double  *U, 
						  double  *L,int SIZE )
{
	int i,valor;
	valor=3;
	
	//compute CSR rows values
	CSR_ROW_PTR_A[0]=1;
	for(i=1;i<ROW_SIZE-1;++i){
		CSR_ROW_PTR_A[i]=valor;
		valor+=3;
	}
	CSR_ROW_PTR_A[ROW_SIZE-1] = CSR_ROW_PTR_A[ROW_SIZE-2] + 2;

	CSR_COL_IND_A[0] = 1;
	CSR_COL_IND_A[1] = 2;
	
	//compute CSR colum values
	valor=1;
	for(i=2;i<COL_SIZE-3;i+=3){
		CSR_COL_IND_A[i] = valor;
		CSR_COL_IND_A[i+1] = valor+1;
		CSR_COL_IND_A[i+2] = valor+2;
		++valor;
	}

	CSR_COL_IND_A[COL_SIZE-2] = valor;
	CSR_COL_IND_A[COL_SIZE-1] = valor+1;

	// Compute CSR values
	CSR_VAL_A[0] = D[0];
	CSR_VAL_A[1] = U[0];

	int d_index=1; 	
	for(i=2;d_index<SIZE-1;i+=3){
		CSR_VAL_A[i] = L[d_index-1];
		CSR_VAL_A[i+1] = D[d_index];
		CSR_VAL_A[i+2] = U[d_index];
		++d_index;
	}
	
	CSR_VAL_A[COL_SIZE-2] = L[d_index-1];
	CSR_VAL_A[COL_SIZE-1] = D[d_index];		

}

 
