 /**
 *
 * 	@file memm.cu
 *
 * 	@brief Memory Management functions
 *
 * 	@author Ivan Martinez-Perez ivan.martinez@bsc.es
 * 	@author Pedro Valero-Lara   pedro.valero@bsc.es 
 * 	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

// #define COMP_PADDING
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <omp.h>
#include "memm.h"
#include "utils.h"
#include "cuThomasVBatch.h"

using namespace std;

double init;

double *u_device;
double *l_device;
double *d_device;
int *p_device;
double *rhs_device;
int *sys_size_device;

double *u_interleave;
double *l_interleave;
double *d_interleave;
int *p_interleave;
double *rhs_interleave;
double *rhs_interleave_input;
double *rhs_interleave_test;
double *rhs_interleave_device;

int chunk;
int inter;
int normal;
double seq_time;

double *u_host;
double *l_host;
double *d_host;
int *p_host;
// Contains a non-modified version of the d array, used to copy the original data into the GPU
double *d_input;
// Contains a non-modified version of the u array, used to copy the original data into the GPU
double *u_input;
// Contains the value calculated for the last version executed (except sequential)
double *rhs_host;
// Contains a non-modified version of the rhs array
double *rhs_input;
// Contains the result of sequential cpu execution only
double *rhs_host_test;


int data_input(string input_path)
{
    ifstream infile;
    infile.open(input_path.c_str(), ios::in);
    if (!infile.is_open())
    {
        cout << "Error: file not exist: " << input_path;
        std::exit(-1);
    }

    int index;
    int maxsize;
    infile >> maxsize;
    cout << input_path << endl;
    for (int i = 0; i < maxsize; i++)
    {
        infile >> index;
        infile >> u_host[index] >> l_host[index] >> rhs_host[index] >> d_host[index];
        infile >> p_host[index];
        u_input[index] = u_host[index];
        rhs_input[index] = rhs_host[index];
        rhs_host_test[index] = rhs_host[index];
    }
    
    infile.close();
    return N;
}

void data_output(string output_path)
{
    ofstream outfile;
    outfile.open(output_path.c_str(), ios::out);
    if (!outfile.is_open())
    {
        cout << "Error: open output file failed." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
    {
        outfile << i << " " << u_host[i] << " " << l_host[i] << " " << rhs_host[i] << " " << d_host[i] << endl;
    }
    outfile.close();
}


// Sequential implementation

void Sequential( double *u,  double *l, double *d,
                 double *rhs,
                 int ncells, int cell_size ) 
{
	int first,last;
	//printf("BATCHCOUNT = %d, SIZE = %d\n",N,n);
    for (int j = 0; j < ncells; ++j)
    {
        first = j*cell_size;
        last = first + cell_size;
        last--; 

        u[first] /= d[first];
        rhs[first] /= d[first];

        for (int i = first+1; i < last; i++) {
            u[i] /= d[i] - l[i]*u[i-1];
            rhs[i] = (rhs[i] - l[i]*rhs[i-1]) / (d[i] - l[i]*u[i-1]);
        }

        rhs[last] = (rhs[last] - l[last]*rhs[last-1]) / (d[last] - l[last]*u[last-1]);

        for (int i = last-1; i >= first; i--) {
            rhs[i] -= u[i]*rhs[i+1];
        }
    }
}

void Sequential_diff( double *u,  double *l, double *d,
                 double *rhs,
                 int ncells, int maxsize, int *sys_size ) 
{
	int first,last;
	//printf("BATCHCOUNT = %d, SIZE = %d\n",N,n);
    for (int j = 0; j < ncells; ++j)
    {
        first = j*maxsize;
        last = first + sys_size[j];
        last--; 

        u[first] /= d[first];
        rhs[first] /= d[first];

        for (int i = first+1; i < last; i++) {
            u[i] /= d[i] - l[i]*u[i-1];
            rhs[i] = (rhs[i] - l[i]*rhs[i-1]) / (d[i] - l[i]*u[i-1]);
        }

        rhs[last] = (rhs[last] - l[last]*rhs[last-1]) / (d[last] - l[last]*u[last-1]);

        for (int i = last-1; i >= first; i--) {
            rhs[i] -= u[i]*rhs[i+1];
        }
    }


}


// OpenMP (multicore) implementation

void Multicore( double *u, double *l, double *d,
                double *rhs,
                int ncells, int cell_size ) 
{
    #pragma omp parallel for shared(l,d,u,rhs)
    for (int j = 0; j < ncells; ++j)
    {
        int first = j*cell_size;
        int last = first + cell_size;

        last--;
		
		u[first] /= d[first];
        rhs[first] /= d[first];



        for (int i = first+1; i < last; i++) {
            u[i] /= d[i] - l[i]*u[i-1];
            rhs[i] = (rhs[i] - l[i]*rhs[i-1]) / (d[i] - l[i]*u[i-1]);
        }

        rhs[last] = (rhs[last] - l[last]*rhs[last-1]) / (d[last] - l[last]*u[last-1]);

        for (int i = last-1; i >= first; i--) {
            rhs[i] -= u[i]*rhs[i+1];
        }
 
		
    }        
	       
}

void Multicore_diff( double *u,  double *l, double *d,
                 double *rhs,
                 int ncells, int maxsize, int *sys_size ) 
{
    #pragma omp parallel for shared(l,d,u,rhs)
    for (int j = 0; j < ncells; ++j)
    {
        int first = j*maxsize;
        int last = first + sys_size[j];

        last--;
		
		u[first] /= d[first];
        rhs[first] /= d[first];

        for (int i = first+1; i < last; i++) {
            u[i] /= d[i] - l[i]*u[i-1];
            rhs[i] = (rhs[i] - l[i]*rhs[i-1]) / (d[i] - l[i]*u[i-1]);
        }

        rhs[last] = (rhs[last] - l[last]*rhs[last-1]) / (d[last] - l[last]*u[last-1]);

        for (int i = last-1; i >= first; i--) {
            rhs[i] -= u[i]*rhs[i+1];
        }
 
		
    }        
	       
}

// Memory Management functions

void allocateGlobalVar(int maxsize){
    cout << "begin allocate Global Var" << endl;
    u_host = (double*) malloc(N*maxsize*sizeof(double));
    l_host = (double*) malloc(N*maxsize*sizeof(double));
    d_host = (double*) malloc(N*maxsize*sizeof(double));
    p_host = (int*) malloc(N*maxsize*sizeof(int));
    d_input = (double*) malloc(N*maxsize*sizeof(double));
    u_input = (double*) malloc(N*maxsize*sizeof(double));
    rhs_host = (double*) malloc(N*maxsize*sizeof(double));
    rhs_input = (double*) malloc(N*maxsize*sizeof(double));
    rhs_host_test = (double*) malloc(N*maxsize*sizeof(double));

    u_interleave = (double*) malloc(N*maxsize*sizeof(double));
    l_interleave = (double*) malloc(N*maxsize*sizeof(double));
    d_interleave = (double*) malloc(N*maxsize*sizeof(double));
    p_interleave = (int*) malloc(N*maxsize*sizeof(int));
    rhs_interleave = (double*) malloc(N*maxsize*sizeof(double));
    rhs_interleave_input = (double*) malloc(N*maxsize*sizeof(double));
    rhs_interleave_test = (double*) malloc(N*maxsize*sizeof(double));
    cout << "end allocate Global Var" << endl;
}

void allocateGPUVar(int maxsize){

    cudaMalloc(&u_device,N*maxsize*sizeof(double));
    cudaMalloc(&l_device,N*maxsize*sizeof(double));
    cudaMalloc(&d_device,N*maxsize*sizeof(double));
    cudaMalloc(&rhs_device,N*maxsize*sizeof(double));
}

void freeGlobalVar(){

    free(u_host);  
    free(l_host);
    free(d_host);
    free(d_input); 
    free(rhs_host); 
    free(rhs_input); 
    free(rhs_host_test); 

    free(u_interleave); 
    free(l_interleave);
    free(d_interleave); 
    free(rhs_interleave);
    free(rhs_interleave_input); 
    free(rhs_interleave_test);

}

void initRandGlobalVar(int *sys_size, int maxsize){

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < maxsize; ++j)
        {
            if (j<sys_size[i])
            {
                u_host[(i * maxsize) + j] = doubleRand((double)-2.,(double)2.);
                u_input[(i * maxsize) + j] = u_host[(i * maxsize) + j];
                l_host[(i * maxsize) + j] = u_host[(i*maxsize)+j];//doubleRand((double)-2.,(double)2.);
                d_host[(i * maxsize) + j] = doubleRand((double)5.,(double)10.);
                d_input[(i * maxsize) + j] = d_host[(i * maxsize) + j];
                rhs_host[(i * maxsize) + j] = doubleRand((double)-2.,(double)2.);
                rhs_input[(i * maxsize) + j] = rhs_host[(i * maxsize) + j];
                rhs_host_test[(i * maxsize) + j] = rhs_host[(i * maxsize) + j];
            }
            else{
                u_host[(i * maxsize) + j] = 0;
                u_input[(i * maxsize) + j] = 0;
                l_host[(i * maxsize) + j] = 0;
                d_host[(i * maxsize) + j] = 0;
                d_input[(i * maxsize) + j] = 0;
                rhs_host[(i * maxsize) + j] = 0;
                rhs_input[(i * maxsize) + j] = 0;
                rhs_host_test[(i * maxsize) + j] = 0;
            }
        }
    }
}


void initGPUGlobalVar(int maxsize){
    for (int i = 0; i < maxsize; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            u_interleave[i*N+j] = u_input[j*maxsize+i];
            l_interleave[i*N+j] = l_host[j*maxsize+i];
            d_interleave[i*N+j] = d_input[j*maxsize+i];
            rhs_interleave[i*N+j] = rhs_input[j*maxsize+i];
            rhs_interleave_test[i*N+j] = rhs_host_test[j*maxsize+i];
            rhs_interleave_input[i*N+j] = rhs_input[j*maxsize+i];
        }
    }
}

void initInterleavedGlobalVar(int maxsize){

    for (int i = 0; i < maxsize; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            u_interleave[i*N+j] = u_input[j*maxsize+i];
            l_interleave[i*N+j] = l_host[j*maxsize+i];
            d_interleave[i*N+j] = d_input[j*maxsize+i];
            rhs_interleave[i*N+j] = rhs_input[j*maxsize+i];
            rhs_interleave_test[i*N+j] = rhs_host_test[j*maxsize+i];
            rhs_interleave_input[i*N+j] = rhs_input[j*maxsize+i];
        }
    }

}

// Wrapers 

void execSeqHines(int maxsize,int *sys_size){

    init = time_wtime();
    #ifdef COMP_PADDING
        //printf("COMPUTING PADDING\n"); 
        Sequential(
            u_host,
            l_host,
            d_host,
            rhs_host_test,
            N,
            maxsize
        );
    #else
        //printf("NOT COMPUTING PADDING\n"); 
        Sequential_diff(
            u_host,
            l_host,
            d_host,
            rhs_host_test,
            N,
            maxsize,
            sys_size

        );
    #endif
    seq_time = time_wtime()-init;
    printf("        CPU SEQ %e\n", seq_time);

}

void execMultiHines(int maxsize, int *sys_size){
    
	// Reset modified pointers
    for (int i = 0; i < N; ++i)
	{
        for (int j = 0; j < maxsize; ++j)
        {
            d_host[(i * maxsize) + j] = d_input[(i * maxsize) + j];
            u_host[(i * maxsize) + j] = u_input[(i * maxsize) + j];
            rhs_host[(i * maxsize) + j] = rhs_input[(i * maxsize) + j];
		}
	}
        
    init = time_wtime();
    #ifdef COMP_PADDING
        Multicore(
            u_host,
            l_host,
            d_host,
            rhs_host,
            N,
            maxsize
        );
    #else
        Multicore_diff(
            u_host,
            l_host,
            d_host,
            rhs_host,
            N,
            maxsize,
            sys_size
		);
	#endif

    double  multi_time = time_wtime()-init;
    //double multi_speedup = seq_time/multi_time;
    printf("        CPU MUL %e          ",multi_time );
    //printf("%1.2f\n",multi_speedup);
    calcError(rhs_host,rhs_host_test,N*maxsize);

}

void execGpuHines(int maxsize, int *sys_size){

    cudaMemcpy(u_device,u_interleave,N*maxsize*sizeof(double),cudaMemcpyHostToDevice);
    check();

    cudaMemcpy(l_device,l_interleave,N*maxsize*sizeof(double),cudaMemcpyHostToDevice);
    check();

    cudaMemcpy(d_device,d_interleave,N*maxsize*sizeof(double),cudaMemcpyHostToDevice);
    check();

    cudaMalloc(&rhs_interleave_device,N*maxsize*sizeof(double));
    cudaMemcpy(rhs_interleave_device,rhs_interleave_input,N*maxsize*sizeof(double),cudaMemcpyHostToDevice);
    check();

    check();

    cudaMalloc(&sys_size_device,N*sizeof(int));
    cudaMemcpy(sys_size_device,sys_size,N*sizeof(int),cudaMemcpyHostToDevice);
    check();

    init = time_wtime();
        
		cuThomasVBatch<<<N/BlockSize, BlockSize>>>
        (
            u_device, l_device, d_device, rhs_interleave_device,
            N, maxsize
        );
        check();
        
    cudaDeviceSynchronize();
    double gpu_time = time_wtime()-init;
    //double gpu_speedup = seq_time/gpu_time; 
    printf("        GPU %e            ", gpu_time);
    //printf("%1.2f\n", gpu_speedup);
    //printf("%e\n", gpu_time);
    cudaMemcpy(rhs_host,rhs_interleave_device,N*maxsize*sizeof(double),cudaMemcpyDeviceToHost);
    calcError(rhs_host,rhs_interleave_test,N*maxsize);

    cudaDeviceReset();

}
