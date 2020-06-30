/**
 *
 * 	@file parallel.cc
 *
 * 	@brief main function
 *
 * 	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <string>
#include "memm.h"

using namespace std;


/* Global Variables */
int N;
int BlockSize;

#define NUM_SIZE  6
#define NUM_TESTS 12

//Value of the minimun size of a system
int Size[NUM_SIZE] = {1,2,4,8,16,32};
//Number of systems
int Nsys[NUM_TESTS] = {76, 76, 305, 322, 698, 699, 74277, 181851, 568449, 921601, 2752001, 5728200};


string input_base_path = "../data/case";
string output_base_path = "../sresult/res";

int main(int argc, char const *argv[])
{
    string input_path;
    string output_path;

    if (argc == 3) {
        input_path = argv[1];
        output_path = argv[2];

        int maxSize = N;
        allocateGlobalVar(maxSize);

		// Input Data
		data_input(input_path);

        // CPU serial
        execSeqHines(maxSize, NULL);

        // cuHinesBatch
        initGPUGlobalVar(maxSize);

        // initInterleavedGlobalVar(maxSize);
        execGpuHines(maxSize,NULL);


        // Free
        freeGlobalVar();

		// Output Data
    	data_output(output_path);

    } else if (argc == 1) {
	BlockSize  = 128;

    for (int i = 1 ; i <= NUM_TESTS; ++i)
    {
		N=Nsys[i-1];
        int *sys_size = (int*) malloc(N*sizeof(int));
        int tmpSize;
        int maxSize;
        printf("N ---- %d\n",N);
       
        for (int minsize = 0; minsize <NUM_SIZE; ++minsize)
        {
            maxSize=0;
            printf("    RANGE --- %d-%d - ",Size[minsize], 32);
			std::default_random_engine generator;
			std::uniform_int_distribution<int> dist(Size[minsize], 32);
            for (int i = 0; i < N; ++i)
            {
                tmpSize= dist(generator);
                sys_size[i] = tmpSize;
                if(maxSize<tmpSize) maxSize=tmpSize;
            }
			
			int mean=0;
			for(int i=0;i<N;++i){
				mean+=sys_size[i];
			}
			printf("MEAN:%.2f\n",(double)mean/N);			

            allocateGlobalVar(maxSize);
            allocateGPUVar(maxSize);

            // Input Data
		    data_input(input_base_path + to_string(i) + ".txt");

            // CPU serial
            execSeqHines(maxSize,sys_size);

            // CPU multicore 
            execMultiHines(maxSize,sys_size);
            
            // cuHinesBatch
            initGPUGlobalVar(maxSize);

           // initInterleavedGlobalVar(maxSize);
            execGpuHines(maxSize,sys_size);

            // Free
            freeGlobalVar();

            // Output Data
    	    data_output(output_base_path + to_string(i) + ".txt");
	        printf("\n");	
            }
        }
    } else {
        cout << "Invalid Argument" << endl;
        exit(-1);
    }

	return 0;
}
