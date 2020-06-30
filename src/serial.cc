/**
 *
 * 	@file serial.cc
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

int main(int argc, char *argv[])
{
    N = 1;
    for (int i = 1; i <= 12; i++) {
		cout << "#Case: " << to_string(i) << endl;
        
        int maxSize = Nsys[i-1];
        allocateGlobalVar(maxSize);

		// Input Data
		data_input(input_base_path + to_string(i) + ".txt");

        // CPU serial
        execSeqHines(maxSize, NULL);

        // Free
        freeGlobalVar();

		// Output Data
    	data_output(output_base_path + to_string(i) + ".txt");
    }

	return 0;
}
