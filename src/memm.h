/**
 *
 * 	@file memm.h
 *
 * 	@brief header of memm.cu  
 *
 * 	@author Ivan Martinez-Perez ivan.martinez@bsc.es
 * 	@author Pedro Valero-Lara   pedro.valero@bsc.es 
 * 	@author Yifan Zhang   yifzhang@pku.edu.cn
 *
 **/

#include <string>

extern int N;
extern int BlockSize;

/* IO */
int data_input(std::string input_path);
void data_output(std::string output_path);

/* MEMM */
void allocateGlobalVar(int maxsize);
void allocateGPUVar(int maxsize);
void initRandGlobalVar(int *sys_size, int maxsize);
void initGPUGlobalVar(int maxsize);
void initInterleavedGlobalVar(int maxsize);
void freeGlobalVar();

/* Wraper */
void execGpuHines(int maxsize, int *sys_size);
void execMultiHines(int maxsizeint, int *sys_size);
void execSeqHines(int maxsizeint, int *sys_size);
 
