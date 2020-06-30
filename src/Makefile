 #
 #
 # @file Makefile
 #
 #

CC		   = g++
NVCC       = nvcc
LD_LIBRARY_PATH= /usr/local/cuda/lib64
CUDA_LIB   = -L /usr/local/cuda/lib64/ -lcuda -lcudart -lcusparse -lcusolver 
ARCH_CUDA  = -arch=sm_30
NVCC_FLAGS = --ptxas-options=-v -Xcompiler -fopenmp -O3 -std=c++11  -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES  
GCC_FLAGS  = -fopenmp -std=c++11 -Wall -pedantic  


all: clean serial parallel remove

serial: memm.o cuThomasVBatch.o serial.cc
	$(CC)   $(GCC_FLAGS) -o serial serial.cc $(CUDA_LIB) memm.o cuThomasVBatch.o	

parallel: memm.o cuThomasVBatch.o parallel.cc
	$(CC)   $(GCC_FLAGS) -o parallel parallel.cc $(CUDA_LIB) memm.o cuThomasVBatch.o	

memm.o: memm.cu

	$(NVCC) -c  $(NVCC_FLAGS) $(ARCH_CUDA) memm.cu

cuThomasVBatch.o: cuThomasVBatch.cu 
	$(NVCC) -c  $(NVCC_FLAGS) $(ARCH_CUDA) cuThomasVBatch.cu	 

clean:
	rm -rf *.o serial parallel

remove:
	rm -rf *.o
