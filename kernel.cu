#include <iostream>
#include <stdio.h>
using namespace std;

#include "modular.h"

__global__ void myKernel(biguint_t a, biguint_t b, biguint_t n,biguint_t c)
{
    __shared__ VOL carry_t b_cy[NB_DIGITS]; 
	__shared__ VOL digit_t r[NB_DIGITS];

	memset((void*)r,0,NB_DIGITS*sizeof(int));
	memset((void*)b_cy,0,NB_DIGITS*sizeof(int));

	Cuda_Mul_mod(c,b_cy,a,b,r,n[threadIdx.x],2047647423);
	Cuda_Add_mod(c,b_cy,b);
}

extern "C" cudaError_t testCuda(biguint_t a,biguint_t b,biguint_t c,biguint_t n)
{
	void* devA = NULL,*devB = NULL,*devC = NULL,*devN = NULL;
    cudaError_t cudaStatus;
    
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return cudaStatus;
    }

    cudaMalloc((void**)&devA, MAX_BYTES);
    cudaMalloc((void**)&devB, MAX_BYTES);
	cudaMalloc((void**)&devC, MAX_BYTES);
	cudaMalloc((void**)&devN, MAX_BYTES);
    
    // Copy input vectors from host memory to GPU buffers.
	cudaMemcpy(devA, (void*)a, MAX_BYTES, cudaMemcpyHostToDevice);
    cudaMemcpy(devB, (void*)b, MAX_BYTES, cudaMemcpyHostToDevice);
    cudaMemcpy(devC, (void*)c, MAX_BYTES, cudaMemcpyHostToDevice);
    cudaMemcpy(devN, (void*)n, MAX_BYTES, cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
	myKernel<<<1, NB_DIGITS>>>((unsigned int*)devA, (unsigned int*)devB, (unsigned int*)devN,(unsigned int*)devC);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess)
      fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess)
	  fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));

    // Copy output vector from GPU buffer to host memory.
	cudaMemcpy((void*)c, devC, MAX_BYTES, cudaMemcpyDeviceToHost);
    
    cudaFree(devC);
    cudaFree(devA);
    cudaFree(devB);
    
    return cudaStatus;
}

//////////////////////////////////////////////////////////////////////////////////////

void printBigInt(biguint_t B)
{
	for (int i = 0;i < NB_DIGITS;++i)
	{
		printf("%#010x",B[i]);
		if (B[i+1] == 0) break;
		if (i != NB_DIGITS-1) printf(",");
	}
	printf("\n");
}

int main()
{
	// A,B jsou v Montgomeryho reprezentaci, A,B,N v bázi 2^32
	// N = 215714093118538583256769
	// A = 21799067859837164737
	// B = 104402829964868711809

	biguint_t A = {0x0e6c2ef9,0x9e4d4f27,0x0000299e,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000};
	biguint_t B = {0x47769205,0xdceddf18,0x00002c89,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000};
	biguint_t N = {0x1b8a2ec1,0xe2695510,0x00002dad,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000};
	biguint_t C = {0};

    cudaError_t cudaStatus = testCuda(A,B,C,N);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
    
	printBigInt(C);
	
	char c;
	cin >> c;
    return 0;
}
