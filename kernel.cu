
#include "def.h"

#include <iostream>
#include <stdio.h>
using namespace std;


extern "C" cudaError_t testCuda(biguint_t a,biguint_t b,biguint_t c,biguint_t n);

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
