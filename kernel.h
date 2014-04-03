#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"

__global__ void myKernel(biguint_t a, biguint_t b, biguint_t n,biguint_t c);
extern "C" cudaError_t testCuda(biguint_t a,biguint_t b,biguint_t c,biguint_t n);

#endif
