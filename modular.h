#ifndef MODULAR_H
#define MODULAR_H

#include "def.h"

/* Compute Rmod <- A + B */ 
/* Input: 0 <= A, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Add_mod(biguint_t Rmod, bigint_t cy, const biguint_t A, const biguint_t B);

/* Compute Rmod <- Rmod + B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 3*N, 0 < B < 7*N ) */ 
/* Ouput: 0 <= Rmod < 6*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 10*N) */ 
__device__ void Cuda_Add_mod(biguint_t Rmod, bigint_t cy, const biguint_t A);

/* Compute Rmod <- Rmod - B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Sub_mod(biguint_t Rmod, bigint_t cy, const biguint_t B, const digit_t N3thdx);

/* Compute r <- 2*a */ 
/* Input: 0 <= a < 3*N */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Dbl_mod(biguint_t r, biguint_t a);

/* Compute r <- A*b */ 
/* Input: 0 < b < 2^SIZE_DIGIT, 0 <= A < 6*N */ 
/* Ouput: 0 <= r < 7*N */ 
__device__ void Cuda_Mulint_mod(biguint_t r, bigint_t cy, biguint_t A, digit_t b, const digit_t Nthdx,const digit_t invN);

/* Compute r <- A*B */ 
/* Input: 0 <= A, B < 6*N */
/* (except when it follows Cuda_Mulint_mod, 0 <= A < 6*N, 0 < B < 10*N ) */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Mul_mod(biguint_t mul, bigint_t cy, const biguint_t A, const biguint_t B, biguint_t r,const digit_t Nthdx, const digit_t invN);

/* Compute r <- A*A */ 
/* Input: 0 <= A < 6*N */
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Square_mod(biguint_t mul, bigint_t cy, const biguint_t A, biguint_t r, const digit_t Nthdx, const digit_t invN);

#endif
