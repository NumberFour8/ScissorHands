#ifndef MODULAR_H
#define MODULAR_H

#include "def.h"

/* Compute Rmod <- A + B  
   Input: 0 <= A, B < 3*N  
   Ouput: 0 <= Rmod < 6*N 
*/ 
__device__ void Cuda_Add_mod(biguint_t Rmod, bigint_t cy, const biguint_t A, const biguint_t B);

/* Compute Rmod <- Rmod + B 
   Input: 0 <= Rmod, B < 3*N 
   (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 3*N, 0 < B < 7*N ) 
   Ouput: 0 <= Rmod < 6*N 
  (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 10*N)
*/ 
__device__ void Cuda_Add_mod(biguint_t Rmod, bigint_t cy, const biguint_t A);

/* Compute Rmod <- Rmod - B  
   Input: 0 <= Rmod, B < 3*N  
   Ouput: 0 <= Rmod < 6*N 
*/ 
__device__ void Cuda_Sub_mod(biguint_t Rmod, bigint_t cy, const biguint_t B, const digit_t N3thdx);

/* Compute r <- 2*a  
   Input: 0 <= a < 3*N  
   Ouput: 0 <= r < 3*N 
*/ 
__device__ void Cuda_Dbl_mod(biguint_t r, biguint_t a);

/* Compute r <- A*b 
   Input: 0 < b < 2^SIZE_DIGIT, 0 <= A < 6*N  
   Ouput: 0 <= r < 7*N 
 
__device__ void Cuda_Mulint_mod(biguint_t r, bigint_t cy, biguint_t A, digit_t b, const digit_t Nthdx,const digit_t invN);
*/

/* Compute r <- A*B 
   Input: 0 <= A, B < 6*N 
   (except when it follows Cuda_Mulint_mod, 0 <= A < 6*N, 0 < B < 10*N )  
   Ouput: 0 <= r < 3*N 
*/ 
__device__ void Cuda_Mul_mod(biguint_t mul, bigint_t cy, const biguint_t A, const biguint_t B, biguint_t r,const digit_t Nthdx, const digit_t invN);

/* Compute r <- A*A  
   Input: 0 <= A < 6*N 
   Ouput: 0 <= r < 3*N 
*/ 
__device__ void Cuda_Square_mod(biguint_t mul, bigint_t cy, const biguint_t A, biguint_t r, const digit_t Nthdx, const digit_t invN);

// Názvy pomocných proměnných
#define _CARRY c_tcy
#define _AUX   c_aux
#define _N	   c_N
#define _3N	   c_3N
#define _INVN  c_invN

// Zjednodušující makra pro modulární aritmetiku

// A = B
#define EQL_MOD(A,B)   A[threadIdx.x] = B[threadIdx.x]

// C = A+B
#define ADD_MOD(C,A,B) Cuda_Add_mod(C,_CARRY,A,B)

// C = A*B
#define MUL_MOD(C,A,B) Cuda_Mul_mod(C,_CARRY,A,B,_AUX,_N,_INVN) 

// C = A-B
#define SUB_MOD(C,A,B) EQL_MOD(C,A); \
					   Cuda_Sub_mod(C,_CARRY,B,_3N)

// C = C-B
#define SUE_MOD(C,B)   Cuda_Sub_mod(C,_CARRY,B,_3N)

// A = 2*A
#define DBL_MOD(A)	   Cuda_Dbl_mod(A,A)

// A = -B
#define NEG_MOD(A,B)   SUB_MOD(A,_N,B)

// A = B^2
#define SQR_MOD(A,B)   Cuda_Square_mod(B,_CARRY,A,_AUX,_N,_INVN)

#endif
