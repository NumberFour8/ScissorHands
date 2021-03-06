#include "modular.h"

// Makra instrukc� v sad� PTX
#define __add_cc(r,a,b) 	ASM ("add.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b)) 
#define __addc_cc(r,a,b) 	ASM ("addc.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b))
#define __sub_cc(r,a,b) 	ASM ("sub.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b)) 

#define __addcy(carry) 		ASM ("addc.s32 %0, 0, 0;": "=r"(carry)) 
#define __addcy2(carry) 	ASM ("addc.cc.s32 %0, %0, 0;": "+r"(carry)) 

#define __subcy(carry) 		ASM ("subc.s32 %0, 0, 0;": "=r"(carry)) 
#define __subcy2(carry) 	ASM ("subc.s32 %0, %0, 0;": "+r"(carry)) 

#define __mul_lo(r,a,b) 	ASM("mul.lo.u32 %0, %1, %2;": "=r"(r): "r"(a),"r"(b)) 
#define __mul_hi(r,a,b) 	ASM("mul.hi.u32 %0, %1, %2;": "=r"(r): "r"(a),"r"(b)) 
#define __mad_lo_cc(r,a,b)  ASM("mad.lo.cc.u32 %0, %1, %2, %0;":\
                                                      "+r"(r): "r"(a),"r"(b)) 
#define __madc_hi_cc(r,a,b) ASM("madc.hi.cc.u32 %0, %1, %2, %0;":\
                                                  "+r"(r):"r"(a),"r"(b)) 
                                                  
__device__ void Cuda_Fully_Normalize (biguint_t A, bigint_t cy)
{
  carry_t cytemp;
  unsigned int thm1;

  while(__any(cy[threadIdx.x])!=0)
  {
    thm1 = (threadIdx.x - 1) % NB_DIGITS;
    cytemp = cy[thm1];

    __add_cc(A[threadIdx.x], A[threadIdx.x], cytemp);
  
    if (cytemp >= 0)
      __addcy(cy[threadIdx.x]);
    else /* if (cytemp < 0) */
      __subcy(cy[threadIdx.x]);
  }
}

/* Spo�ti: Rmod <- A + B */ 
/* Vstup: 0 <= A, B < 3*N */ 
/* V�stup: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, bigint_t cy, const biguint_t A, const biguint_t B)
{
  unsigned int thp1 = (threadIdx.x + 1) % NB_DIGITS;
  __add_cc (Rmod[threadIdx.x], A[threadIdx.x], B[threadIdx.x]);
  __addcy2(Rmod[thp1]); 
  __addcy (cy[thp1]);
  Cuda_Fully_Normalize (Rmod, cy); 
}

/* Spo�ti: Rmod <- Rmod + B */ 
/* Vstup: 0 <= Rmod, B < 3*N */ 
/* (vyjma p�i spln�n� p�epoklad� Cuda_Mulint_mod, 0 <= Rmod < 3*N, 0 < B < 7*N ) */ 
/* V�stup: 0 <= Rmod < 6*N */ 
/* (vyjma p�i spln�n� p�edpoklad� Cuda_Mulint_mod, 0 <= Rmod < 10*N) */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, bigint_t cy, const biguint_t A)
{
  unsigned int thp1 = (threadIdx.x + 1) % NB_DIGITS;
  __add_cc (Rmod[threadIdx.x], Rmod[threadIdx.x], A[threadIdx.x]);
  //__addcy (cy[threadIdx.x]);
  __addcy2(Rmod[thp1]); 
  __addcy (cy[thp1]);
  Cuda_Fully_Normalize (Rmod, cy);
}

/* Spo�ti: Rmod <- Rmod - B */ 
/* Vstup: 0 <= Rmod, B < 3*N */ 
/* V�stup: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, bigint_t cy, const biguint_t B, const digit_t N3thdx)
{
  digit_t reg_Rmod = Rmod[threadIdx.x];
  carry_t reg_cy = 0; 
  
  __add_cc (reg_Rmod, reg_Rmod, N3thdx);
  __addcy (reg_cy);
  __sub_cc (reg_Rmod, reg_Rmod, B[threadIdx.x]);
  __subcy2 (reg_cy);

  Rmod[threadIdx.x] = reg_Rmod;
  cy[threadIdx.x] = reg_cy;
  Cuda_Fully_Normalize (Rmod, cy); 
}

// Jeden krok REDC 
__device__ void Cuda_Mulmod_step
(biguint_t r, bigint_t cy, digit_t a, digit_t b, const digit_t Nthdx,
 const digit_t invN)
{
  digit_t t;
  digit_t reg_hi = 0;
  unsigned int thp1= (threadIdx.x + 1) % NB_DIGITS;
  carry_t reg_cy = cy[thp1];

  __mad_lo_cc(r[threadIdx.x],a,b);
  __madc_hi_cc(reg_hi,a,b);
  __addcy2(reg_cy);

  __mul_lo(t, invN, r[0]);
  __mad_lo_cc(r[threadIdx.x],t,Nthdx);
  __madc_hi_cc(reg_hi,t,Nthdx);
  __addcy2(reg_cy);

  // Normalizace a shift doprava
  __add_cc(r[threadIdx.x],r[thp1],reg_hi);
  __addc_cc(r[thp1],r[thp1],reg_cy);
  __addcy(cy[thp1]); 
}

/* Spo�ti: r <- 2*a */ 
/* Vstup: 0 <= a < 3*N */ 
/* V�stup: 0 <= r < 3*N */ 
__device__ void Cuda_Dbl_mod
(biguint_t r, biguint_t a)
{
  unsigned int thp1= (threadIdx.x + 1) % NB_DIGITS;
  asm ("add.cc.u32 %0, %1, %1;" : "=r"(r[threadIdx.x]) : "r"(a[threadIdx.x]));
  __addcy2(r[thp1]);
}


// NEFUNGUJE!
/* Spo�ti: r <- A*b */ 
/* Vstup: 0 < b < 2^SIZE_DIGIT, 0 <= A < 6*N */ 
/* V�stup: 0 <= r < 7*N */ 
/*__device__ void Cuda_Mulint_mod
(biguint_t r, bigint_t cy, biguint_t A, digit_t b, const digit_t Nthdx,const digit_t invN)
{
  digit_t t;
  digit_t reg_hi;
  unsigned int thp1= (threadIdx.x + 1) % NB_DIGITS;
  digit_t reg_A = A[threadIdx.x];
  carry_t reg_cy;

  __mul_lo(r[threadIdx.x],reg_A,b);
  __mul_hi(reg_hi,reg_A,b);

  __mul_lo(t, invN, r[0]);
  __mad_lo_cc(r[threadIdx.x],t,Nthdx);
  __madc_hi_cc(reg_hi,t,Nthdx);
  __addcy(reg_cy);

  __add_cc(r[threadIdx.x],r[thp1],reg_hi);
  __addc_cc(r[thp1],r[thp1],reg_cy);
  __addcy(cy[thp1]); 

  Cuda_Fully_Normalize(r,cy); 
}
*/

/* Spo�ti: r <- A*B */ 
/* Vstup: 0 <= A, B < 6*N */
/* (vyjma p�i spln�n� Cuda_Mulint_mod, 0 <= A < 6*N, 0 < B < 10*N ) */ 
/* V�stup: 0 <= r < 3*N */ 
__device__ void Cuda_Mul_mod 
(biguint_t mul, bigint_t cy, const biguint_t A, const biguint_t B, biguint_t r,
 const digit_t Nthdx, const digit_t invN)
{

  int i;
  digit_t temp=A[threadIdx.x];

  r[threadIdx.x]=0;
  
  for (i=0; i< NB_DIGITS; i++)
    Cuda_Mulmod_step (r, cy, temp, B[i], Nthdx, invN);

  
  Cuda_Fully_Normalize (r, cy);
  mul[threadIdx.x]=r[threadIdx.x];
}

__device__ void Cuda_Square_mod 
(biguint_t mul, bigint_t cy, const biguint_t A, biguint_t r, 
 const digit_t Nthdx, const digit_t invN)
{
  Cuda_Mul_mod (mul, cy, A, A, r, Nthdx, invN);
}

/////////////////////////////////////////////////////////////
