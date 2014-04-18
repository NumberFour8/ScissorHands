#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"
#include "coords.h"
#include "helpers.h"

#define CURVE_MEMORY_SIZE 1536
#define CURVES_PER_BLOCK 16

// Hlavní výpočetní metoda 
cudaError_t compute(const ComputeConfig h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff);

//////////////////////////////////////////// INTERNÍ MAKRA KERNELŮ ////////////////////////////////////////////

// Inicializace kernelu
#define PREPARE() ComputeConfig *ax = (ComputeConfig*)swAux;			\
	__shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];			\
	VOL digit_t* c_x1 = x1[threadIdx.y];								\
	VOL digit_t* c_y1 = y1[threadIdx.y];								\
	VOL digit_t* c_z1 = z1[threadIdx.y];								\
	VOL digit_t* c_t1 = t1[threadIdx.y];								\
	VOL digit_t* c_x2 = x2[threadIdx.y];								\
	VOL digit_t* c_y2 = y2[threadIdx.y];								\
	VOL digit_t* c_z2 = z2[threadIdx.y];								\
	VOL digit_t* c_t2 = t2[threadIdx.y];								\
	VOL digit_t* c_tt0  = temp0[threadIdx.y];							\
	VOL digit_t* c_tt1  = temp1[threadIdx.y];							\
	VOL carry_t* _CARRY = carry[threadIdx.y];							\
	VOL digit_t* _AUX   = temp2[threadIdx.y];							\
	const digit_t _N    = ax->N[threadIdx.x];							\
	const digit_t _3N   = ax->N3[threadIdx.x];							\
	const digit_t _INVN = ax->invN;										\
	const digit_t idx = 4*NB_DIGITS*(blockIdx.x*blockDim.y + threadIdx.y);

// Vyčištění dočasných proměnných
#define CLEAR_TEMP() c_tcy[threadIdx.x] = 0;	\
					 c_tt0[threadIdx.x] = 0;	\
					 c_tt1[threadIdx.x] = 0;	\
					 _AUX[threadIdx.x]  = 0; 
 
// Nakopírování pracovních dat bodu Q (celkem 32*4 = 128 cifer) 
#define COPY_Q()	c_x2[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS);  \
					c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS);  \
					c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS);  \
					c_t2[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); 


#endif
