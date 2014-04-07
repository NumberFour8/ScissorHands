#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"
#include "edwards.h"

#define NUM_CURVES 1

// Paralelní metody pro zdvojnásovení a součet bodů v Extended souřadnicích, a = -1
__global__ void edwardsAdd(ExtendedPoint* R,ExtendedPoint *P,ExtendedPoint *Q);
__global__ void edwardsDbl(ExtendedPoint* R,ExtendedPoint *P);

/*
  Non-adjacent Form
 */
struct NAF {
	char* bits;
	unsigned int l;
	unsigned char w;
};

/*
 Pomocná struktura pro N, 3*N a inverzi N modulo velikost báze
*/
struct h_Aux {
	biguint_t N,N3;
	digit_t invN;
};

// Hlavní výpočetní metoda pro Extended souřadnice
extern "C" cudaError_t computeExtended(const h_Aux input,h_ExtendedPoint* initPoints,const NAF coeff);

#endif
