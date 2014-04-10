#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"
#include "coords.h"
#include "helpers.h"

#define NUM_BLOCKS 1
#define CURVES_PER_BLOCK 1


/*
 Pomocná struktura pro N, 3*N a inverzi N modulo velikost báze
*/
struct h_Aux {
	biguint_t N,N3;
	digit_t invN;
};

// Hlavní výpočetní metoda pro Extended souřadnice
cudaError_t computeExtended(const h_Aux input,h_ExtendedPoint* initPoints,const NAF coeff);

#endif
