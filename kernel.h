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
struct Aux {
	biguint_t N;
	biguint_t N3;
	digit_t invN;
};

// Hlavní výpočetní metoda 
cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF coeff);

#endif
