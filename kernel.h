#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"
#include "coords.h"
#include "helpers.h"

#define NUM_BLOCKS 10
#define CURVES_PER_BLOCK 16

#define NUM_CURVES NUM_BLOCKS*CURVES_PER_BLOCK


// Hlavní výpočetní metoda 
cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff,const unsigned int WS);

#endif
