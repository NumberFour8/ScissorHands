#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>

#include "modular.h"
#include "coords.h"
#include "helpers.h"

#define NUM_BLOCKS 5
#define CURVES_PER_BLOCK 8

#define NUM_CURVES NUM_BLOCKS*CURVES_PER_BLOCK

// Hlavní výpočetní metoda 
cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff);

#endif
