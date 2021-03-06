#ifndef DEF_H
#define DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#ifndef _MSC_VER
	#include <unistd.h>
	#include <sys/resource.h>
	#define ASM asm __volatile__
#else
	#include <windows.h>
	#define ASM asm volatile
#endif

#include <gmp.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <device_launch_parameters.h>

#define VOL volatile

// Konfigurace kernel�
#define CURVE_MEMORY_SIZE 1536
#define CURVES_PER_BLOCK 16

// Konfigurace aritmetiky
#define NB_DIGITS 32 
#define SIZE_DIGIT 32

typedef unsigned int digit_t;
typedef int carry_t;

#define MAX_BITS SIZE_DIGIT*NB_DIGITS
#define MAX_BYTES MAX_BITS/8

typedef digit_t VOL biguint_t[NB_DIGITS];
typedef carry_t VOL bigint_t[NB_DIGITS];

#endif
