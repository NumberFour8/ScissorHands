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
	#include <gmp.h>
	#define ASM asm __volatile__
#else
	#include <windows.h>
	#include <mpir.h>
	#define ASM asm volatile
#endif

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <device_launch_parameters.h>

#define NB_DIGITS 32 
#define SIZE_DIGIT 32

typedef unsigned int digit_t;
typedef int carry_t;

#define VOL volatile

#define MAX_BITS SIZE_DIGIT*NB_DIGITS
#define MAX_BYTES MAX_BITS/8

typedef digit_t VOL biguint_t[NB_DIGITS];
typedef carry_t VOL bigint_t[NB_DIGITS];

#endif
