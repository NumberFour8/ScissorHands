#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#ifndef _MSC_VER
#include <unistd.h>
#include <sys/resource.h>
#else
#include <windows.h>
#endif

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define VOL volatile 

#define NB_DIGITS 32 
#define SIZE_DIGIT 32

typedef unsigned int digit_t;
typedef int carry_t;

#define MAX_BITS SIZE_DIGIT * NB_DIGITS

typedef digit_t VOL biguint_t[NB_DIGITS];
typedef carry_t VOL bigint_t[NB_DIGITS];

#define MAX_BYTES MAX_BITS/8

#if defined(_MSC_VER)
#  define ASM asm volatile
#else
#  define ASM asm __volatile__
#endif