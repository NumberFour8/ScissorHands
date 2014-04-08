#ifndef HELPERS_H
#define HELPERS_H

#include "def.h"

/*
   Struktura obsahující výsledek pøevodu z Extended souøadnic
*/
struct ExtResult {
	bool factorFound;
	mpz_t factor;
};

/*
  Non-adjacent form
 */
class NAF {
public:
	char* bits;
	unsigned int l;
	const unsigned char w;
	NAF(unsigned char w,mpz_t number);
	virtual ~NAF();
	
	int build(unsigned int start,unsigned int end) const;
};

void to_mont_repr (mpz_t x, mpz_t n);
void from_mont_repr (mpz_t x, mpz_t n,mpz_t invB);

void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);

void printmpz(const char* format,mpz_t number);
void printBigUInt(biguint_t B);


// Pomocná makra na zachycení chyb v CUDA
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define cuda_Malloc(A,B) gpuErrchk(cudaMalloc(A,B))
#define cuda_Memcpy(A,B,C,D) gpuErrchk(A,B,C,D)
#define cuda_Free(A) gpuErrchk(A)
#define cuda_MemcpyToSymbol(A,B,C) gpuErrchk(cudaMemcpyToSymbol(A,B,C))

#endif
