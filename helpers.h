#ifndef HELPERS_H
#define HELPERS_H

#include "def.h"

#define mpz_intz(...) mpz_inits(__VA_ARGS__,NULL)
#define mpz_clrs(...) mpz_clears(__VA_ARGS__,NULL)

// Třída obsahující výsledek převodu z Extended souøadnic
class ExtResult {
public:
	// Indikátor, zda byl nalezen faktor
	bool factorFound;
	
	// Nalezený faktor
	mpz_t factor;

	// Výchozí konstruktor
	ExtResult() : factorFound(false)
	{ mpz_init_set_ui(factor,0); }
	
	~ExtResult()
	{ mpz_clear(factor); }
};

// Rozvoj čísla do non-adjacent form (NAF)
class NAF {
public:
	// Koeficienty rozvoje
	char* bits;
	
	// Délka rozvoje
	unsigned int l;
	
	// Šířka rozvoje
	const unsigned char w;
	
	// Zkontruuje rozvoj šířky W z čísla N
	NAF(unsigned char W,mpz_t N);
	
	virtual ~NAF();
	
	// Vypíše rozvoj
	void print() const;
	
	// Vyčíslí výsek z rozvoje daný parametry start a end
	int build(unsigned int start,unsigned int end) const;
};

// Převede X do Montgomeryho reprezentace modulo N
void to_mont_repr(mpz_t x, mpz_t n);

// Převede X z Montgomeryho reprezentace modulo N
void from_mont_repr(mpz_t x, mpz_t n,mpz_t invB);

// Převede MPZ číslo do báze 2^32
void mpz_to_biguint(biguint_t a, mpz_t b);

// Převede číslo z báze 2^32 do MPZ
void biguint_to_mpz(mpz_t a, biguint_t b);

// Vypíše číslo MPZ
void printmpz(const char* format,mpz_t number);

// Vypíše číslo v bázi 2^32
void printBigInt(const char* tag,biguint_t B);

// Spočítá LCM(1,2,3...,n)
void lcmToN(mpz_t res,const unsigned int n);

// Vynuluje číslo v bázi 2^32
inline void reset(biguint_t n)
{
	memset((void*)n,0,MAX_BYTES);
}

// Pomocná třída pro N, 3*N a inverzi N modulo velikost báze
class Aux {
public:
	biguint_t N;
	biguint_t N3;
	digit_t invN;
	
	Aux(mpz_t N)
	{
		reset(N);
		reset(N3);
		invN = 0;
	
		mpz_t z3N,zW,zInvN;
		mpz_intz(z3,z2,zInvN);
	
		mpz_mul_ui(z3N,zN,3);
	
		mpz_ui_pow_ui (zW, 2, SIZE_DIGIT); 
    
		mpz_invert(zInvN, zN, zW);
		mpz_sub(zInvN, zW, zInvN);

		mpz_to_biguint(N,zN);
		mpz_to_biguint(N3,z3N);
		invN = (digit_t)mpz_get_ui(zInvN);
		
		mpz_clrs(z3N,zW,zInvN);
	}
};

// Vypíše chybu a přeruší program
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//////////////////////// POMOCNÁ MAKRA NA ZACHYCENÍ CHYB V CUDA ////////////////////////

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define cuda_Malloc(A,B) gpuErrchk(cudaMalloc(A,B))
#define cuda_Memcpy(A,B,C,D) gpuErrchk(cudaMemcpy(A,B,C,D))
#define cuda_Free(A) gpuErrchk(cudaFree(A))
#define cuda_MemcpyToSymbol(A,B,C) gpuErrchk(cudaMemcpyToSymbol(A,B,C,cudaMemcpyHostToDevice))

#endif
