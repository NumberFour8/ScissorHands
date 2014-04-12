#ifndef HELPERS_H
#define HELPERS_H

#include "def.h"

#define mpz_intz(...) mpz_inits(__VA_ARGS__,NULL)
#define mpz_clrs(...) mpz_clears(__VA_ARGS__,NULL)

// T��da obsahuj�c� vysledek p�evodu z Extended sou�adnic
class ExtResult {
public:
	// Indik�tor, zda byl nalezen faktor
	bool factorFound;
	
	// Nalezen� faktor
	mpz_t factor;

	// V�choz� konstruktor
	ExtResult() : factorFound(false)
	{ mpz_init_set_ui(factor,0); }
	
	~ExtResult()
	{ mpz_clear(factor); }
};

// Rozvoj �sla do non-adjacent form (NAF)
class NAF {
public:
	// Koeficienty rozvoje
	char* bits;
	
	// D�lka rozvoje
	unsigned int l;
	
	// ���ka rozvoje
	const unsigned char w;
	
	// Zkontruuje rozvoj ���ky W z ��sla N
	NAF(unsigned char W,mpz_t N);
	
	virtual ~NAF();
	
	// Vyp�e rozvoj
	void print() const;
	
	// V�sek z rozvoje dan� parametry start a end
	int build(unsigned int start,unsigned int end) const;
};

// P�evede X do Montgomeryho reprezentace modulo N
void to_mont_repr(mpz_t x, mpz_t n);

// P�evede X z Montgomeryho reprezentace modulo N
void from_mont_repr(mpz_t x, mpz_t n,mpz_t invB);

// P�evede MPZ ��slo do b�ze 2^32
void mpz_to_biguint(biguint_t a, mpz_t b);

// P�evede ��slo z b�ze 2^32 do MPZ
void biguint_to_mpz(mpz_t a, biguint_t b);

// Vyp�e �slo MPZ
void printmpz(const char* format,mpz_t number);

// Vyp�e ��slo v b�zi 2^32
void printBigInt(const char* tag,biguint_t B);

// Spo��t� LCM(1,2,3...,n)
void lcmToN(mpz_t res,const unsigned int n);

// Vynuluje ��slo v b�zi 2^32
inline void reset(biguint_t n)
{
	memset((void*)n,0,MAX_BYTES);
}

// Pomocn� t��da pro N, 3*N a inverzi N modulo velikost b�ze
class Aux {
public:
	biguint_t N;
	biguint_t N3;
	digit_t invN;
	
	Aux(mpz_t zN)
	{
		reset(N);
		reset(N3);
		invN = 0;
	
		mpz_t z3N,zW,zInvN;
		mpz_intz(z3N,zW,zInvN);
	
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

// Vyp�e chybu GPU a p�eru�� program
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//////////////////////// POMOCNÁ MAKRA NA ZACHYCEN� CHYB V CUDA ////////////////////////

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define cuda_Malloc(A,B) gpuErrchk(cudaMalloc(A,B))
#define cuda_Memcpy(A,B,C,D) gpuErrchk(cudaMemcpy(A,B,C,D))
#define cuda_Free(A) gpuErrchk(cudaFree(A))
#define cuda_MemcpyToSymbol(A,B,C) gpuErrchk(cudaMemcpyToSymbol(A,B,C,cudaMemcpyHostToDevice))

#endif
