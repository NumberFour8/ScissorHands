#ifndef HELPERS_H
#define HELPERS_H

#include "def.h"

#define mpz_intz(...) mpz_inits(__VA_ARGS__,NULL)
#define mpz_clrs(...) mpz_clears(__VA_ARGS__,NULL)

#define mpq_intz(...) mpq_inits(__VA_ARGS__,NULL)
#define mpq_clrs(...) mpq_clears(__VA_ARGS__,NULL)

// Rozvoj è­sla do non-adjacent form (NAF)
class NAF {
public:
	// Koeficienty rozvoje
	char* bits;
	
	// Délka rozvoje
	unsigned int l;
	
	// Šíøka rozvoje
	const unsigned char w;
	
	// Výchozí konstruktor
	NAF(unsigned char W) : w(W), bits(NULL), l(0)
	{ }
	
	virtual ~NAF();

	// Vytvoøí NAF rozvoj èísla N dané šíøky W
	void initialize(mpz_t N);
	
	// Vypíše rozvoj
	void print() const;
	
	// Výsek z rozvoje daný parametry start a end
	int build(unsigned int start,unsigned int end) const;
};

// Pøevede X do Montgomeryho reprezentace modulo N
void to_mont_repr(mpz_t x, mpz_t n);

// Pøevede X z Montgomeryho reprezentace modulo N
void from_mont_repr(mpz_t x, mpz_t n,mpz_t invB);

// Pøevede MPZ èíslo do báze 2^32
void mpz_to_biguint(biguint_t a, mpz_t b);

// Pøevede èí­slo z báze 2^32 do MPZ
void biguint_to_mpz(mpz_t a, biguint_t b);

// Vypíše è­slo MPZ
void printmpz(const char* format,mpz_t number);

// Vypíše èíslo v bázi 2^32
void printBigInt(const char* tag,biguint_t B);

// Spoèítá LCM(1,2,3...,n)
void lcmToN(mpz_t res,const unsigned int n);

/* Pokusí se invertovat X modulo N.
   Pokud inverze neexistuje vrací false a v invX je faktor N.
   Je-li gcd(X,N) = N, pak vrací 0.
 */ 
bool try_invert_mod(mpz_t invx,mpz_t x,mpz_t N);

/* Vypoèítá redukci racionálního èísla q modulo n.
   Pøi chybì vrací false a v r je faktor èísla N nebo 0. 
*/
bool reduce_mod(mpz_t r,mpq_t q,mpz_t n);


// Vynuluje èí­slo v bázi 2^32
inline void reset(biguint_t n)
{
	memset((void*)n,0,MAX_BYTES);
}

// Pomocná tøí­da pro N, 3*N a inverzi N modulo velikost báze
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

// Vypíše chybu GPU a pøeruší­ program
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//////////////////////// POMOCNÁ MAKRA PRO CUDA ///////////////////////////

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define cuda_Malloc(A,B) gpuErrchk(cudaMalloc(A,B))
#define cuda_Memcpy(A,B,C,D) gpuErrchk(cudaMemcpy(A,B,C,D))
#define cuda_Free(A) gpuErrchk(cudaFree(A))
#define cuda_MemcpyToSymbol(A,B,C) gpuErrchk(cudaMemcpyToSymbol(A,B,C,cudaMemcpyHostToDevice))


#define START_MEASURE(start) gpuErrchk(cudaEventRecord(start,0))
#define STOP_MEASURE(name,start,stop) {gpuErrchk(cudaEventRecord(stop,0));\
								  gpuErrchk(cudaEventSynchronize(stop));\
								  float time = 0;\
								  gpuErrchk(cudaEventElapsedTime(&time,start,stop));\
								  time > 2000.0f ? printf("%s : %.3f s\n",name,time/1000.0f) :\
												   printf("%s : %.3f ms\n",name,time); }
#endif
