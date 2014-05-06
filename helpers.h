#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "def.h"

#define mpz_intz(...) mpz_inits(__VA_ARGS__,NULL)
#define mpz_clrs(...) mpz_clears(__VA_ARGS__,NULL)

#define mpq_intz(...) mpq_inits(__VA_ARGS__,NULL)
#define mpq_clrs(...) mpq_clears(__VA_ARGS__,NULL)

// Rozvoj cisla do non-adjacent form (NAF)
class NAF {
public:
	// Koeficienty rozvoje
	char* bits;
	
	// Delka rozvoje
	unsigned long l;
	
	// Sirka rozvoje
	unsigned char w;
	
	// Vychozi­ konstruktor
	NAF() : w(2), bits(NULL), l(0)
	{ }
	
	virtual ~NAF();

	// Vytvori­ NAF rozvoj ci­sla N dane delky W
	void initialize(mpz_t N,unsigned char W = 2);
	
	// Vypise rozvoj
	void print() const;
	
};

// Preevede X do Montgomeryho reprezentace modulo N
void to_mont_repr(mpz_t x, mpz_t n);

// Preevede X z Montgomeryho reprezentace modulo N
void from_mont_repr(mpz_t x, mpz_t n,mpz_t invB);

// Preevede MPZ cislo do baze 2^32
void mpz_to_biguint(biguint_t a, mpz_t b);

// Preevede ci­slo z baze 2^32 do MPZ
void biguint_to_mpz(mpz_t a, biguint_t b);

// Prevede MPZ cislo do retezce
std::string mpz_to_string(mpz_t number);

// Vypise ci­slo v bazi 2^32
void printBigInt(const char* tag,biguint_t B);

// Spocita LCM(1,2,3...,n)
void lcmToN(mpz_t res,const unsigned int n);

/* Pokusi se invertovat X modulo N.
   Pokud inverze neexistuje vraci false a v invX je faktor N.
   Je-li gcd(X,N) = N, pak vraci­ 0.
 */ 
bool try_invert_mod(mpz_t invx,mpz_t x,mpz_t N);

/* Vypocita redukci racionalni­ho ci­sla q modulo n.
   Pri chybe vyhodi vyjimku­ s faktorem ci­sla N nebo 0. 
*/
void reduce_mod(mpz_t r,mpq_t q,mpz_t n);

// Redukuje racionalni souradnice modulo N. Stejne chovani jako reduce_mod
void reduce_rational_point(mpz_t X,mpz_t Y,mpq_t Qx,mpq_t Qy,const mpz_t N);

// Vrati true, je-li n prvocislo s pravdepobnosti mensi nez 2^(-50)
inline bool is_almost_surely_prime(mpz_t n)
{
	return (mpz_probab_prime_p(n,25) != 0);
}

// Vynuluje cislo v bazi 2^32
inline void reset(biguint_t n)
{
	memset((void*)n,0,MAX_BYTES);
}

// Pomocna tri­da konfigurace vypoctu
class ComputeConfig {
public:
	biguint_t N;   // N
	biguint_t N3;  // 3*N
	digit_t invN;  // N^(-1) mod W
	
	unsigned short windowSz;	 // Velikost okna
	unsigned long  nafLen;		 // Sirka NAF
	unsigned long  numCurves;	 // Pocet krivek
	int			   deviceId;	 // ID zarizeni pro vypocet
	bool		   minus1;		 // Ktere krivky pouzit v single rezimu?

	mutable float  cudaRunTime;  // Po dobehnuti je hodnota nastavena na cas
	void initialize(mpz_t zN);
};

// Vypise chybu GPU a prerusi program
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//////////////////////// POMOCNA MAKRA PRO CUDA ///////////////////////////

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define cuda_Malloc(A,B) gpuErrchk(cudaMalloc(A,B))
#define cuda_Memcpy(A,B,C,D) gpuErrchk(cudaMemcpy(A,B,C,D))
#define cuda_Free(A) gpuErrchk(cudaFree(A))
#define cuda_MemcpyToSymbol(A,B,C) gpuErrchk(cudaMemcpyToSymbol(A,B,C,cudaMemcpyHostToDevice))


#define START_MEASURE(start) gpuErrchk(cudaEventRecord(start,0))
#define STOP_MEASURE(name,start,stop,tt) {gpuErrchk(cudaEventRecord(stop,0));\
								  gpuErrchk(cudaEventSynchronize(stop));\
								  float time = 0;\
								  gpuErrchk(cudaEventElapsedTime(&time,start,stop));\
								  tt += time; \
								  time > 2000.0f ? printf("%s : %.3f s\n",name,time/1000.0f) :\
												   printf("%s : %.3f ms\n",name,time); } 
#endif
