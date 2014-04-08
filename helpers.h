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

#endif
