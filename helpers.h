#ifndef HELPERS_H
#define HELPERS_H

#include <mpir.h>

void to_mont_repr (mpz_t x, mpz_t n);
void from_mont_repr (mpz_t x, mpz_t n,mpz_t invB);

void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);

void reset(biguint_t r);
void printmpz(const char* format,mpz_t number);

#endif