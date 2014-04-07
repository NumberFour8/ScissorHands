#include "def.h"
#include "helpers.h"

void to_mont_repr (mpz_t x, mpz_t n)
{
  mpz_mul_2exp (x, x, MAX_BITS);
  mpz_mod(x, x, n);
}

void from_mont_repr (mpz_t x, mpz_t n,mpz_t invB)
{
  mpz_mul(x, x, invB);
  mpz_mod(x, x, n);
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
  int i;

  for (i = 0;i < NB_DIGITS;i++)
  {
	#if mp_bits_per_limb == 32
		a[i]=mpz_getlimbn(b, i);  
	#else 
		if (i%2 == 0)
		  a[i]=(mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
		else
		  a[i]=(mpz_getlimbn(b, i/2) >> 32);  
	#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui(a, 0);

  for (i=NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp(a, a, 32);
	mpz_add_ui(a , a, b[i]);
  }
}

void reset(biguint_t r)
{
	memset((void*)r,0,MAX_BITS/8);
}

void printmpz(const char* format,mpz_t number)
{
	char out[2048];
	memset((void*)out,0,2048);
	printf(format,mpz_get_str((char*)out,10,number));
}