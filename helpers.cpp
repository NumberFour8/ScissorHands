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
  for (int i = 0;i < NB_DIGITS;i++)
  {
	#if mp_bits_per_limb == 32
		a[i] = mpz_getlimbn(b, i);  
	#else 
		if (i%2 == 0)
		  a[i] = (mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
		else
		  a[i] = (mpz_getlimbn(b, i/2) >> 32);  
	#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui(a, 0);

  for (i = NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp(a, a, 32);
	mpz_add_ui(a, a, b[i]);
  }
}

void printmpz(const char* format,mpz_t number)
{
	char out[2048];
	memset((void*)out,0,2048);
	printf(format,mpz_get_str((char*)out,10,number));
}

void printBigInt(biguint_t B)
{
	for (int i = 0;i < NB_DIGITS;++i)
	{
		printf("%#010x",B[i]);
		if (B[i+1] == 0) break;
		if (i != NB_DIGITS-1) printf(",");
	}
	printf("\n");
}

NAF::NAF(unsigned char W,mpz_t number) : w(W)
{
	size_t sz = mpz_sizeinbase(exp,2);
	bits = (char*)malloc(sz*8+1);
	l = 0;

	unsigned int pow2 = 1 << w,rem,pww = 1 << (w-1);

	for (int i = 0;mpz_sgn(number) > 0;++i)
	{
		if (mpz_odd_p(number))
		{
			rem = (unsigned int)mpz_fdiv_ui(number,pow2);
			if (rem&pww)
			{
				bits[i] = (unsigned char)(rem - pow2);
			}
			else {
				bits[i] = (unsigned char)rem;
			}
			if (bits[i] < 0)
			  mpz_add_ui(number,number,(int)(bits[i]*(-1)));
			else mpz_sub_ui(number,number,(int)(bits[i]));
			l = i;
		}
		else bits[i] = 0;
		mpz_div_2exp(number,number,1);
	}

	l += 1;
	bits = (char*)realloc((void*)bits,(size_t)l*8);
}

int NAF::build(unsigned int start,unsigned int end) const
{
	int ret = 0;
	for (int i = start;i <= end;i++)
	{
		ret += bits[i]*(1 << (i-start));
	}

	return ret;
}

NAF::~NAF()
{
	free((void*)bits);
}
