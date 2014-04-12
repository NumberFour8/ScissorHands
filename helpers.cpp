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
	#if GMP_NUMB_BITS == 32
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

void printBigInt(const char* tag,biguint_t B)
{
	printf("%s: {",tag);
	for (int i = 0;i < NB_DIGITS;++i)
	{
		printf("%#010x",B[i]);
		if (B[i+1] == 0) break;
		if (i != NB_DIGITS-1) printf(",");
	}
	printf("}\n");
}

void lcmToN(mpz_t res,const unsigned int n)
{
	mpz_set_ui(res,1);
	for (unsigned int i = 2;i <= n;++i)
	{
		mpz_lcm_ui(res,res,i);
	}
}

bool try_invert_mod(mpz_t invx,mpz_t x,mpz_t N)
{
	bool ret = false;
		 
	mpz_t b,r;
	mpz_intz(b,r);
         
	mpz_gcdext(r,invx,b,x,N);
    if (mpz_cmp_ui(r,1) == 0)
	{
		if (mpz_sgn(invx) == -1) 
		  mpz_add(invx,invx,N);
		ret = true;
	}
    else if (mpz_cmp(r,N) == 0)
	{
		mpz_set_si(invx,0);
		printmpz("Error during inversion of: %s\n",x);
	}
	else 
	{
	    mpz_set(invx,r);
	    printmpz("Factor found: %s\n",invx);
	}

	mpz_clrs(r,b);
	
	return ret;
}

bool reduce_mod(mpz_t r,mpq_t q,mpz_t n)
{
	mpz_t den;
	mpz_init_set(den,mpq_denref(q));
	
	bool s = try_invert_mod(r,den,n);
	mpz_clear(den);

	if (!s) return false;
	
	mpz_mul(r,r,mpq_numref(q));
	mpz_mod(r,r,n);

	if (mpz_sgn(r) == -1) mpz_add(r,r,n);

	return true;
}

NAF::NAF(unsigned char W,mpz_t N) : w(W)
{
	mpz_t number;
	mpz_init_set(number,N);
	
	size_t sz = mpz_sizeinbase(number,2);
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
	mpz_clear(number);

	l += 1;
	bits = (char*)realloc((void*)bits,(size_t)l*8);
}

int NAF::build(unsigned int start,unsigned int end) const
{
	int ret = 0;
	for (unsigned int i = start;i <= end;i++)
	{
		ret += bits[i]*(1 << (i-start));
	}

	return ret;
}

void NAF::print() const
{
	printf("NAF%d: ",(int)w);
	for (int i = l-1;i >= 0;--i)
	{
		printf("%d",(int)bits[i]);
	}
	printf("\n");
}

NAF::~NAF()
{
	free((void*)bits);
}
