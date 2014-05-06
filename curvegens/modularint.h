#ifndef MODULARINT_H
#define MODULARINT_H


#include "def.h"

class Zmod {
private:
		mpz_t X,N;
public:
	Zmod(mpz_t A,mpz_t n)
	{
		mpz_init_set(X,A);
		mpz_init_set(N,n);
	}
	
	~Zmod()
	{
		mpz_clrs(A,n);
	}
	
	
	
};

#endif
