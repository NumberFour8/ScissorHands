#ifndef ZINT_H
#define ZINT_H


#include "def.h"

class Zint {
private:
		mpz_t X;
public:
	Zint()
	{
		mpz_init_set_ui(X,0);
	}

	Zint(const unsigned int A)
	{
		mpz_init_set_ui(X,A);
	}

	Zint(const mpz_t A)
	{
		mpz_init_set(X,A);
	}
	
	~Zint()
	{
		mpz_clear(X);
	}
	
	mpz_t get()
	{
		return X;
	}
	
	Zint& invert_mod(const Zint& N)
	{
		if (!try_invert_mod(X,X,N.X))
		{
			throw inv;
		}
		return *this;
	}
	
	Zint& operator=(const Zint& Y)
	{
		mpz_set(X,Y.X);
	}	
	
	
	Zint& operator+=(const Zint& Y)
	{
		mpz_add(X,X,Y.X);
		return *this;
	}
	
	Zint& operator-=(const Zint& Y)
	{
		mpz_sub(X,X,Y.X);
		return *this;
	}
	
	Zint& operator*=(const Zint& Y)
	{
		mpz_mul(X,X,Y.X);
		return *this;
	}
	
	Zint& operator%=(const Zint& Y)
	{
		mpz_mod(X,X,Y.X);
		return *this;
	}
	
	Zint& operator+=(const unsigned Y)
	{
		mpz_add_ui(X,X,Y);
		return *this;
	}
	
	Zint& operator-=(const unsigned Y)
	{
		mpz_sub_ui(X,X,Y);
		return *this;
	}
	
	Zint& operator*=(const unsigned Y)
	{
		mpz_mul_ui(X,X,Y);
		return *this;
	}
	
	const Zint operator+(const Zint& Y)
	{
		Zint R = *this;
		R += Y;
		return R;
	}
	
	const Zint operator-(const Zint& Y)
	{
		Zint R = *this;
		R -= Y;
		return R;
	}
	
	const Zint operator*(const Zint& Y)
	{
		Zint R = *this;
		R *= Y;
		return R;
	}
	
	const Zint operator%(const Zint& Y)
	{
		Zint R = *this;
		R %= Y;
		return R;
	}
	
	const Zint operator+(const unsigned Y)
	{
		Zint R = *this;
		R += Y;
		return R;
	}
	
	const Zint operator-(const unsigned Y)
	{
		Zint R = *this;
		R -= Y;
		return R;
	}
	
	const Zint operator*(const unsigned Y)
	{
		Zint R = *this;
		R *= Y;
		return R;
	}
	
};

#endif
