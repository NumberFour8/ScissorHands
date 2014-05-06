#ifndef ZINT_H
#define ZINT_H


#include "../def.h"
#include "../helpers.h"

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
	
	void get(mpz_t x)
	{
		mpz_set(x,X);
	}
	
	Zint& invert_mod(Zint& N)
	{
		if (!try_invert_mod(X,X,N.X))
		{
			throw X;
		}
		return *this;
	}
	
	Zint& operator=(const Zint& Y)
	{
		mpz_set(X,Y.X);
		return *this;
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
};

	inline Zint operator+(Zint R,const Zint& Y)
	{
		R += Y;
		return R;
	}
	
	inline Zint operator-(Zint R,const Zint& Y)
	{
		R -= Y;
		return R;
	}
	
	inline Zint operator*(Zint R,const Zint& Y)
	{
		R *= Y;
		return R;
	}
	
	inline Zint operator%(Zint R,const Zint& Y)
	{
		R %= Y;
		return R;
	}
	
	inline Zint operator+(Zint R,const unsigned Y)
	{
		R += Y;
		return R;
	}
	
	inline Zint operator-(Zint R,const unsigned Y)
	{
		R -= Y;
		return R;
	}
	
	inline Zint operator*(Zint R,const unsigned Y)
	{
		R *= Y;
		return R;
	}

#endif
