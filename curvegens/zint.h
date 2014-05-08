#ifndef ZINT_H
#define ZINT_H

#include "../helpers.h"

class Zint {
private:
		mpz_t X;
public:
	
	// Constructors
	Zint()
	{
		mpz_init_set_ui(X,0);
	}
	
	Zint(const Zint& Y)
	{
		mpz_init_set(X,Y.X);
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
	
	
	// Functions
	
	void get(mpz_t x)
	{
		mpz_set(x,X);
	}
	
	Zint& invert_mod(const Zint& N)
	{
		if (!try_invert_mod(X,X,N.X))
		{
			throw_factor(X);
		}
		return *this;
	}

	Zint& invert_mod(const mpz_t N)
	{
		if (!try_invert_mod(X,X,N))
		{
			throw_factor(X);
		}
		return *this;
	}
	
	string str() const 
	{
		return mpz_to_string(X);
	}
	
	// Zint operators
	
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
	
	// MPZ_T operators
	
	Zint& operator=(const mpz_t Y)
	{
		mpz_set(X,Y);
		return *this;
	}	
	
	
	Zint& operator+=(const mpz_t Y)
	{
		mpz_add(X,X,Y);
		return *this;
	}
	
	Zint& operator-=(const mpz_t Y)
	{
		mpz_sub(X,X,Y);
		return *this;
	}
	
	Zint& operator*=(const mpz_t Y)
	{
		mpz_mul(X,X,Y);
		return *this;
	}
	
	Zint& operator%=(const mpz_t Y)
	{
		mpz_mod(X,X,Y);
		return *this;
	}
	
	// Uint operators
	
	Zint& operator=(const unsigned int Y)
	{
		mpz_set_ui(X,Y);
		return *this;
	}	
	
	Zint& operator+=(const unsigned int Y)
	{
		mpz_add_ui(X,X,Y);
		return *this;
	}
	
	Zint& operator-=(const unsigned int Y)
	{
		mpz_sub_ui(X,X,Y);
		return *this;
	}
	
	Zint& operator*=(const unsigned int Y)
	{
		mpz_mul_ui(X,X,Y);
		return *this;
	}
	
	Zint& operator^=(const unsigned int Y)
	{
		mpz_pow_ui(X,X,Y);
		return *this;
	}

	Zint operator-()
	{
		mpz_t n;
		mpz_init(n);
		mpz_neg(n,X);
	
		Zint ret(n);
		mpz_clear(n);
		return ret;
	}
	
	friend bool operator< (const Zint& lhs, const Zint& rhs);
	friend bool operator< (const Zint& lhs, const unsigned int rhs);
	friend bool operator> (const Zint& lhs, const Zint& rhs);
	friend bool operator> (const Zint& lhs, const unsigned int rhs);
};

// Comparison operators

inline bool operator< (const Zint& lhs, const Zint& rhs)
{
	return mpz_cmp(lhs.X,rhs.X) < 0;
}

inline bool operator< (const Zint& lhs, const unsigned int rhs)
{
	return mpz_cmp_ui(lhs.X,rhs) < 0;
}

inline bool operator> (const Zint& lhs, const Zint& rhs)
{
	return rhs < lhs;
}

inline bool operator> (const Zint& lhs, const unsigned int rhs)
{
	return mpz_cmp_ui(lhs.X,rhs) > 0;
}

// Zint operators

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

// MPZ_T operators

inline Zint operator+(Zint R,const mpz_t Y)
{
	R += Y;
	return R;
}

inline Zint operator-(Zint R,const mpz_t Y)
{
	R -= Y;
	return R;
}

inline Zint operator*(Zint R,const mpz_t Y)
{
	R *= Y;
	return R;
}

inline Zint operator%(Zint R,const mpz_t Y)
{
	R %= Y;
	return R;
}

// Uint operators

inline Zint operator+(Zint R,const unsigned int Y)
{
	R += Y;
	return R;
}

inline Zint operator-(Zint R,const unsigned int Y)
{
	R -= Y;
	return R;
}

inline Zint operator*(Zint R,const unsigned int Y)
{
	R *= Y;
	return R;
}

inline Zint operator^(Zint R,const unsigned int Y)
{
	R ^= Y;
	return R;
}

#endif
