#ifndef QRAC_H
#define QRAC_H

#include "../helpers.h"

class Qrac {
private:
		mpq_t X;
public:
	
	// Constructors
	Qrac()
	{
		mpq_init(X);
	}
	
	Qrac(const Qrac& Y)
	{
		mpq_init_set(X,Y.X);
	}

	Qrac(const string& A)
	{
		mpq_init_set_str(X,A);
		mpq_canonicalize(X);
	}

	Qrac(const mpq_t A)
	{
		mpq_init_set(X,A);
		mpq_canonicalize(X);
	}
	
	~Qrac()
	{
		mpq_clear(X);
	}
	
	
	// Functions
	
	void get(mpq_t x)
	{
		mpq_set(x,X);
	}
	
	string str() const 
	{
		return mpq_to_string(X);
	}
	
	// Qrac operators
	
	Qrac& operator=(const Qrac& Y)
	{
		mpq_set(X,Y.X);
		return *this;
	}	
	
	
	Qrac& operator+=(const Qrac& Y)
	{
		mpq_add(X,X,Y.X);
		return *this;
	}
	
	Qrac& operator-=(const Qrac& Y)
	{
		mpq_sub(X,X,Y.X);
		return *this;
	}
	
	Qrac& operator*=(const Qrac& Y)
	{
		mpq_mul(X,X,Y.X);
		return *this;
	}
	
	Qrac& operator/=(const Qrac& Y)
	{
		mpq_div(X,X,Y.X);
		return *this;
	}
	
	// MPQ_T operators
	
	Qrac& operator=(const mpq_t Y)
	{
		mpq_set(X,Y);
		return *this;
	}	
	
	
	Qrac& operator+=(const mpq_t Y)
	{
		mpq_add(X,X,Y);
		return *this;
	}
	
	Qrac& operator-=(const mpq_t Y)
	{
		mpq_sub(X,X,Y);
		return *this;
	}
	
	Qrac& operator*=(const mpq_t Y)
	{
		mpq_mul(X,X,Y);
		return *this;
	}
	
	Qrac& operator/=(const mpq_t Y)
	{
		mpq_div(X,X,Y);
		return *this;
	}
	
	// Uint operators
	
	Qrac& operator=(const unsigned int Y)
	{
		mpq_set_ui(X,Y,1);
		return *this;
	}	
	
	Qrac& operator+=(const unsigned int Y)
	{
		mpq_t B;
		mpq_init_set_ui(B,Y,1);
		
		mpq_add(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator-=(const unsigned int Y)
	{
		mpq_t B;
		mpq_init_set_ui(B,Y,1);
		
		mpq_sub(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator*=(const unsigned int Y)
	{
		mpq_t B;
		mpq_init_set_ui(B,Y,1);
		
		mpq_mul(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator/=(const unsigned int Y)
	{
		mpq_t B;
		mpq_init_set_ui(B,Y,1);
		
		mpq_div(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator^=(const unsigned int Y)
	{
		mpz_t n = mpq_numref(X);
		mpz_t d = mpq_denref(X);
		
		mpz_pow_ui(b,b,Y);
		mpz_pow_ui(d,d,Y);
		
		return *this;
	}

	Qrac operator-()
	{
		mpq_t n;
		mpq_init(n);
		mpq_neg(n,X);
	
		Qrac ret(n);
		mpq_clear(n);
		return ret;
	}
	
	friend bool operator< (const Qrac& lhs, const Qrac& rhs);
	friend bool operator< (const Qrac& lhs, const unsigned int rhs);
	friend bool operator> (const Qrac& lhs, const Qrac& rhs);
	friend bool operator> (const Qrac& lhs, const unsigned int rhs);
};

// Comparison operators

inline bool operator< (const Qrac& lhs, const Qrac& rhs)
{
	return mpq_cmp(lhs.X,rhs.X) < 0;
}

inline bool operator< (const Qrac& lhs, const unsigned int rhs)
{
	return mpq_cmp_ui(lhs.X,rhs) < 0;
}

inline bool operator> (const Qrac& lhs, const Qrac& rhs)
{
	return rhs < lhs;
}

inline bool operator> (const Qrac& lhs, const unsigned int rhs)
{
	return mpq_cmp_ui(lhs.X,rhs) > 0;
}

// Qrac operators

inline Qrac operator+(Qrac R,const Qrac& Y)
{
	R += Y;
	return R;
}

inline Qrac operator-(Qrac R,const Qrac& Y)
{
	R -= Y;
	return R;
}

inline Qrac operator*(Qrac R,const Qrac& Y)
{
	R *= Y;
	return R;
}

inline Qrac operator/(Qrac R,const Qrac& Y)
{
	R /= Y;
	return R;
}

// MPZ_T operators

inline Qrac operator+(Qrac R,const mpq_t Y)
{
	R += Y;
	return R;
}

inline Qrac operator-(Qrac R,const mpq_t Y)
{
	R -= Y;
	return R;
}

inline Qrac operator*(Qrac R,const mpq_t Y)
{
	R *= Y;
	return R;
}

inline Qrac operator/(Qrac R,const mpq_t Y)
{
	R %= Y;
	return R;
}

// Uint operators

inline Qrac operator+(Qrac R,const unsigned int Y)
{
	R += Y;
	return R;
}

inline Qrac operator-(Qrac R,const unsigned int Y)
{
	R -= Y;
	return R;
}

inline Qrac operator*(Qrac R,const unsigned int Y)
{
	R *= Y;
	return R;
}


inline Qrac operator/(Qrac R,const unsigned int Y)
{
	R *= Y;
	return R;
}

inline Qrac operator^(Qrac R,const unsigned int Y)
{
	R ^= Y;
	return R;
}

#endif
