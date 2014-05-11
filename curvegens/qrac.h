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
		mpq_init(X);
		mpq_set(X,Y.X);
		mpq_canonicalize(X);
	}

	Qrac(const string& A)
	{
		mpq_init(X);
		mpq_set_str(X,A.c_str(),10);
		mpq_canonicalize(X);
	}

	Qrac(const mpq_t A)
	{
		mpq_init(X);
		mpq_set(X,A);
		mpq_canonicalize(X);
	}
	
	~Qrac()
	{
		mpq_clear(X);
	}
	
	
	// Functions
	
	mpq_ptr get()
	{
		return X;
	}

	void invert()
	{
		mpq_inv(X,X);
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
	
	Qrac& operator=(const char* Y)
	{
		mpq_set_str(X,Y,10);
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
	
	// int operators
	
	Qrac& operator=(const int Y)
	{
		mpq_set_si(X,Y,1);
		return *this;
	}	
	
	Qrac& operator+=(const int Y)
	{
		mpq_t B;
		mpq_init(B);
		mpq_set_si(B,Y,1);
		
		mpq_add(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator-=(const int Y)
	{
		mpq_t B;
		mpq_init(B);
		mpq_set_si(B,Y,1);
		
		mpq_sub(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator*=(const int Y)
	{
		mpq_t B;
		mpq_init(B);
		mpq_set_si(B,Y,1);
		
		mpq_mul(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator/=(const int Y)
	{
		mpq_t B;
		mpq_init(B);
		mpq_set_si(B,Y,1);
		
		mpq_div(X,X,B);
		
		mpq_clear(B);
		return *this;
	}
	
	Qrac& operator^=(const unsigned int Y)
	{	
		mpz_pow_ui(mpq_numref(X),mpq_numref(X),Y);
		mpz_pow_ui(mpq_denref(X),mpq_denref(X),Y);
		
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
	return mpq_cmp_ui(lhs.X,rhs,1) < 0;
}

inline bool operator> (const Qrac& lhs, const Qrac& rhs)
{
	return rhs < lhs;
}

inline bool operator> (const Qrac& lhs, const unsigned int rhs)
{
	return mpq_cmp_ui(lhs.X,rhs,1) > 0;
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
	R /= Y;
	return R;
}

// Int operators

inline Qrac operator+(Qrac R,const int Y)
{
	R += Y;
	return R;
}

inline Qrac operator+(const int Y,Qrac R)
{
	R += Y;
	return R;
}

inline Qrac operator-(Qrac R,const int Y)
{
	R -= Y;
	return R;
}

inline Qrac operator-(const int Y,Qrac R)
{
	R = -R + Y;
	return R;
}

inline Qrac operator*(Qrac R,const int Y)
{
	R *= Y;
	return R;
}

inline Qrac operator*(const int Y,Qrac R)
{
	R *= Y;
	return R;
}

inline Qrac operator/(Qrac R,const int Y)
{
	R /= Y;
	return R;
}

inline Qrac operator/(const int Y,Qrac R)
{
	R.invert();
	R *= Y;
	return R;
}

inline Qrac operator^(Qrac R,const unsigned int Y)
{
	R ^= Y;
	return R;
}

#endif
