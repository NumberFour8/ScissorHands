#infdef ELLIPTIC_H
#define ELLIPTIC_H

#include "../def.h"

class ReducedPoint {
public:
	mpz_t X,Y;
	
	ReducedPoint(mpz_t x,mpz_t y)
	{
		mpz_init_set(X,x);
		mpz_init_set(Y,y);
	}
	
	~ReducedPoint()
	{
		mpz_clrs(X,Y);
	}
};

class EllipticCurve {
private:
	mpz_t a1,a2,a3,a4,a6,N;
	
public:

	EllipticCurve(mpz_t A1,mpz_t A2,mpz_t A3,mpz_t A4,mpz_t A6,mpz_t ch);
	~EllipticCurve()
	{ mpz_clrs(a1,a2,a3,a4,a6,N); }
	
	ReducedPoint addPoints(ReducedPoint P,ReducedPoint Q);
	ReducedPoint doublePoint(ReducedPoint P);
	
};

#endif
