#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "zint.h"
#include "../helpers.h"

// Racionalni bod elipticke krivky
class RationalPoint {
public:
	Qrac X,Y;
	
	RationalPoint() 
	  : X(), Y()
	{ }
	
	RationalPoint(string x,string y)
	  : X(x), Y(y)
	{ }
	
	inline void set(mpq_t x,mpq_t y)
	{
		X = x;
		Y = y;
	}
	
	inline void set(const RationalPoint& K)
	{
		X = K.X;
		Y = K.Y;
	}	
}

// Bod elipticke krivky redukovany z racionalniho bodu modulo N
class ReducedPoint {
public:
	mpz_t X,Y;
	
	ReducedPoint()
    {
		mpz_intz(X,Y);
	}

	void set(mpz_t x,mpz_t y)
	{
		mpz_init_set(X,x);
		mpz_init_set(Y,y);
	}
	
	void set(Zint x,Zint y)
	{
		x.get(X);
		y.get(Y);
	}

	ReducedPoint(const RationalPoint& Q,const mpz_t N)
	{
		mpz_intz(X,Y);
		reduce_rational_point(X,Y,Q.X,Q.Y,N);	
	}

	ReducedPoint(mpz_t x,mpz_t y)
	{
		set(x,y);
	}
	
	~ReducedPoint()
	{
		mpz_clrs(X,Y);
	}
};


// Predstavuje obecnou eliptickou krivku nad Q
class EllipticCurve {
private:
	Qrac a1,a2,a3,a4;
	
	void add(RationalPoint& R,const Qrac& L,const RationalPoint& P,const RationalPoint& Q) const
	{
		Qrac x3 = L*L+a1*L-a2-P.X-Q.X;
		Qrac y3 = L*(P.X-x3)-P.Y-a1*x3-a3;
		
		R.set(x3,y3);
	}
	
public:

	EllipticCurve(string A1,string A2,string A3,string A4)
	  : a1(A1),a2(A2),a3(A3),a4(A4)
	{ 
	
	}
	
	void addPoints(RationalPoint& R,RationalPoint& P,RationalPoint& Q) const
	{
		Qrac L = (P.Y-Q.Y)/(P.X-Q.X);
		add(R,L,P,Q);
	}
	
	void doublePoint(RationalPoint& R,RationalPoint& P) const
	{
		Qrac L = (P.X*P.X*3+a2*P.X*2+a4-a1*P.Y)/(P.Y*2+a1*P.X+a3);
		add(R,L,P,P);
	}
	
};

#endif
