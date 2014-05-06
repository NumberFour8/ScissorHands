#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "../def.h"
#include "zint.h"
#include "../helpers.h"

class ReducedPoint {
public:
	mpz_t X,Y;
	
	ReducedPoint()
    {
		mpz_init_set_ui(X,0);
		mpz_init_set_ui(Y,0);
	}

	ReducedPoint(mpz_t x,mpz_t y)
	{
		mpz_init_set(X,x);
		mpz_init_set(Y,y);
	}
	
	ReducedPoint(string x,string y,mpz_t N)
	{
		mpq_t qx,qy;
		mpq_intz(qx,qy);
		mpq_set_str(qx,x.c_str(),10);
		mpq_set_str(qy,y.c_str(),10);
		
		reduce_rational_point(X,Y,qx,qy,N);
	}
	
	ReducedPoint(Zint x,Zint y)
	{
		mpz_intz(X,Y);
		x.get(X);
		y.get(Y);
	}
	
	~ReducedPoint()
	{
		mpz_clrs(X,Y);
	}
};

class EllipticCurve {
private:
	Zint a1,a2,a3,a4,N;
	
	ReducedPoint add(Zint L,ReducedPoint P,ReducedPoint Q)
	{
		Zint x1 = Zint(P.X);
		Zint x3 = L*L+a1*L-a2-x1-Zint(Q.X);
		Zint y3 = L*(x1-x3)-Zint(P.Y)-a1*x3-a3;
		
		return ReducedPoint(x3 % N,y3 % N);
	}
	
public:

	EllipticCurve(string A1,string A2,string A3,string A4,mpz_t n)
	 : N(n)
	{ 
			mpq_t T;
			mpz_t I;
			mpq_init(T);
			mpz_init(I);
			
			mpq_set_str(T,A1.c_str(),10);
			reduce_mod(I,T,n);
			a1 = Zint(I);
			
			mpq_set_str(T,A2.c_str(),10);
			reduce_mod(I,T,n);
			a2 = Zint(I);
			
			mpq_set_str(T,A3.c_str(),10);
			reduce_mod(I,T,n);
			a3 = Zint(I);
			
			mpq_set_str(T,A4.c_str(),10);
			reduce_mod(I,T,n);
			a4 = Zint(I);
			
			mpq_clear(T);
			mpz_clear(I);
	}
	
	ReducedPoint addPoints(ReducedPoint P,ReducedPoint Q)
	{
		Zint L = Zint(P.X)-Zint(Q.X);
		L.invert_mod(N);
		L *= (Zint(P.Y)-Zint(Q.Y));
		
		return add(L,P,Q);
	}
	
	ReducedPoint doublePoint(ReducedPoint P)
	{
		Zint x1(P.X),y1(P.Y);
		
		Zint L = y1*2+a1*x1+a3;
		L.invert_mod(N);
		L *= x1*x1*3+a2*x1*2+a4-a1*y1;
		
		return add(L,P,P);
	}
	
};

#endif
