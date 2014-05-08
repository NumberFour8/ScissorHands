#ifndef ELLIPTIC_H
#define ELLIPTIC_H

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

	ReducedPoint(mpz_t x,mpz_t y)
	{
		set(x,y);
	}
	
	ReducedPoint(string x,string y,const mpz_t N)
	{
		mpq_t qx,qy;
		mpq_intz(qx,qy);
		mpz_intz(X,Y);
		
		mpq_set_str(qx,x.c_str(),10);
		mpq_set_str(qy,y.c_str(),10);
		
		reduce_rational_point(X,Y,qx,qy,N);
		mpq_clrs(qx,qy);
	}
	
	~ReducedPoint()
	{
		mpz_clrs(X,Y);
	}
};

class EllipticCurve {
private:
	Zint a1,a2,a3,a4,N;
	
	void add(ReducedPoint& R,const Zint& L,const ReducedPoint& P,const ReducedPoint& Q) const
	{
		Zint x1 = Zint(P.X);
		Zint x3 = (L^2)+a1*L-a2-x1-Zint(Q.X);
		Zint y3 = L*(x1-x3)-Zint(P.Y)-a1*x3-a3;
		
		x3 %= N;
		y3 %= N;
		
		if (x3 < 0) x3 += N;
		if (y3 < 0) y3 += N;
		
		R.set(x3,y3);
	}
	
public:

	EllipticCurve(string A1,string A2,string A3,string A4,const mpz_t n)
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
	
	void addPoints(ReducedPoint& R,ReducedPoint& P,ReducedPoint& Q) const
	{
		Zint L = Zint(P.X)-Zint(Q.X);
		L.invert_mod(N);
		L *= (Zint(P.Y)-Zint(Q.Y));
		
		add(R,L % N,P,Q);
	}
	
	void doublePoint(ReducedPoint& R,ReducedPoint& P) const
	{
		Zint x1(P.X),y1(P.Y);
		
		Zint L = y1*2+a1*x1+a3;
		L.invert_mod(N);
		L *= (x1^2)*3+a2*x1*2+a4-a1*y1;
		
		add(R,L % N,P,P);
	}
	
};

#endif
