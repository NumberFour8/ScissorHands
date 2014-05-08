#include "generators.h"

EdwardsGenerator::EdwardsGenerator(mpz_t n,Torsion t,unsigned int from,unsigned int b)
 : CurveGenerator(t,from,b)
 {
	switch (T)
	{
		case Z12:
			C = new EllipticCurve("0","0","0","-12",n);
			G = new ReducedPoint("6","-12",n);
			A = 1;
		break;
		case Z2xZ8:
			C = new EllipticCurve("0","0","0","-8",n);
			G = new ReducedPoint("12","40",n);
			A = 1;
		break;
		case Z6:
			C = new EllipticCurve("0","-1/2304","0","-5/221184",n);
			G = new ReducedPoint("1/192","1/4608",n);
			A = -1;
		break;
		case Z8:
			C = new EllipticCurve("0","0","0","48",n);
			G = new ReducedPoint("4","-16",n);
			A = -1;
		break;
		case Z2xZ4:
			C = new EllipticCurve("0","0","0","-11664/25",n);
			G = new ReducedPoint("36","-864/5",n);
			A = -1;
		break;
	};
 }
 
void EdwardsGenerator::generate_base_point(ReducedPoint& P,mpz_t zN)
{	
	Zint s(Q.X),t(Q.Y),x,y;
	
	if (T == Z12)
	{
		x = (s-2)*(s+6)*((s^2)+12*s-12);
		x.invert_mod(N);
		x *= 8*t*((s^2)+12);
		
		y = (s-2)*(s+6)*((s^2)-12);
		y.invert_mod(N);
		y *= 4*s*((s^2)-12*s-12);
		y  = -y;	
	}
	else if (T == Z2xZ8)
	{
		Zint a = t+25;
		a.invert_mod(zN);
		a *= s-9;
		a += 1;
		a.invert_mod(zN);
		
		Zint b = 8*(a^2)-1;
		b.invert_mod(zN);
		b *= 2*a*(4*a+1);
	
		x  = 6*b-5;
		x.invert_mod(zN);
		x *= (2*b-1)*(4*b-3);
		
		y  = (t+3*s-2)*(t+s+16);
		y.invert_mod(zN);
		y *= (2*b-1)*((t^2)+50*t-2*(s^3)+27*(s^2)-104);
	}
	else if (T == Z6)
	{
		Zint sig = 96*s;
		sig.invert_mod(zN);
		sig *= (1-96*s);
		
		Zint r  = s^2;
		r.invert_mod(zN);
		r *= t;
		
		Zint al = (sig^2)-5;
		Zint bt = 4*sig;

		x = (sig-1)*(sig+5)*(sig*sig+5);
		x.invert_mod(zN);
		x *= 2*r*sig;
		
		y = (al^3)+(bt^3);
		y.invert_mod(zN);
		y *= (al^3)-(bt^3);
	}
	else if (T == Z8)
	{
		Zint u = t;
		u.invert_mod(zN);
		u *= 2*s;
		
		Zint v = t^2;
		v.invert_mod(zN);
		v *= 2*(s^3)-(t^2);
	
		x = 2*(u^2);
		
		y = v;
		y.invert_mod(zN);
		y *= 4*(u^4)-1;
	}
	else if (T == Z2xZ4)
	{
		x = 3;
		x.invert_mod(zN);
		
		y = 625*(t^2);
		y = -y + 77760*t+1250*(s^3)+24300*(s^2)-3779136;
		y.invert_mod(zN);
		y *= (30*s-25*t+1944)^2;
	}
	else return;
	
	y %= zN;
	x %= zN;
		
	if (x < 0) x += zN;
	if (y < 0) y += zN;
		
	P.set(x,y);
}
