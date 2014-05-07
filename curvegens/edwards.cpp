#include "generators.h"
#include "../helpers.h"

EdwardsGenerator::EdwardsGenerator(mpz_t n,Torsion T,unsigned int from,unsigned int to)
 : CurveGenerator(n), tor(T), start(from), end(to)
 {
	A = 1;
	if (T == Torsion::Z12)
	{
		C = new EllipticCurve("0","0","0","-12",n);
		G = new ReducedPoint("6","-12",n);
	}
	else if (T == Torsion::Z2xZ8)
	{
		C = new EllipticCurve("0","0","0","-8",n);
		G = new ReducedPoint("12","40",n);
	}
 }
 
 EdwardsGenerator::~EdwardsGenerator()
 {
	delete C;
	delete G;
 }

void EdwardsGenerator::reset()
{
	s = 1; 
}

bool EdwardsGenerator::next(ReducedPoint& P)
{
	if (s > end) return false;
	
	if (s <= 1)
	{
	  C->doublePoint(Q,*G);
	  s = 2;
	}
	for (;s < start;s++)
	   C->addPoints(Q,Q,*G);
	
	if (tor == Torsion::Z12)
	{
		Zint u(Q.X);
		Zint t(Q.Y);
		
		try {
			Zint x = ((u-2)*(u+6)*((u^2)+12*u-12));
			x.invert_mod(N);
			x *= 8*t*((u^2)+12);
			
			Zint y = ((u-2)*(u+6)*((u^2)-12));
			y.invert_mod(N);
			y *= -4*u*((u^2)-12*u-12);
			
			P.set(x,y);
			++s;
		}
		catch (mpz_t f)
		{
			cout << "ERROR: Cannot reduce on curve #" << s << endl;
			if (mpz_cmp_ui(f,0) != 0) // Byl nalezen faktor?
			{
				cout << "Factor found: " << mpz_to_string(f) << endl;
			}
			return false;
		}
		return true;
	}
}
