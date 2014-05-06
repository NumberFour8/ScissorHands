#include "generators.h"
#include "../helpers.h"

EdwardsGenerator::EdwardsGenerator(mpz_t n,Torsion T)
 : CurveGenerator(n)
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


bool EdwardsGenerator::next(ReducedPoint& P)
{
	return false; 
}
