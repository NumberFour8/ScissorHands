#include "generators.h"
#include "helpers.h"

EdwardsGenerator::EdwardsGenerator(mpz_t n,EdwardsTorsion T)
 : CurveGenerator(n), A(1)
 {
	if (T == EdwardsTorsion::Z12)
	{
		E = new EllipticCurve("0","0","0","-12","0",n);
		G = new ReducedPoint("6","-12");
	}
	else if (T == EdwardsTorsion::Z2xZ8)
	{
		E = new EllipticCurve("0","0","0","-8","-32",n);
		G = new ReducedPoint("12","40");
	}
 }
 
 EdwardsGenerator::~EdwardsGenerator()
 {
	delete E;
	delete G;
 }


bool EdwardsGenerator::next(ReducedPoint& P)
{
	
}
