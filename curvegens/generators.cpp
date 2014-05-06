#include "generators.h"

CurveGenerator::CurveGenerator(mpz_t n) 
 : s(1), edwards(0), twisted(0), A(0) 
{ 
	mpz_init_set(N,n);
}

CurveGenerator::~CurveGenerator()
{
	mpz_clear(N);
}

int CurveGenerator::getCoeff() 
{
	return s;
}

int CurveGenerator::getA()  
{
	return A;
}

int CurveGenerator::countEdwards() 
{
	return edwards;
}

int CurveGenerator::countTwisted() 
{ 
	return twisted;
}

void CurveGenerator::restart()
{
	A = s = edwards = twisted = 0;
	reset();
}

bool CurveGenerator::next_base_point(ReducedPoint& P)
{
	if (next(P) && A*A == 1)
	{
		if (A == 1) edwards++;
		else 		twisted++;
		
		return true;
	}
	return false;
}
