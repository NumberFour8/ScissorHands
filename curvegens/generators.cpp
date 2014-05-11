#include "generators.h"

Generator::Generator() 
 : S(1), edwards(0), twisted(0), A(0) 
{ 
	
}

int Generator::getCoeff() 
{
	return S;
}

int Generator::getA()  
{
	return A;
}

int Generator::countEdwards() 
{
	return edwards;
}

int Generator::countTwisted() 
{ 
	return twisted;
}

int Generator::countCurves()
{
	return edwards+twisted;
}

void Generator::restart()
{
	edwards = twisted = 0;
	reset();
}

bool Generator::next_base_point(ReducedPoint& P,const mpz_t zN)
{
	if (next(P,zN) && A*A == 1)
	{
		if (A == 1) edwards++;
		else 		twisted++;
			
		return true;
	}
	else return false;
}

CurveGenerator::CurveGenerator(Torsion t,unsigned int from,unsigned int b)
	: Generator(), T(t), burst(b), C(NULL), G(NULL), curveCounter(0)
{
	C->doublePoint(Q,*G);
	for (S = 2;S < from-1;S++)
	   C->addPoints(Q,Q,*G);
	R.set(Q);
}

CurveGenerator::~CurveGenerator()
{
	if (C != NULL) delete C;
	if (G != NULL) delete G;
	C = NULL;
	G = NULL;
}

void CurveGenerator::reset()
{
	curveCounter = 0;
}

void CurveGenerator::revert()
{
	R.set(Q);
}

int CurveGenerator::getCurveNumber()
{
	return curveCounter;
}

bool CurveGenerator::next(ReducedPoint& P,const mpz_t zN)
{
	if (G == NULL || C == NULL) return false;
	
	if (curveCounter == burst) return false;
		
	C->addPoints(R,R,*G);
	
	generate_base_point(P,zN);
	curveCounter++;
	S++;
	return true;
}
