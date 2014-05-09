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
	: Generator(), T(t), start(from), burst(b), C(NULL), G(NULL), curveCounter(0)
{
	
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
	if (S > 1) S = 2;
	curveCounter = 0;
}

bool CurveGenerator::next(ReducedPoint& P,const mpz_t zN)
{
	if (G == NULL || C == NULL) return false;
	
	if (curveCounter == burst) return false;
	
	if (S <= 1)
	{
	  C->doublePoint(Q,*G);
	  S = 2;
	}
	
	for (;S < start-1;S++)
	   C->addPoints(Q,Q,*G);
		
	if (S > 1) C->addPoints(Q,Q,*G);
	
	generate_base_point(P,zN);
	S++;
	curveCounter++;
	return true;
}
