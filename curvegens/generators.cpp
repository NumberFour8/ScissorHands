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

void Generator::new_point_set()
{
	edwards = twisted = 0;
	reset();
}

void Generator::revert()
{
	edwards = twisted = 0;
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

CurveGenerator::CurveGenerator(Torsion t,unsigned int b)
	: Generator(), T(t), burst(b), C(NULL), G(NULL), curveCounter(0)
{ 
		cache = new RationalPoint*[burst];
		for (unsigned int i = 0;i < burst;i++)
			cache[i] = NULL;

		fromCache = false;
}

void CurveGenerator::initialize(unsigned int from)
{
	C->doublePoint(Q,*G);
	for (S = 2;S < from-1;S++)
	   C->addPoints(Q,Q,*G);
	
	R.set(Q);
	origS = S;
}

CurveGenerator::~CurveGenerator()
{
	if (C != NULL) delete C;
	if (G != NULL) delete G;
	C = NULL;
	G = NULL;
	
	for (unsigned int i = 0;i < burst;i++){
	  if (cache[i] != NULL) delete cache[i];
	}
	delete[] cache;
	fromCache = false;
}

void CurveGenerator::reset()
{
	curveCounter = 0;
	
	for (unsigned int i = 0;i < burst;i++){
	  if (cache[i] != NULL) delete cache[i];
	  cache[i] = NULL;
	}
	fromCache = false;
}

void CurveGenerator::revert()
{
	R.set(Q);
	S = origS;
	curveCounter = 0;
	Generator::revert();
}

int CurveGenerator::getCurveNumber()
{
	return curveCounter;
}

bool CurveGenerator::next(ReducedPoint& P,const mpz_t zN)
{
	if (G == NULL || C == NULL) return false;
	
	if (curveCounter == burst)
	{
		 fromCache = true;
		 return false;
	}
	
	if (!fromCache) 
	{
	   cache[curveCounter] = new RationalPoint();
	   
	   C->addPoints(R,R,*G); 
	   generate_base_point(*cache[curveCounter]);
	}
	
	P.reduce(*cache[curveCounter],zN);
	
	curveCounter++;
	S++;
	return true;
}
