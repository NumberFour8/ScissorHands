#include "generators.h"

Generator::Generator(unsigned int b) 
 : S(1), edwards(0), twisted(0), A(0), curveCounter(0),burst(b) 
{ 
	cache = new RationalPoint*[burst];
	for (unsigned int i = 0;i < burst;i++)
		cache[i] = NULL;

	fromCache = false;
}

Generator::~Generator()
{
	for (unsigned int i = 0;i < burst;i++){
	  if (cache[i] != NULL) delete cache[i];
	}
	delete[] cache;
	fromCache = false;
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
	return curveCounter;
}

void Generator::new_point_set()
{
	edwards = twisted = 0;
	curveCounter = 0;
	
	for (unsigned int i = 0;i < burst;i++){
	  if (cache[i] != NULL) delete cache[i];
	  cache[i] = NULL;
	}
	fromCache = false;
	
	reset();
}

void Generator::revert()
{
	edwards = twisted = 0;
	curveCounter = 0;
}

bool Generator::next_base_point(ReducedPoint& P,const mpz_t zN)
{	
	if (curveCounter == burst)
	{
		 fromCache = true;
		 return false;
	}
	
	if (!fromCache) 
	{
	   cache[curveCounter] = new RationalPoint();
	   next(*(cache[curveCounter]));
	}
	P.reduce(*cache[curveCounter],zN);
	
	curveCounter++;
		
	if (A == 1) edwards++;
	else 		twisted++;
		
	return true;
}

CurveGenerator::CurveGenerator(Torsion t,unsigned int b)
	: Generator(b), T(t), C(NULL), G(NULL)
{ 
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
}

void CurveGenerator::reset()
{	
	return;
}

void CurveGenerator::revert()
{
	R.set(Q);
	S = origS;
	Generator::revert();
}

bool CurveGenerator::next(RationalPoint& P)
{
	if (G == NULL || C == NULL) return false;
	
	C->addPoints(R,R,*G);
	generate_base_point(P);
	S++;
	return true;
}
