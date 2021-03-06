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

bool Generator::next_base_point(RationalPoint& P)
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

	P.set(*(cache[curveCounter]));
	A = cache[curveCounter]->minus1 ? -1 : 1;
	curveCounter++;		

	if (A == 1) edwards++;
	else 		twisted++;
		
	return true;
}

bool Generator::next_base_point(ReducedPoint& P,const mpz_t n)
{
	RationalPoint pt; 
	bool r = next_base_point(pt);
	P.reduce(pt,n);
	return r;
}

CurveGenerator::CurveGenerator(Torsion t,unsigned int b)
	: Generator(b), T(t), C(NULL), G(NULL)
{ 
}

void CurveGenerator::initialize(unsigned int from)
{
	C->doublePoint(Q,*G);
	for (S = 2;S < from;S++)
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
	
	generate_base_point(P);
	C->addPoints(R,R,*G);

	S++;
	return true;
}
