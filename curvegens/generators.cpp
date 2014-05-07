#include "generators.h"

Generator::Generator(mpz_t n) 
 : s(1), edwards(0), twisted(0), A(0) 
{ 
	mpz_init_set(N,n);
}

Generator::~Generator()
{
	mpz_clear(N);
}

int Generator::getCoeff() 
{
	return s;
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

void Generator::restart()
{
	A = s = edwards = twisted = 0;
	reset();
}

bool Generator::next_base_point(ReducedPoint& P)
{
	if (next(P) && A*A == 1)
	{
		if (A == 1) edwards++;
		else 		twisted++;
		
		return true;
	}
	return false;
}

CurveGenerator::CurveGenerator(mpz_t n,Torsion t,unsigned int from,unsigned int b)
	: Generator(n), T(t), start(from), burst(b+2), C(NULL), G(NULL)
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
	s = 1;
}

bool CurveGenerator::next(ReducedPoint& P)
{
	if (G == NULL || C == NULL) return false;
	
	if (s == burst) return false;
	
	if (s <= 1)
	{
	  C->doublePoint(Q,*G);
	  s = 2;
	}
	else {
		for (;s < start-1;s++)
		   C->addPoints(Q,Q,*G);
		C->addPoints(Q,Q,*G);
	}
	
	try
	{
	  generate_base_point(P);
	  ++s;
	  return true;
	}	
	catch (mpz_t f)
	{
		cout << "ERROR: Cannot reduce on curve #" << (s+1) << endl;
		if (mpz_cmp_ui(f,0) != 0) // Byl nalezen faktor?
		{
			cout << "Factor found: " << mpz_to_string(f) << endl;
		}
		return false;
	}
}
