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
	try {
		if (next(P,zN) && A*A == 1)
		{
			if (A == 1) edwards++;
			else 		twisted++;
			
			return true;
		}
		else return false;
	}	
	catch (mpz_t f)
	{
		cout << "ERROR: Cannot reduce on curve #" << (S+1) << endl;
		if (mpz_cmp_ui(f,0) != 0) // Byl nalezen faktor?
		{
			cout << "Factor found: " << mpz_to_string(f) << endl;
		}
		return false;
	}
}

CurveGenerator::CurveGenerator(Torsion t,unsigned int from,unsigned int b)
	: Generator(), T(t), start(from), burst(b+2), C(NULL), G(NULL)
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
}

bool CurveGenerator::next(ReducedPoint& P,const mpz_t zN)
{
	if (G == NULL || C == NULL) return false;
	
	if (S%burst == 0) return false;
	
	if (S <= 1)
	{
	  C->doublePoint(Q,*G);
	  S = 2;
	}
	else {
		for (;S < start-1;S++)
		   C->addPoints(Q,Q,*G);
		C->addPoints(Q,Q,*G);
	}
	
	generate_base_point(P,zN);
	S++;
	return true;
}
