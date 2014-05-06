#ifndef GENERATORS_H
#define GENERATORS_H

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#include "elliptic.h"

class CurveGenerator {
private:
	unsigned int edwards,twisted;

protected:
	unsigned int s;
	int A;
	virtual bool next(ReducedPoint& P) = 0;
	mpz_t N;

public: 
	BasicGenerator(mpz_t n) : s(1), edwards(0), twisted(0), A(0) 
	{ 
		mpz_init_set(N,n);
	};
	
	~BasicGenerator()
	{
		mpz_clear(N);
	}
	
	int getCoeff() { return s; }
	mpz_t getN() { return N; }
	
	int countEdwards() { return edwards; }
	int countTwisted() { return twisted; }
	int getCurrentA()  { return A; }
	
	bool next_base_point(ReducedPoint& P)
	{
		if (next(P))
		{
			if (A == 1){
			  edwards++;
			  return true;
			}
			else if (A == -1)
			{
			  twisted++;
			  return true;
			}
		}
			
		return false;
	}
	
};
	
class FileGenerator : public CurveGenerator {

private:
	ifstream fp;
	
public:
	FileGenerator(mpz_t n,string filename);
	~FileGenerator();
	
protected:
	bool next(ReducedPoint& P);
};

class EdwardsGenerator : public CurveGenerator {

private:
	EllipticCurve C;
	ReducedPoint G;
	
public:
	enum EdwardsTorsion { C12,C2x8 };

	EdwardsGenerator(mpz_t n,EdwardsTorsion T);
	~EdwardsGenerator();
	
protected:
	bool next(ReducedPoint& P);
};

class TwistedGenerator : public CurveGenerator {

private:
	EllipticCurve C;
	ReducedPoint G;
	
public:
	enum TwistedTorsion { C6,C2x4,C8 };

	TwistedGenerator(mpz_t n,TwistedTorsion T);
	~TwistedGenerator();
	
protected:
	bool next(ReducedPoint& P);
};

#endif
