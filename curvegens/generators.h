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
	virtual void reset() = 0;

public:
	mpz_t N;

	CurveGenerator(mpz_t n);
	~CurveGenerator();
	
	int getCoeff();
	int getA();
	int countEdwards();
	int countTwisted();
	
	void restart();
	bool next_base_point(ReducedPoint& P);
};
	
class FileGenerator : public CurveGenerator {

private:
	ifstream fp;
	
public:
	FileGenerator(mpz_t n,string filename);
	~FileGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
};

enum Torsion { Z12 = 1,Z2xZ8 = 2,Z6 = 3,Z8 = 4,Z2xZ4 = 5 };

class EdwardsGenerator : public CurveGenerator {

private:
	EllipticCurve* C;
	ReducedPoint* G;
	Torsion tor;
	
public:

	EdwardsGenerator(mpz_t n,Torsion T);
	~EdwardsGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
};

class TwistedGenerator : public CurveGenerator {

private:
	EllipticCurve* C;
	ReducedPoint* G;
	Torsion tor;
	
public:
	
	TwistedGenerator(mpz_t n,Torsion T);
	~TwistedGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
};

class MixedGenerator : public CurveGenerator {
private:
	CurveGenerator** gens;

public:
	MixedGenerator(mpz_t n);
	~MixedGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
	
};

#endif
