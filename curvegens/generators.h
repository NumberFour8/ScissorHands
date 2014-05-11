#ifndef GENERATORS_H
#define GENERATORS_H

#include <fstream>

#include "elliptic.h"
#include "../helpers.h"

typedef enum { Z12,Z2xZ8,Z6,Z8,Z2xZ4 } Torsion;

// Generickı generátor køivek
class Generator {
private:
	unsigned int edwards,twisted;

protected:
	unsigned int S;
	int A;
		
	virtual bool next(ReducedPoint& P,const mpz_t zN) = 0;
	virtual void reset() = 0;

public:
	Generator();
	
	int getCoeff();
	int getA();
	int countEdwards();
	int countTwisted();
	int countCurves();
	void getN(mpz_t r);
	
	void new_point_set();
	bool next_base_point(ReducedPoint& P,const mpz_t zN);

	virtual void revert() = 0;
};
	
// Generátor køivek z nekonèenıch rodin
class CurveGenerator: public Generator {
private:
	unsigned int curveCounter;
protected:
	EllipticCurve *C;
	RationalPoint* G,Q,R; // G je generator, Q je startovni a R pracovni  
	Torsion T;
	unsigned int burst;

	bool next(ReducedPoint& P,const mpz_t zN);
	void reset();
	void initialize(unsigned int from);

	virtual void generate_base_point(ReducedPoint& P,const mpz_t zN) = 0;
public:
	CurveGenerator(Torsion t,unsigned int b);
	~CurveGenerator();
	
	int getCurveNumber();
	void revert();
};
	
///////////////////////// KØIVKOVÉ GENERÁTORY //////////////////////////
class EdwardsGenerator : public CurveGenerator {

public:
	EdwardsGenerator(Torsion T,unsigned int from,unsigned int b);
	
protected:
	void generate_base_point(ReducedPoint& P,const mpz_t zN);
};

///////////////////////////// JINÉ GENERÁTORY //////////////////////////

// Generátor køivek ze souboru
class FileGenerator : public Generator {

private:
	ifstream fp;
	
public:
	FileGenerator(string filename);
	~FileGenerator();
	
	void revert();

protected:
	bool next(ReducedPoint& P,const mpz_t zN);
	void reset();
};

// Smíšenı generátor køivek
class MixedGenerator : public Generator {
private:
	CurveGenerator** gens;
	unsigned int burst,ctr;
	const int num_gens;
	
public:
	MixedGenerator(unsigned int start,unsigned int b);
	~MixedGenerator();
	
	void revert();

protected:
	bool next(ReducedPoint& P,const mpz_t zN);
	void reset();
	
};

#endif
