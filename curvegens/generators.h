#ifndef GENERATORS_H
#define GENERATORS_H

#include <fstream>

#include "elliptic.h"
#include "../helpers.h"

typedef enum { Z12,Z2xZ8,Z6,Z8,Z2xZ4 } Torsion;

// Generickı generátor køivek
class Generator {
private:
	unsigned int edwards,twisted,curveCounter;
	bool fromCache;
	RationalPoint** cache;
	
protected:
	const unsigned int burst;
		
	unsigned int S;
	int A;
		
	virtual bool next(RationalPoint& P) = 0;
	virtual void reset() = 0;

public:
	Generator(unsigned int b);
	~Generator();
	
	int getCoeff();
	int getA();
	int countEdwards();
	int countTwisted();
	int countCurves();
	void getN(mpz_t r);
	
	void new_point_set();
	bool next_base_point(ReducedPoint& P,const mpz_t zN);

	virtual void revert();
};
	
// Generátor køivek z nekonèenıch rodin
class CurveGenerator: public Generator {
private:
	unsigned int origS;
	
protected:
	EllipticCurve *C;
	RationalPoint* G,Q,R; // G je generator, Q je startovni a R pracovni  
	Torsion T;

	bool next(RationalPoint& P);
	void reset();
	void initialize(unsigned int from);

	virtual void generate_base_point(RationalPoint& P) = 0;
public:
	CurveGenerator(Torsion t,unsigned int b);
	~CurveGenerator();
	
	int getCurveNumber();
	void revert();
	
	friend class MixedGenerator;
};
	
///////////////////////// KØIVKOVÉ GENERÁTORY //////////////////////////
class EdwardsGenerator : public CurveGenerator {

public:
	EdwardsGenerator(Torsion T,unsigned int from,unsigned int b);
	
protected:
	void generate_base_point(RationalPoint& P);
};

///////////////////////////// JINÉ GENERÁTORY //////////////////////////

// Generátor køivek ze souboru
class FileGenerator : public Generator {

private:
	ifstream fp;
	
public:
	FileGenerator(string filename,unsigned int cacheSize);
	~FileGenerator();
	
	void revert();

protected:
	bool next(RationalPoint& P);
	void reset();
};

// Smíšenı generátor køivek
class MixedGenerator : public Generator {
private:
	CurveGenerator** gens;
	unsigned int ctr;
	const int num_gens;
	
public:
	MixedGenerator(unsigned int start,unsigned int b);
	~MixedGenerator();
	
	void revert();

protected:
	bool next(RationalPoint& P);
	void reset();
	
};

#endif
