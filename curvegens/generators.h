#ifndef GENERATORS_H
#define GENERATORS_H

#include <fstream>

#include "elliptic.h"
#include "../helpers.h"

typedef enum { Z12,Z2xZ8,Z6,Z8,Z2xZ4 } Torsion;

// Generický generátor křivek z libovolného zdroje
class Generator {
private:
	unsigned int edwards,twisted;

protected:
	unsigned int s;
	int A;
	
	mpz_t N;
	
	virtual bool next(ReducedPoint& P,mpz_t zN) = 0;
	virtual void reset() = 0;

public:

	Generator();
	
	int getCoeff();
	int getA();
	int countEdwards();
	int countTwisted();
	void getN(mpz_t r);
	
	void restart();
	bool next_base_point(ReducedPoint& P,mpz_t zN);
};
	
// Generátor křivek z nekonečných rodin
class CurveGenerator: public Generator {
protected:
	EllipticCurve *C;
	ReducedPoint* G,Q;
	Torsion T;
	unsigned int start,burst;

	bool next(ReducedPoint& P,mpz_t zN);
	void reset();
	
	virtual void generate_base_point(ReducedPoint& P,mpz_t zN) = 0;
public:
	CurveGenerator(Torsion t,unsigned int from,unsigned int b);
	~CurveGenerator();

};
	
///////////////////////// KŘIVKOVÉ GENERÁTORY //////////////////////////
class EdwardsGenerator : public CurveGenerator {

public:
	EdwardsGenerator(mpz_t n,Torsion T,unsigned int from,unsigned int b);
	
protected:
	void generate_base_point(ReducedPoint& P,mpz_t zN);
};

///////////////////////////// JINÉ GENERÁTORY //////////////////////////

// Generátor křivek skrz čtení ze souboru
class FileGenerator : public Generator {

private:
	ifstream fp;
	
public:
	FileGenerator(string filename);
	~FileGenerator();
	
protected:
	bool next(ReducedPoint& P,mpz_t zN);
	void reset();
};

// Generátor smíšených druhů křivek
class MixedGenerator : public Generator {
private:
	CurveGenerator** gens;
	unsigned int start,end;
	int ctr;

public:
	MixedGenerator(unsigned int from,unsigned int b);
	~MixedGenerator();
	
protected:
	bool next(ReducedPoint& P,mpz_t zN);
	void reset();
	
};

#endif
