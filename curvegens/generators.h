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
	
	virtual bool next(ReducedPoint& P) = 0;
	virtual void reset() = 0;

public:
	mpz_t N;

	Generator(mpz_t n);
	~Generator();
	
	int getCoeff();
	int getA();
	int countEdwards();
	int countTwisted();
	
	void restart();
	bool next_base_point(ReducedPoint& P);
};
	
// Generátor křivek z nekonečných rodin
class CurveGenerator: public Generator {
protected:
	EllipticCurve *C;
	ReducedPoint* G,Q;
	Torsion T;
	unsigned int start,burst;

	bool next(ReducedPoint& P);
	void reset();
	
	virtual void generate_base_point(ReducedPoint& P) = 0;
public:
	CurveGenerator(mpz_t n,Torsion t,unsigned int from,unsigned int b);
	~CurveGenerator();

};
	
///////////////////////// KŘIVKOVÉ GENERÁTORY //////////////////////////
class EdwardsGenerator : public CurveGenerator {

public:
	EdwardsGenerator(mpz_t n,Torsion T,unsigned int from,unsigned int b);
	
protected:
	void generate_base_point(ReducedPoint& P);
};

///////////////////////////// JINÉ GENERÁTORY //////////////////////////

// Generátor křivek skrz čtení ze souboru
class FileGenerator : public Generator {

private:
	ifstream fp;
	
public:
	FileGenerator(mpz_t n,string filename);
	~FileGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
};

// Generátor smíšených druhů křivek
class MixedGenerator : public Generator {
private:
	CurveGenerator** gens;
	unsigned int start,end;

public:
	MixedGenerator(mpz_t n,unsigned int from,unsigned int b);
	~MixedGenerator();
	
protected:
	bool next(ReducedPoint& P);
	void reset();
	
};

#endif
