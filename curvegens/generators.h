#ifndef GENERATORS_H
#define GENERATORS_H

#include <fstream>

#include "elliptic.h"
#include "../helpers.h"

typedef enum { Z12,Z2xZ8,Z6,Z8,Z2xZ4 } Torsion;

// Generick� gener�tor k�ivek
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
	
// Gener�tor k�ivek z nekon�en�ch rodin
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
	
///////////////////////// K�IVKOV� GENER��TORY //////////////////////////
class EdwardsGenerator : public CurveGenerator {

public:
	EdwardsGenerator(mpz_t n,Torsion T,unsigned int from,unsigned int b);
	
protected:
	void generate_base_point(ReducedPoint& P,mpz_t zN);
};

///////////////////////////// JIN� GENER��TORY //////////////////////////

// Gener�tor k�ivek ze souboru
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

// Sm�en� gener�tor k�ivek
class MixedGenerator : public Generator {
private:
	CurveGenerator** gens;
	unsigned int start,end;
	int ctr;

public:
	MixedGenerator(mpz_t n,unsigned int from,unsigned int b);
	~MixedGenerator();
	
protected:
	bool next(ReducedPoint& P,mpz_t zN);
	void reset();
	
};

#endif
