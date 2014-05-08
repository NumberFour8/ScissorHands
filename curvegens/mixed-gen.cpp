#include "generators.h"

MixedGenerator::MixedGenerator(mpz_t n,unsigned int start,unsigned int b)
	: Generator(), ctr(-1)
{
	gens = new CurveGenerators*[5];
	unsigned int bb = b/5;
	
	gens[0] = new EdwardsGenerator(n,Z6,start,bb);
	gens[1] = new EdwardsGenerator(n,Z12,start,bb);
	gens[2] = new EdwardsGenerator(n,Z8,start,bb);
	gens[3] = new EdwardsGenerator(n,Z2xZ8,start,bb);
	gens[4] = new EdwardsGenerator(n,Z2xZ4,start,bb);
}

MixedGenerator::~MixedGenerator()
{
	delete gens[0];
	delete gens[1];
	delete gens[2];
	delete gens[3];
	delete gens[4];
	delete[] gens;
}

void MixedGenerator::reset()
{
	gens[0]->reset();
	gens[1]->reset();
	gens[2]->reset();
	gens[3]->reset();
	gens[4]->reset();
}

bool MixedGenerator::next(ReducedPoint& P,mpz_t zN)
{
	ctr = (ctr+1) % 5;
	return gens[ctr]->next(P,zN); 
}
