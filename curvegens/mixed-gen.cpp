#include "generators.h"

MixedGenerator::MixedGenerator(mpz_t n,unsigned int start,unsigned int b)
	: Generator(), ctr(0), burst(b)
{
	gens = new CurveGenerator*[5];
	
	gens[0] = new EdwardsGenerator(n,Z6,start,burst);
	gens[1] = new EdwardsGenerator(n,Z12,start,burst);
	gens[2] = new EdwardsGenerator(n,Z8,start,burst);
	gens[3] = new EdwardsGenerator(n,Z2xZ8,start,burst);
	gens[4] = new EdwardsGenerator(n,Z2xZ4,start,burst);
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
	gens[0]->restart();
	gens[1]->restart();
	gens[2]->restart();
	gens[3]->restart();
	gens[4]->restart();
	ctr = 0;
}

bool MixedGenerator::next(ReducedPoint& P,mpz_t zN)
{
	if (ctr == burst) return false; 

	bool r =  gens[ctr % 5]->next_base_point(P,zN); 
	ctr++;
	return r; 
}
