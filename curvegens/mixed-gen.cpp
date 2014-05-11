#include "generators.h"

MixedGenerator::MixedGenerator(unsigned int start,unsigned int b)
	: Generator(), ctr(0), burst(b), num_gens(4)
{
	gens = new CurveGenerator*[num_gens];
	
	gens[0] = new EdwardsGenerator(Z6,start,burst);
	gens[1] = new EdwardsGenerator(Z12,start,burst);
	gens[2] = new EdwardsGenerator(Z8,start,burst);
	gens[3] = new EdwardsGenerator(Z2xZ8,start,burst);
	//gens[4] = new EdwardsGenerator(n,Z2xZ4,start,burst);
}

MixedGenerator::~MixedGenerator()
{
	delete gens[0];
	delete gens[1];
	delete gens[2];
	delete gens[3];
	//delete gens[4];
	delete[] gens;
}

void MixedGenerator::reset()
{
	gens[0]->restart();
	gens[1]->restart();
	gens[2]->restart();
	gens[3]->restart();
	//gens[4]->restart();
	ctr = 0;
}

bool MixedGenerator::next(ReducedPoint& P,const mpz_t zN)
{
	if (ctr == burst) return false; 

	bool r =  gens[ctr % num_gens]->next_base_point(P,zN); 
	A = gens[ctr % num_gens]->getA();
	ctr++;
	return r; 
}
