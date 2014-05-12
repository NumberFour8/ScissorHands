#include "generators.h"

MixedGenerator::MixedGenerator(unsigned int start,unsigned int b)
	: Generator(b), ctr(0), num_gens(4)
{
	gens = new CurveGenerator*[num_gens];
	
	gens[0] = new EdwardsGenerator(Z6,start,burst);
	gens[1] = new EdwardsGenerator(Z12,start,burst);
	gens[2] = new EdwardsGenerator(Z8,start,burst);
	gens[3] = new EdwardsGenerator(Z2xZ8,start,burst);
	//gens[4] = new EdwardsGenerator(Z2xZ4,start,burst);
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
	gens[0]->new_point_set();
	gens[1]->new_point_set();
	gens[2]->new_point_set();
	gens[3]->new_point_set();
	//gens[4]->new_point_set();
	ctr = 0;
}

void MixedGenerator::revert()
{
	gens[0]->revert();
	gens[1]->revert();
	gens[2]->revert();
	gens[3]->revert();
	Generator::revert();
}

bool MixedGenerator::next(RationalPoint& P)
{
	if (ctr == burst) return false; 

	bool r =  gens[ctr % num_gens]->next(P); 
	A = gens[ctr % num_gens]->getA();
	ctr++;
	return r; 
}
