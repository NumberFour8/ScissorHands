#include "generators.h"

MixedGenerator::MixedGenerator(unsigned int start,unsigned int b,vector<GeneratorSetup> s)
	: Generator(b), setup(s)
{
	gens = new CurveGenerator*[setup.size()];
	for (unsigned int i = 0; i < setup.size();i++){
	  gens[i] = new EdwardsGenerator(setup[i].second,start,setup[i].first);
	}
	
}

MixedGenerator::~MixedGenerator()
{
	for (unsigned int i = 0;i < setup.size();i++)
	  delete gens[i];
	
	delete[] gens;
}

void MixedGenerator::reset()
{
	for (unsigned int i = 0; i < setup.size();i++){
	  gens[i]->reset();
	}
}

void MixedGenerator::revert()
{
	for (unsigned int i = 0; i < setup.size();i++){
	  gens[i]->revert();
	}
	Generator::revert();
}

bool MixedGenerator::next(RationalPoint& P)
{
	for (unsigned int i = 0;i < setup.size();i++)
	{
		if (gens[i]->next_base_point(P))
		{
		   return true;
		}
	} 
	return false;
}
