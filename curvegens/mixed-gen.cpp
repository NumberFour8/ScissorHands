#include "generators.h"

MixedGenerator::MixedGenerator(unsigned int start,vector<GeneratorSetup> s)
	: Generator(b), ctr(0), setup(s)
{
	
	burst = std::accumulate(s.begin(),s.end(),[](const pair<GeneratorSetup> p){ return p.first; }));
	
	gens = new CurveGenerator*[setup.size()];
	for (unsigned int i = 0; i < setup.size();i++){
	  gens[i] = new EdwardsGenerator(setup[i].second,start,setup[i].first+1);
	  origSetup.push_back(setup[i].first);
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
	  setup[i].first = origSetup[i];
	}
}

void MixedGenerator::revert()
{
	for (unsigned int i = 0; i < setup.size();i++){
	  gens[i]->revert();
	  setup[i].first = origSetup[i];
	}
	Generator::revert();
}

bool MixedGenerator::next(RationalPoint& P)
{
	for (unsigned int i = 0;i < setup.size();i++)
	{
	   if (setup[i].first > 0)
	   {
			bool r =  gens[i]->next(P); 
			A = gens[i]->getA();
		    setup[i].first -= 1;
			return r;
	   }
	} 
	return false;
}
