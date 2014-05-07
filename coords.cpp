#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include "coords.h"

void ExtendedPoint::initAll(bool minusOne)
{
	reset(X);
	reset(Y);
	reset(Z);
	reset(T);
	isMinus1 = minusOne;
}

ExtendedPoint::ExtendedPoint(bool minusOne)
{
	initAll(minusOne);
}
	
ExtendedPoint::ExtendedPoint(mpz_t x,mpz_t y,mpz_t N,bool minusOne) 
{
	initAll(minusOne);
	fromAffine(x,y,N);
}

ExtendedPoint::ExtendedPoint(mpz_t N,bool minusOne) 
{
	initAll(minusOne);
	infinity(N);
}

void ExtendedPoint::infinity(mpz_t N)
{
	mpz_t x,y;
	mpz_init_set_ui(x,0);
	mpz_init_set_ui(y,1);

	fromAffine(x,y,N);
	mpz_clrs(x,y);
}

void ExtendedPoint::fromAffine(mpz_t x,mpz_t y,mpz_t N)
{
	mpz_t z,t;
	mpz_init_set_ui(z,1);
	mpz_init(t);
		
	mpz_mul(t,x,y);
		
	to_mont_repr(x,N);
	to_mont_repr(y,N);
	to_mont_repr(z,N);
	to_mont_repr(t,N);

	mpz_to_biguint(X,x);
	mpz_to_biguint(Y,y);
	mpz_to_biguint(Z,z);
	mpz_to_biguint(T,t);

	mpz_clrs(t,z);
}

bool ExtendedPoint::toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,mpz_t fact)
{
	mpz_t z;
	mpz_init(z);
	mpz_set_ui(fact,0);

	biguint_to_mpz(x,X);
	biguint_to_mpz(y,Y);
	biguint_to_mpz(z,Z);

	from_mont_repr(x,N,invB);
	from_mont_repr(y,N,invB);
	from_mont_repr(z,N,invB);
	
	bool ret = try_invert_mod(fact,z,N);
	if (ret) 
	{
		mpz_mul(x,x,fact);
		mpz_mul(y,y,fact);
		mpz_mod(x,x,N);
		mpz_mod(y,y,N);
		
		mpz_set_ui(fact,0);
	}
		
	mpz_clear(z);
	return ret;
}

// Zvoli vhodnou strategii vypoctu podle poctu nactenych typu krivek 
computeStrategy chooseStrategy(int edwardsRead,int twistedRead,int& usableCurves)
{

	int curvesRead 	    = edwardsRead+twistedRead;
	if (curvesRead <= 0)
	{
		cout << "ERROR: No curves read." << endl;
		usableCurves = 0;
		return computeStrategy::csNone;
	}
	
	if (edwardsRead > 0 && twistedRead < edwardsRead && edwardsRead%CURVES_PER_BLOCK == 0)
	{
		cout << "INFO: Using " << edwardsRead << " Edwards curves." << endl;
		usableCurves = edwardsRead;
		return computeStrategy::csEdwards;
	}
	
	if (twistedRead > 0 && edwardsRead < twistedRead && twistedRead%CURVES_PER_BLOCK == 0)
	{
		cout << "INFO: Using " << twistedRead << " twisted Edwards curves." << endl;
		usableCurves = twistedRead;
		return computeStrategy::csTwisted;
	}
	
	if (twistedRead*edwardsRead > 0 && twistedRead == edwardsRead && curvesRead%(2*CURVES_PER_BLOCK) == 0)
	{
		cout << "INFO: Using " << edwardsRead << " Edwards curves and " << twistedRead << " twisted Edwards curves." << endl;
		usableCurves = curvesRead;
		return computeStrategy::csMixed; 
	}
	
	cout << "ERROR: Inappropriate number of curves: " << edwardsRead << " Edwards curves, " << twistedRead << " twisted Edwards curves." << endl;
	usableCurves = 0;
	return computeStrategy::csNone;
}

computeStrategy readCurves(Generator* source,ExtendedPoint** pInit,int& edwards,int& twisted,int &usableCurves)
{
	computeStrategy strat = computeStrategy::csNone;
	vector<ExtendedPoint> v;
	
	cout << "Loading curves..." << endl;
	
	ReducedPoint P;
	while (source->next_base_point(P))
	{
		// Vytvor bod v Extended souradnicích z redukovanych afinnich bodu modulo N
		v.push_back(ExtendedPoint(P.X,P.Y,source->N,source->getA() == -1)); 
	}

	twisted = source->countTwisted();
	edwards = source->countEdwards();
	cout << "Curve generation finished." << endl;

	// Prekroucene Edwardsovy krivky prijdou na zacatek
    std::sort(v.begin(), v.end(), [](const ExtendedPoint& a, const ExtendedPoint & b) -> bool { return a.isMinus1 && !b.isMinus1; });
	
	usableCurves = 0;
	strat = chooseStrategy(edwards,twisted,usableCurves);
	if (strat == computeStrategy::csNone) goto read_finish;

	// Prekopiruj body z vektoru do pameti
	*pInit = new ExtendedPoint[usableCurves];
	std::copy((strat == computeStrategy::csEdwards ? v.begin()+twisted : v.begin()),v.begin()+usableCurves,*pInit);

	// Jsme hotovi
	read_finish:
	v.clear();
	return strat;
}
